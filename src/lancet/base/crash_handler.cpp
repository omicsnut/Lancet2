#include "lancet/base/crash_handler.h"

#include "lancet/base/types.h"

#include <ucontext.h>

// ============================================================================
// Crash handler implementation — POSIX signal-based crash diagnostics with
// per-thread window context and all-thread enumeration.
//
// On SIGSEGV or SIGABRT the handler prints to stderr:
//   1. Faulting address (si_addr) — null deref vs heap vs stack corruption
//   2. Signal code (si_code)     — SEGV_MAPERR (1) or SEGV_ACCERR (2)
//   3. Thread ID                 — identifies the crashing worker thread
//   4. Instruction pointer (RIP) — the exact crashing instruction
//   5. Stack backtrace           — best-effort, may be empty if stack is corrupt
//   6. All registered thread contexts — which window each worker was processing
//   7. /proc/self/task enumeration — instruction pointer of every OS thread
//
// The handler runs on a dedicated 64KB alternate signal stack (sigaltstack)
// so it can function even when the crash is caused by stack overflow.
// After printing, it re-raises the signal with SIG_DFL to produce a core dump.
//
// All output uses write() and manual hex formatting — no heap allocation.
// Pipe output through `c++filt` to demangle C++ symbols in the backtrace.
// Resolve the RIP with: addr2line -e ./Lancet2 -f <RIP address>
// ============================================================================

#if defined(__linux__) || defined(__APPLE__)

#include <execinfo.h>
#include <fcntl.h>
#include <pthread.h>
#include <unistd.h>

#ifdef __linux__
#include <sys/syscall.h>
#endif

#include <array>

#include <csignal>
#include <cstring>

namespace {

// Maximum number of stack frames to capture in the backtrace.
constexpr int MAX_BACKTRACE_FRAMES = 64;

// Minimum alternate stack size — at least 64KB or SIGSTKSZ, whichever is larger.
constexpr usize ALT_STACK_SIZE = SIGSTKSZ > 65'536 ? SIGSTKSZ : 65'536;

// ============================================================================
// Async-signal-safe output helpers — no heap allocation, no stdio.
// ============================================================================

// Async-signal-safe write of a NUL-terminated string to a file descriptor.
void WriteStr(int file_desc, char const* str) {
  usize len = 0;
  while (str[len] != '\0') {
    len++;
  }
  // NOLINTNEXTLINE(bugprone-signal-handler,cert-sig30-c) — write() is async-signal-safe per POSIX
  static_cast<void>(write(file_desc, str, len));
}

// Async-signal-safe hex formatter — writes "0x" + 16 hex digits, no heap allocation.
void WriteHex(int file_desc, u64 val) {
  // NOLINTNEXTLINE(modernize-avoid-c-arrays) — raw byte buffer for signal-safe formatting
  char buf[18] = {'0', 'x'};
  for (int idx = 17; idx >= 2; --idx) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
    buf[idx] = "0123456789abcdef"[val & 0xFU];
    val >>= 4U;
  }
  // NOLINTNEXTLINE(bugprone-signal-handler,cert-sig30-c) — write() is async-signal-safe per POSIX
  static_cast<void>(write(file_desc, buf, sizeof(buf)));
}

// Async-signal-safe decimal integer writer for small values.
void WriteInt(int file_desc, int val) {
  // NOLINTNEXTLINE(modernize-avoid-c-arrays) — raw byte buffer for signal-safe formatting
  char buf[16];
  int pos = 0;
  if (val < 0) {
    buf[pos++] = '-';
    val = -val;
  }

  // Format digits in reverse, then reverse the digit portion.
  int const digit_start = pos;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-do-while)
  do {
    buf[pos++] = static_cast<char>('0' + (val % 10));
    val /= 10;
  } while (val > 0);

  // Reverse the digits in place.
  for (int lo = digit_start, hi = pos - 1; lo < hi; ++lo, --hi) {
    char const tmp = buf[lo];
    buf[lo] = buf[hi];
    buf[hi] = tmp;
  }

  // NOLINTNEXTLINE(bugprone-signal-handler,cert-sig30-c) — write() is async-signal-safe per POSIX
  static_cast<void>(write(file_desc, buf, pos));
}

// Async-signal-safe decimal writer for u64 values (used for genome index).
void WriteU64(int file_desc, u64 val) {
  // NOLINTNEXTLINE(modernize-avoid-c-arrays) — raw byte buffer for signal-safe formatting
  char buf[20];
  int pos = 0;

  // NOLINTNEXTLINE(cppcoreguidelines-avoid-do-while)
  do {
    buf[pos++] = static_cast<char>('0' + (val % 10));
    val /= 10;
  } while (val > 0);

  // Reverse in place.
  for (int lo = 0, hi = pos - 1; lo < hi; ++lo, --hi) {
    char const tmp = buf[lo];
    buf[lo] = buf[hi];
    buf[hi] = tmp;
  }

  // NOLINTNEXTLINE(bugprone-signal-handler,cert-sig30-c) — write() is async-signal-safe per POSIX
  static_cast<void>(write(file_desc, buf, pos));
}

// ============================================================================
// Per-thread crash context slot array.
//
// Each worker thread claims a slot via atomic compare-exchange on mInUse.
// The crash handler reads all slots — no locks needed since writes are
// atomic stores and reads tolerate tearing (worst case: stale region string).
//
// Memory layout per slot (cache-line aligned to prevent false sharing):
//   mInUse       — 0 = free, 1 = registered thread, read atomically
//   mActive      — 0 = idle, 1 = currently processing a window
//   mGenomeIdx   — Window::GenomeIndex() of the window being processed
//   mThreadId    — pthread_t of the owning thread (for identification)
//   mRegion      — human-readable region string ("chr4:12345-67890")
// ============================================================================

struct alignas(64) CrashSlot {
  int volatile mInUse;      // 4B — atomic flag: 0 = free, 1 = claimed by a thread
  int volatile mActive;     // 4B — atomic flag: 0 = idle, 1 = processing a window
  u64 volatile mGenomeIdx;  // 8B — genome-wide window ordinal being processed
  pthread_t mThreadId;      // 8B — owning thread, used to correlate with crashing thread
  // NOLINTNEXTLINE(modernize-avoid-c-arrays) — fixed buffer for async-signal-safe region string
  char mRegion[lancet::base::MAX_REGION_LEN];  // 128B
} /* 192B, padded to 192B by alignas(64) */;

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables) — mutable by design
std::array<CrashSlot, lancet::base::MAX_CRASH_SLOTS> g_crash_slots = {};

// ============================================================================
// DumpAllThreadContexts — called from the crash handler to print crash context.
//
// Iterates all crash slots and prints each registered thread's state:
//   - Thread ID (pthread_t as hex)
//   - Active/idle status
//   - If active: genome index + region string
// ============================================================================
void DumpAllThreadContexts(int file_desc) {
  WriteStr(file_desc, "\n── Worker Thread Contexts ──────────────────────\n");

  bool found_any = false;
  for (auto const& slot : g_crash_slots) {
    if (slot.mInUse == 0) continue;
    found_any = true;

    WriteStr(file_desc, "  Thread ");
    WriteHex(file_desc, static_cast<u64>(slot.mThreadId));

    if (slot.mActive != 0) {
      WriteStr(file_desc, " [ACTIVE] window_idx=");
      WriteU64(file_desc, slot.mGenomeIdx);
      WriteStr(file_desc, " region=");
      WriteStr(file_desc, slot.mRegion);
    } else {
      WriteStr(file_desc, " [idle]");
    }
    WriteStr(file_desc, "\n");
  }

  if (!found_any) {
    WriteStr(file_desc, "  (no worker threads registered)\n");
  }
}

#ifdef __linux__
// ============================================================================
// DumpAllThreadIPs — read instruction pointers for every OS thread via procfs.
//
// /proc/self/task/<tid>/syscall contains kernel-level register state:
//   Format: "syscall_nr arg0 arg1 ... SP IP"
//   The last two hex fields are stack pointer and instruction pointer.
//   If the thread is in userspace (not in a syscall), the line is "running".
//
// This gives us the IP of threads that did NOT register a crash slot
// (e.g., the main thread, profiler thread).
// ============================================================================

// Linux dirent64 layout for raw getdents64 parsing.  Must match kernel ABI
// exactly — field sizes and offsets are prescribed by the syscall interface.
// NOLINTNEXTLINE(readability-identifier-naming) — mirrors kernel struct name for clarity
struct linux_dirent64 {
  u64 mIno;                // 8B — inode number
  u64 mOff;                // 8B — offset to next dirent
  unsigned short mReclen;  // 2B — length of this record
  unsigned char mType;     // 1B — file type
  // NOLINTNEXTLINE(modernize-avoid-c-arrays) — variable-length kernel ABI field
  char mName[1];  // variable length — NUL-terminated name
};

// Dump the syscall register state for a single thread via /proc/self/task/<tid>/syscall.
void DumpOneThreadIP(int file_desc, char const* tid_str) {
  // Build /proc/self/task/<tid>/syscall path on the stack.
  // NOLINTNEXTLINE(modernize-avoid-c-arrays) — stack buffer for path assembly, no heap
  char path[64] = "/proc/self/task/";
  constexpr usize PREFIX_LEN = 16;
  usize plen = PREFIX_LEN;

  for (usize idx = 0; tid_str[idx] != '\0' && plen < 48; ++idx) {
    path[plen++] = tid_str[idx];
  }

  // Append "/syscall\0" — memcpy is async-signal-safe per POSIX.
  memcpy(path + plen, "/syscall", 9);

  // Read the syscall file — format: "nr a0 a1 a2 a3 a4 a5 SP IP"
  // NOLINTNEXTLINE(modernize-avoid-c-arrays) — stack buffer for procfs read, no heap
  char syscall_buf[256] = {};
  int const sfd = open(path, O_RDONLY);
  if (sfd < 0) return;
  auto const nread = read(sfd, syscall_buf, sizeof(syscall_buf) - 1);
  close(sfd);
  if (nread <= 0) return;

  WriteStr(file_desc, "  tid=");
  WriteStr(file_desc, tid_str);
  WriteStr(file_desc, " syscall: ");

  // Truncate to first newline for clean output.
  for (long idx = 0; idx < nread; ++idx) {
    if (syscall_buf[idx] == '\n') {
      syscall_buf[idx] = '\0';
      break;
    }
  }
  WriteStr(file_desc, syscall_buf);
  WriteStr(file_desc, "\n");
}

void DumpAllThreadIPs(int file_desc) {
  WriteStr(file_desc, "\n── All Thread IPs (via /proc/self/task) ────────\n");

  // Open /proc/self/task to enumerate thread IDs (tids).
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg) — open() is a POSIX variadic C function
  int const dir_fd = open("/proc/self/task", O_RDONLY | O_DIRECTORY);
  if (dir_fd < 0) {
    WriteStr(file_desc, "  (could not open /proc/self/task)\n");
    return;
  }

  // getdents64 buffer — stack-allocated, no heap.
  // NOLINTNEXTLINE(modernize-avoid-c-arrays) — raw byte buffer for syscall interop
  char dirents_buf[4096];

  while (true) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg) — syscall() is a POSIX variadic C function
    auto const nbytes = syscall(SYS_getdents64, dir_fd, dirents_buf, sizeof(dirents_buf));
    if (nbytes <= 0) break;

    for (long offset = 0; offset < nbytes;) {
      // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-bounds-pointer-arithmetic)
      auto const* entry = reinterpret_cast<linux_dirent64 const*>(dirents_buf + offset);
      offset += entry->mReclen;

      // Skip "." and ".." entries.
      if (entry->mName[0] == '.') continue;
      DumpOneThreadIP(file_desc, entry->mName);
    }
  }

  close(dir_fd);
}
#endif  // __linux__

// NOLINTNEXTLINE(cert-err33-c,cert-dcl03-c) — signal handler signature prescribed by POSIX
void CrashHandler(int sig, siginfo_t* info, void* ctx) {
  // SA_RESETHAND already reset us to SIG_DFL — the raise() at the end
  // produces a proper core dump via the default handler.

  int const file_desc = STDERR_FILENO;
  WriteStr(file_desc, "\n============================================\n");
  WriteStr(file_desc, sig == SIGSEGV ? "*** SIGSEGV ***\n" : "*** SIGABRT ***\n");

  // Faulting address — null deref (0x0), heap corruption, or stack overflow
  WriteStr(file_desc, "Faulting address (si_addr): ");
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast) — required to print void* as hex
  WriteHex(file_desc, reinterpret_cast<u64>(info->si_addr));
  WriteStr(file_desc, "\n");

  // Signal code — SEGV_MAPERR (1) = unmapped, SEGV_ACCERR (2) = no permission
  WriteStr(file_desc, "Signal code    (si_code):  ");
  WriteInt(file_desc, info->si_code);
  WriteStr(file_desc, "\n");

  // Thread ID — identifies which worker thread crashed
  WriteStr(file_desc, "Thread ID:                 ");
  WriteHex(file_desc, static_cast<u64>(pthread_self()));
  WriteStr(file_desc, "\n");

  // Instruction pointer from ucontext — the exact crashing instruction.
  // Resolve with: addr2line -e ./Lancet2 -f <RIP address>
#if defined(__linux__) && defined(__x86_64__)
  auto* uctx = static_cast<ucontext_t*>(ctx);
  auto const rip = static_cast<u64>(uctx->uc_mcontext.gregs[REG_RIP]);
  WriteStr(file_desc, "Instruction ptr (RIP):     ");
  WriteHex(file_desc, rip);
  WriteStr(file_desc, "\n");
#elif defined(__linux__) && defined(__aarch64__)
  auto* uctx = static_cast<ucontext_t*>(ctx);
  auto const prog_ctr = static_cast<u64>(uctx->uc_mcontext.pc);
  WriteStr(file_desc, "Instruction ptr (PC):      ");
  WriteHex(file_desc, prog_ctr);
  WriteStr(file_desc, "\n");
#elif defined(__APPLE__) && defined(__x86_64__)
  auto* uctx = static_cast<ucontext_t*>(ctx);
  auto const rip = static_cast<u64>(uctx->uc_mcontext->__ss.__rip);
  WriteStr(file_desc, "Instruction ptr (RIP):     ");
  WriteHex(file_desc, rip);
  WriteStr(file_desc, "\n");
#elif defined(__APPLE__) && defined(__aarch64__)
  auto* uctx = static_cast<ucontext_t*>(ctx);
  auto const prog_ctr = static_cast<u64>(uctx->uc_mcontext->__ss.__pc);
  WriteStr(file_desc, "Instruction ptr (PC):      ");
  WriteHex(file_desc, prog_ctr);
  WriteStr(file_desc, "\n");
#else
  static_cast<void>(ctx);
  WriteStr(file_desc, "Instruction ptr:           (unavailable on this arch)\n");
#endif

  // Best-effort backtrace — may produce 0 frames if the stack is too corrupted
  // for the DWARF unwinder to walk. The RIP above is always available.
  WriteStr(file_desc, "\nBacktrace:\n");
  std::array<void*, MAX_BACKTRACE_FRAMES> frames{};
  int const depth = backtrace(frames.data(), MAX_BACKTRACE_FRAMES);
  if (depth > 0) {
    backtrace_symbols_fd(frames.data(), depth, file_desc);
  } else {
    WriteStr(file_desc, "  (no frames — stack is corrupted)\n");
  }

  // Per-thread crash context — which window each worker was processing.
  DumpAllThreadContexts(file_desc);

  // All-thread instruction pointers via procfs (Linux only).
#ifdef __linux__
  DumpAllThreadIPs(file_desc);
#endif

  WriteStr(file_desc, "\n============================================\n");
  WriteStr(file_desc, "Resolve addresses with:\n");
  WriteStr(file_desc, "  addr2line -e ./Lancet2 -f <address>\n");
  WriteStr(file_desc, "============================================\n");

  // Re-raise with default handler to produce a core dump.
  // NOLINTNEXTLINE(cert-err33-c) — raise() return value is irrelevant during crash
  static_cast<void>(raise(sig));
}

}  // namespace

namespace lancet::base {

void InstallCrashHandler() {
  // Alternate signal stack (64KB) — the handler runs here even if the crash
  // is caused by stack overflow on the thread's normal stack.
  static std::array<char, ALT_STACK_SIZE> alt_stack_mem{};
  stack_t alt_stack{};
  alt_stack.ss_sp = static_cast<void*>(alt_stack_mem.data());
  alt_stack.ss_size = alt_stack_mem.size();
  alt_stack.ss_flags = 0;
  sigaltstack(&alt_stack, nullptr);

  // SA_SIGINFO  — gives us siginfo_t with faulting address and ucontext
  // SA_ONSTACK  — use the alternate signal stack
  // SA_RESETHAND — auto-reset to SIG_DFL before entering the handler
  struct sigaction act{};
  act.sa_sigaction = CrashHandler;
  act.sa_flags = SA_ONSTACK | SA_RESETHAND | SA_SIGINFO;
  sigemptyset(&act.sa_mask);
  sigaction(SIGSEGV, &act, nullptr);
  sigaction(SIGABRT, &act, nullptr);
}

auto RegisterThreadSlot() -> CrashSlotIdx {
  auto const tid = pthread_self();
  for (usize idx = 0; idx < MAX_CRASH_SLOTS; ++idx) {
    auto& slot = g_crash_slots[idx];
    int expected = 0;
    // Atomic compare-exchange: claim a free slot.
    if (__atomic_compare_exchange_n(&slot.mInUse, &expected, 1, false, __ATOMIC_SEQ_CST,
                                    __ATOMIC_SEQ_CST)) {
      slot.mThreadId = tid;
      __atomic_store_n(&slot.mActive, 0, __ATOMIC_RELEASE);
      slot.mGenomeIdx = 0;
      slot.mRegion[0] = '\0';
      return idx;
    }
  }
  return INVALID_CRASH_SLOT;
}

void UnregisterThreadSlot(CrashSlotIdx const slot_idx) {
  if (slot_idx >= MAX_CRASH_SLOTS) return;
  auto& slot = g_crash_slots[slot_idx];
  __atomic_store_n(&slot.mActive, 0, __ATOMIC_RELEASE);
  __atomic_store_n(&slot.mInUse, 0, __ATOMIC_RELEASE);
}

void SetSlotWindowInfo(CrashSlotIdx const slot_idx, u64 const genome_idx, char const* region_str) {
  if (slot_idx >= MAX_CRASH_SLOTS) return;
  auto& slot = g_crash_slots[slot_idx];

  // Write the genome index and region string before marking active.
  // The fence ensures the crash handler sees consistent data.
  slot.mGenomeIdx = genome_idx;

  // Copy region string with truncation — no heap allocation.
  usize pos = 0;
  if (region_str != nullptr) {
    for (; pos < MAX_REGION_LEN - 1 && region_str[pos] != '\0'; ++pos) {
      slot.mRegion[pos] = region_str[pos];
    }
  }
  slot.mRegion[pos] = '\0';

  // Mark active AFTER data is written — release fence.
  __atomic_store_n(&slot.mActive, 1, __ATOMIC_RELEASE);
}

void ClearSlotWindowInfo(CrashSlotIdx const slot_idx) {
  if (slot_idx >= MAX_CRASH_SLOTS) return;
  __atomic_store_n(&g_crash_slots[slot_idx].mActive, 0, __ATOMIC_RELEASE);
}

}  // namespace lancet::base

#else  // Non-POSIX platforms

namespace lancet::base {

void InstallCrashHandler() {
  // No-op on platforms without POSIX signal support.
}

auto RegisterThreadSlot() -> CrashSlotIdx {
  return INVALID_CRASH_SLOT;
}
void UnregisterThreadSlot(CrashSlotIdx /*slot_idx*/) {}
void SetSlotWindowInfo(CrashSlotIdx /*slot_idx*/, u64 /*genome_idx*/, char const* /*region_str*/) {}
void ClearSlotWindowInfo(CrashSlotIdx /*slot_idx*/) {}

}  // namespace lancet::base

#endif
