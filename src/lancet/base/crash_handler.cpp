#include "lancet/base/crash_handler.h"

#include "lancet/base/types.h"

#include <ucontext.h>

// ============================================================================
// Crash handler implementation — POSIX signal-based crash diagnostics.
//
// On SIGSEGV or SIGABRT the handler prints to stderr:
//   1. Faulting address (si_addr) — null deref vs heap vs stack corruption
//   2. Signal code (si_code)     — SEGV_MAPERR (1) or SEGV_ACCERR (2)
//   3. Thread ID                 — identifies the crashing worker thread
//   4. Instruction pointer (RIP) — the exact crashing instruction
//   5. Stack backtrace           — best-effort, may be empty if stack is corrupt
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
#include <pthread.h>
#include <unistd.h>

#include <csignal>

namespace {

// Maximum number of stack frames to capture in the backtrace.
constexpr int MAX_BACKTRACE_FRAMES = 64;

// Minimum alternate stack size — at least 64KB or SIGSTKSZ, whichever is larger.
constexpr usize ALT_STACK_SIZE = SIGSTKSZ > 65'536 ? SIGSTKSZ : 65'536;

// Async-signal-safe write of a NUL-terminated string to a file descriptor.
void WriteStr(int file_desc, char const* str) {
  usize len = 0;
  while (str[len] != '\0') {
    len++;
  }
  // NOLINTNEXTLINE(bugprone-signal-handler,cert-sig30-c)
  static_cast<void>(write(file_desc, str, len));
}

// Async-signal-safe hex formatter — writes "0x" + 16 hex digits, no heap allocation.
void WriteHex(int file_desc, u64 val) {
  // NOLINTNEXTLINE(modernize-avoid-c-arrays)
  char buf[18] = {'0', 'x'};
  for (int idx = 17; idx >= 2; --idx) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
    buf[idx] = "0123456789abcdef"[val & 0xFU];
    val >>= 4U;
  }
  // NOLINTNEXTLINE(bugprone-signal-handler,cert-sig30-c)
  static_cast<void>(write(file_desc, buf, sizeof(buf)));
}

// Async-signal-safe decimal integer writer for small values.
void WriteInt(int file_desc, int val) {
  // NOLINTNEXTLINE(modernize-avoid-c-arrays)
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

  // NOLINTNEXTLINE(bugprone-signal-handler,cert-sig30-c)
  static_cast<void>(write(file_desc, buf, pos));
}

// NOLINTNEXTLINE(cert-err33-c,cert-dcl03-c,cppcoreguidelines-avoid-non-const-global-variables)
void CrashHandler(int sig, siginfo_t* info, void* ctx) {
  // SA_RESETHAND already reset us to SIG_DFL — the raise() at the end
  // produces a proper core dump via the default handler.

  int const file_desc = STDERR_FILENO;
  WriteStr(file_desc, "\n============================================\n");
  WriteStr(file_desc, sig == SIGSEGV ? "*** SIGSEGV ***\n" : "*** SIGABRT ***\n");

  // Faulting address — null deref (0x0), heap corruption, or stack overflow
  WriteStr(file_desc, "Faulting address (si_addr): ");
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
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
  // NOLINTNEXTLINE(modernize-avoid-c-arrays)
  void* frames[MAX_BACKTRACE_FRAMES];
  int const depth = backtrace(frames, MAX_BACKTRACE_FRAMES);
  if (depth > 0) {
    backtrace_symbols_fd(frames, depth, file_desc);
  } else {
    WriteStr(file_desc, "  (no frames — stack is corrupted)\n");
  }

  WriteStr(file_desc, "============================================\n");
  WriteStr(file_desc, "Use: addr2line -e ./Lancet2 -f <RIP address>\n");
  WriteStr(file_desc, "============================================\n");

  // Re-raise with default handler to produce a core dump.
  // NOLINTNEXTLINE(cert-err33-c)
  static_cast<void>(raise(sig));
}

}  // namespace

namespace lancet::base {

void InstallCrashHandler() {
  // Alternate signal stack (64KB) — the handler runs here even if the crash
  // is caused by stack overflow on the thread's normal stack.
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables,modernize-avoid-c-arrays)
  static char alt_stack_mem[ALT_STACK_SIZE];
  stack_t alt_stack{};
  alt_stack.ss_sp = static_cast<void*>(alt_stack_mem);
  alt_stack.ss_size = sizeof(alt_stack_mem);
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

}  // namespace lancet::base

#else  // Non-POSIX platforms

namespace lancet::base {

void InstallCrashHandler() {
  // No-op on platforms without POSIX signal support.
}

}  // namespace lancet::base

#endif
