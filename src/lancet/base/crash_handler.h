#ifndef SRC_LANCET_BASE_CRASH_HANDLER_H_
#define SRC_LANCET_BASE_CRASH_HANDLER_H_

namespace lancet::base {

/// Install a process-wide SIGSEGV/SIGABRT handler that prints diagnostic info
/// to stderr on crash: faulting address, signal code, thread ID, instruction
/// pointer (RIP), and a best-effort stack backtrace.
///
/// Call once from main() before spawning worker threads. On platforms without
/// POSIX signals (non-Linux, non-macOS) this is a no-op.
void InstallCrashHandler();

}  // namespace lancet::base

#endif  // SRC_LANCET_BASE_CRASH_HANDLER_H_
