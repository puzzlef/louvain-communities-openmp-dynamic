#pragma once
#include <cassert>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include "_iostream.hxx"
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=1
#include <unistd.h>
#include <signal.h>
#include <execinfo.h>
#endif




// BUILD
// -----
// Build modes.

#ifndef BUILD_RELEASE
#define BUILD_RELEASE 0
#define BUILD_ERROR   1
#define BUILD_WARNING 2
#define BUILD_INFO    3
#define BUILD_DEBUG   4
#define BUILD_TRACE   5
#endif




// ASSERT
// ------

#ifndef ASSERT
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_ERROR
#define ASSERT(exp)           assert(exp)
#define ASSERT_THAT(exp, msg) assert((exp) && (msg))
#else
#define ASSERT(exp)
#define ASSERT_THAT(exp, msg)
#endif
#endif




// PERFORM
// -------

#ifndef PEFORME
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_ERROR
#define PERFORME(...) __VA_ARGS__
#else
#define PERFORME(...)
#endif
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_WARNING
#define PERFORMW(...) __VA_ARGS__
#else
#define PERFORMW(...)
#endif
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_INFO
#define PERFORMI(...) __VA_ARGS__
#else
#define PERFORMI(...)
#endif
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_DEBUG
#define PERFORMD(...) __VA_ARGS__
#else
#define PERFORMD(...)
#endif
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_TRACE
#define PERFORMT(...) __VA_ARGS__
#else
#define PERFORMT(...)
#endif
#endif




// PRINT
// -----

#ifndef FPRINTFE
#define FPRINTFE(...) PERFORME(fprintf(__VA_ARGS__))
#define FPRINTFW(...) PERFORMW(fprintf(__VA_ARGS__))
#define FPRINTFI(...) PERFORMI(fprintf(__VA_ARGS__))
#define FPRINTFD(...) PERFORMD(fprintf(__VA_ARGS__))
#define FPRINTFT(...) PERFORMT(fprintf(__VA_ARGS__))
#endif

#ifndef PRINTFE
#define PRINTFE(...) PERFORME(printf(__VA_ARGS__))
#define PRINTFW(...) PERFORMW(printf(__VA_ARGS__))
#define PRINTFI(...) PERFORMI(printf(__VA_ARGS__))
#define PRINTFD(...) PERFORMD(printf(__VA_ARGS__))
#define PRINTFT(...) PERFORMT(printf(__VA_ARGS__))
#endif

#ifndef WRITEE
#define WRITEE(...) PERFORME(write(__VA_ARGS__))
#define WRITEW(...) PERFORMW(write(__VA_ARGS__))
#define WRITEI(...) PERFORMI(write(__VA_ARGS__))
#define WRITED(...) PERFORMD(write(__VA_ARGS__))
#define WRITET(...) PERFORMT(write(__VA_ARGS__))
#endif

#ifndef PRINTE
#define PRINTE(...) PERFORME(print(__VA_ARGS__))
#define PRINTW(...) PERFORMW(print(__VA_ARGS__))
#define PRINTI(...) PERFORMI(print(__VA_ARGS__))
#define PRINTD(...) PERFORMD(print(__VA_ARGS__))
#define PRINTT(...) PERFORMT(print(__VA_ARGS__))
#endif

#ifndef PRINTLNE
#define PRINTLNE(...) PERFORME(println(__VA_ARGS__))
#define PRINTLNW(...) PERFORMW(println(__VA_ARGS__))
#define PRINTLNI(...) PERFORMI(println(__VA_ARGS__))
#define PRINTLND(...) PERFORMD(println(__VA_ARGS__))
#define PRINTLNT(...) PERFORMT(println(__VA_ARGS__))
#endif




// LOG
// ---

#ifndef LOG
#define LOG(...) do { print(timeNow()); printf(" " __VA_ARGS__); } while (0)
#endif

#ifndef LOGE
#define LOGE(...) PERFORME(LOG(__VA_ARGS__))
#define LOGW(...) PERFORMW(LOG(__VA_ARGS__))
#define LOGI(...) PERFORMI(LOG(__VA_ARGS__))
#define LOGD(...) PERFORMD(LOG(__VA_ARGS__))
#define LOGT(...) PERFORMT(LOG(__VA_ARGS__))
#endif




// ON SIGNAL
// ---------

#define STACK_TRACE_SIZE 32

void on_sigsegv(int sig) {
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=1
  void *entries[STACK_TRACE_SIZE];
  size_t n = backtrace(entries, STACK_TRACE_SIZE);
  fprintf(stderr, "ERROR: SIGNAL %d:\n", sig);
  backtrace_symbols_fd(entries, n, STDERR_FILENO);
  exit(1);
#endif
}
// - https://stackoverflow.com/a/77336/1413259

void install_sigsegv() {
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=1
  signal(SIGSEGV, on_sigsegv);
#endif
}
// - https://stackoverflow.com/a/77336/1413259
