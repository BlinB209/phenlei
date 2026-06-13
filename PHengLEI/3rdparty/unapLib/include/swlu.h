#ifndef LIB_SWLU_H
#define LIB_SWLU_H

#include <stdint.h>

#ifdef __cplusplus

extern "C" int swlu_backtrace(uint64_t*, int);
extern "C" int swlu_backtrace_symbols(uint64_t*, int);
extern "C" void swlu_show_callstack(int);

extern "C" void swlu_build_symtable();
extern "C" void swlu_debug_init();

extern "C" void swlu_prof_init();
extern "C" void swlu_prof_print();
extern "C" void swlu_prof_start();
extern "C" void swlu_prof_stop();

extern "C" void samprof_init();
extern "C" void samprof_print();
extern "C" void samprof_start();
extern "C" void samprof_stop();

#else

int swlu_backtrace(uint64_t*, int);
int swlu_backtrace_symbols(uint64_t*, int);
void swlu_show_callstack(int);

void swlu_build_symtable();
void swlu_debug_init();

void swlu_prof_init();
void swlu_prof_print();
void swlu_prof_start();
void swlu_prof_stop();

void samprof_init();
void samprof_print();
void samprof_start();
void samprof_stop();

#endif

#endif
