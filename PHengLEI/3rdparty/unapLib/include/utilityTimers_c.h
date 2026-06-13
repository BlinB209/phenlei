#ifndef UTILITY_TIMERS_C_H
#define UTILITY_TIMERS_C_H

#ifdef __cplusplus
extern "C" {
#endif

void TimersStart(const char *name);
void TimersEnd(const char *name);
void TimersPrintAll();
void TimersPrint(const char *name);

//- interface for Fortran
void timersstart_(const char *name);
void timersend_(const char *name);
void timersprintall_();
void timersprint_(const char *name);

#ifdef __cplusplus
}
#endif

#endif
