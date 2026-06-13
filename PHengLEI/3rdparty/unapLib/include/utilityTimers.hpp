#ifndef UTILITY_TIMERS_HPP
#define UTILITY_TIMERS_HPP

#include "utilityContainer.hpp"

using UTILITY::Table;
using UTILITY::Word;

class Timers {
 public:
  static Table<Word, scalar64> timeSum;
  static Table<Word, scalar64> timeStart;
  static Table<Word, scalar64> timeEnd;
  static Table<Word, label32> count;

  void static startTimer(Word in);
  static void endTimer(Word in);
  void static startTimer(const char* in);
  static void endTimer(const char* in);

  ~Timers();

  static void printTimer();
  static void printTimer(Word in);
  static void printTimer(const char* in);

  static void maxTimeSum();
};

#if (TIMERS)
#define TIMERS_START(name) Timers::startTimer(name);
#else
#define TIMERS_START(name)
#endif

#if (TIMERS)
#define TIMERS_END(name) Timers::endTimer(name);
#else
#define TIMERS_END(name)
#endif

#if (TIMERS)
#define TIMERS_PRINT() Timers::printTimer();
#else
#define TIMERS_PRINT()
#endif

#endif
