#ifndef UTILITY_LEAK_CHECK_HPP
#define UTILITY_LEAK_CHECK_HPP

#ifdef __cplusplus
extern "C" {
#endif

void leak_check_print();
void leak_check_init();
void leak_check_end();

//- fortran interfaces
void leak_check_print_();
void leak_check_init_();
void leak_check_end_();

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
class leakCheck {
 public:
  leakCheck() { leak_check_init(); }

  ~leakCheck() {
    leak_check_end();
    leak_check_print();
  }
};
#endif

#endif
