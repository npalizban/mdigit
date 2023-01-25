
#ifndef _MATH_LIBRARY_H
#define _MATH_LIBRARY_H

#include <stdint.h>

#define DIV_KNUTH_MAX_DIGITS  10

void negate_m(uint32_t* a, uint32_t m);
void lshift_m(uint32_t* n, uint32_t m, uint32_t shift);
void rshift_m(uint32_t* a, uint32_t m, uint32_t shift, uint32_t sext);
void add_m(uint32_t* a, uint32_t* b, uint32_t m, uint32_t* res);
void add3_m(uint32_t* a, uint32_t* b, uint32_t* c, uint32_t m, uint32_t* res);
void sub_m(uint32_t* a, uint32_t* b, uint32_t m, uint32_t* res);
void mul_m_n(uint32_t* u, uint32_t* v, uint32_t m, uint32_t n, uint32_t sign, uint32_t* res);
void mad_m_n(uint32_t* u, uint32_t* v, uint32_t m, uint32_t n, uint32_t sign, uint32_t* res);
void div_m_n(uint32_t* u, uint32_t* v, uint32_t m, uint32_t n, uint32_t* q, uint32_t* r);
void div_m_n_fract_k_f(uint32_t* u, uint32_t* v, uint32_t m, uint32_t n, uint32_t k, uint32_t f, uint32_t* q);

#endif
