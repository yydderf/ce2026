#ifndef FE25519_H
#define FE25519_H

#include <stdint.h>
#include <stddef.h>

#if defined(FE25519_R25P5) && defined(FE25519_R8)
#error "define only one of FE25519_R25P5 / FE25519_R8"
#elif defined(FE25519_R25P5)
#include "fe25519_r25p5.h"
#else
#include "fe25519_r8.h"
#endif

extern const fe25519 fe25519_zero;
extern const fe25519 fe25519_one;
extern const fe25519 fe25519_two;

extern const fe25519 fe25519_sqrtm1;
extern const fe25519 fe25519_msqrtm1;
extern const fe25519 fe25519_m1;

void fe25519_unpack(fe25519 *r, const unsigned char x[32]);
void fe25519_pack(unsigned char r[32], const fe25519 *x);
void fe25519_add(fe25519 *r, const fe25519 *x, const fe25519 *y);
void fe25519_sub(fe25519 *r, const fe25519 *x, const fe25519 *y);
void fe25519_mul(fe25519 *r, const fe25519 *x, const fe25519 *y);
void fe25519_mul_small(fe25519 *r, const fe25519 *x, int32_t c); /* 0 <= c < 2^17 */
void fe25519_square(fe25519 *r, const fe25519 *x);
void fe25519_cmov(fe25519 *r, const fe25519 *x, unsigned char b);
void _fe25519_cmov(fe25519 *r, const fe25519 *x, unsigned char b);
void fe25519_cswap(fe25519 *r, fe25519 *x, unsigned char b);

void fe25519_freeze(fe25519 *r);

int fe25519_iszero(const fe25519 *x);
int fe25519_isone(const fe25519 *x);
int fe25519_isnegative(const fe25519 *x);
int fe25519_iseq(const fe25519 *x, const fe25519 *y);
int fe25519_iseq_vartime(const fe25519 *x, const fe25519 *y);
void fe25519_setone(fe25519 *r);
void fe25519_setzero(fe25519 *r);
void fe25519_neg(fe25519 *r, const fe25519 *x);
void fe25519_double(fe25519 *r, const fe25519 *x);
void fe25519_triple(fe25519 *r, const fe25519 *x);
void fe25519_invert(fe25519 *r, const fe25519 *x);
void fe25519_pow2523(fe25519 *r, const fe25519 *x);
void fe25519_invsqrt(fe25519 *r, const fe25519 *x);
void fe25519_print(const fe25519 *x);
unsigned char fe25519_getparity(const fe25519 *x);

#endif
