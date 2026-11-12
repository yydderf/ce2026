#ifndef SMULT_H
#define SMULT_H

#include "group.h"
#include "table.h"

extern const group_ge group_ge_neutral;

void clamp_scalar(unsigned char *t, const unsigned char *s);

// ss = sk * pk
int crypto_scalarmult(unsigned char *ss, const unsigned char *sk, const unsigned char *pk);

// pk = sk * base_point
int crypto_scalarmult_base(unsigned char *pk, const unsigned char *sk);

// k = t * p
void double_and_add(group_ge *k, const group_ge *p, const unsigned char *t);

void double_and_add_always(group_ge *k, const group_ge *p, const unsigned char *t);

// fixed-window_scalar_mult -> get rid of doubling at runtime
// fixed-window_scalar_mult -> comp time optimization? -> for fixed-base points
void fixed_scalar_mult(group_ge *k, const unsigned char *t);

// sliding-window_scalar_mult

// montgomery_ladder -> constant time
void ladder_step(fe25519 *xr, fe25519 *zr, fe25519 *xs, fe25519 *zs);
void ladder_mult(group_ge *k, const group_ge *p, const unsigned char *t);

#endif
