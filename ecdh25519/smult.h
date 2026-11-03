#ifndef SMULT_H
#define SMULT_H

#include "group.h"

int crypto_scalarmult(unsigned char *ss, const unsigned char *sk, const unsigned char *pk);

int crypto_scalarmult_base(unsigned char *pk, const unsigned char *sk);

void double_and_add(group_ge *k, const group_ge *p, const unsigned char t[32]);

#endif
