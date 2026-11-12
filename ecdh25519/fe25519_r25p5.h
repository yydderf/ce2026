#ifndef FE25519_R25P5_H
#define FE25519_R25P5_H

#ifndef FE25519_H
#error "include fe25519.h instead of fe25519_r25_5.h"
#endif

#include <stdint.h>

#define FE25519_LIMBS 10

typedef struct
{
    int32_t v[FE25519_LIMBS];
} fe25519;

#endif
