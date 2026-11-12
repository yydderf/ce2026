#ifndef FE25519_R8_H
#define FE25519_R8_H

#ifndef FE25519_H
#error "include fe25519.h instead of fe25519_r8.h"
#endif

#include <stdint.h>

#define FE25519_LIMBS 32

typedef struct
{
    uint32_t v[FE25519_LIMBS];
} fe25519;

#endif
