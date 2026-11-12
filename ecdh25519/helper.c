#include <stddef.h> // for size_t
#include <stdint.h> // for int32_t
#include "helper.h"

uint32_t equal(uint32_t a, uint32_t b) /* 16-bit inputs */
{
    uint32_t x = a ^ b; /* 0: yes; 1..65535: no */
    x -= 1; /* 4294967295: yes; 0..65534: no */
    x >>= 31; /* 1: yes; 0: no */
    return x;
}

uint32_t ge(uint32_t a, uint32_t b) /* 16-bit inputs */
{
    uint32_t x = a;
    x -= (uint32_t) b; /* 0..65535: yes; 4294901761..4294967295: no */
    x >>= 31; /* 0: yes; 1: no */
    x ^= 1; /* 1: yes; 0: no */
    return x;
}

uint32_t times19(uint32_t a)
{
    return (a << 4) + (a << 1) + a;
}

uint32_t times38(uint32_t a)
{
    return (a << 5) + (a << 2) + (a << 1);
}

void nbyte_cmov(unsigned char *dst, const unsigned char *src, size_t n, unsigned char b)
{
    b = -b;
    size_t i;
    for (i = 0; i < n; i++) {
        dst[i] = (b & src[i]) ^ (~b & dst[i]);
    }
}

void nbyte_cswap(unsigned char *c, unsigned char *d, size_t n, unsigned char b)
{
    b = -b;
    size_t i;
    unsigned char tmp;
    for (i = 0; i < n; i++) {
        tmp = (unsigned char)(b & (c[i] ^ d[i]));
        c[i] ^= tmp;
        d[i] ^= tmp;
    }
}
