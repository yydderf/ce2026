#include <stdio.h>
#include "fe25519.h"
#include "helper.h"

void fe25519_cmov(fe25519 *r, const fe25519 *x, unsigned char b)
{
    nbyte_cmov((unsigned char *)r, (const unsigned char *)x, sizeof(fe25519), b);
}

void fe25519_cswap(fe25519 *r, fe25519 *x, unsigned char b)
{
    // nbyte_cmov(u, v, sizeof(fe25519), b);
    nbyte_cswap((unsigned char *)r, (unsigned char *)x, sizeof(fe25519), b);
}

int fe25519_iseq(const fe25519 *x, const fe25519 *y)
{
    // fe25519 t1,t2;
    // t1 = *x;
    // t2 = *y;
    // fe25519_freeze(&t1);
    // fe25519_freeze(&t2);
    // for(i=0;i<32;i++)
    //     r |= (1-equal(t1.v[i],t2.v[i]));
    // return 1-r;
    
    // for arbitrary limb sizes
    unsigned char a[32], b[32];
    uint32_t d = 0;
    int i;

    fe25519_pack(a, x);
    fe25519_pack(b, y);
    for (i = 0; i < 32; i++) {
        d |= (uint32_t)(a[i] ^ b[i]);
    }
    return (int)equal(d, 0);
}

int fe25519_iszero(const fe25519 *x)
{
    return fe25519_iseq(x, &fe25519_zero);
}

int fe25519_isone(const fe25519 *x)
{
    return fe25519_iseq(x, &fe25519_one);
}

// return true if x has LSB set
int fe25519_isnegative(const fe25519 *x)
{
    unsigned char b[32];
    fe25519_pack(b, x);
    return b[0] & 1;
    // fe25519 t = *x;
    // fe25519_freeze(&t);
    // return t.v[0] & 1;
}

void fe25519_double(fe25519 *r, const fe25519 *x)
{
    fe25519_add(r, x, x);        /* reduction policy comes from fe25519_add */
}

void fe25519_pow2523(fe25519 *r, const fe25519 *x)
{
    fe25519 z2;
    fe25519 z9;
    fe25519 z11;
    fe25519 z2_5_0;
    fe25519 z2_10_0;
    fe25519 z2_20_0;
    fe25519 z2_50_0;
    fe25519 z2_100_0;
    fe25519 t;
    int i;

    /* 2 */ fe25519_square(&z2,x);
    /* 4 */ fe25519_square(&t,&z2);
    /* 8 */ fe25519_square(&t,&t);
    /* 9 */ fe25519_mul(&z9,&t,x);
    /* 11 */ fe25519_mul(&z11,&z9,&z2);
    /* 22 */ fe25519_square(&t,&z11);
    /* 2^5 - 2^0 = 31 */ fe25519_mul(&z2_5_0,&t,&z9);

    /* 2^6 - 2^1 */ fe25519_square(&t,&z2_5_0);
    /* 2^10 - 2^5 */ for (i = 1;i < 5;i++) { fe25519_square(&t,&t); }
    /* 2^10 - 2^0 */ fe25519_mul(&z2_10_0,&t,&z2_5_0);

    /* 2^11 - 2^1 */ fe25519_square(&t,&z2_10_0);
    /* 2^20 - 2^10 */ for (i = 1;i < 10;i++) { fe25519_square(&t,&t); }
    /* 2^20 - 2^0 */ fe25519_mul(&z2_20_0,&t,&z2_10_0);

    /* 2^21 - 2^1 */ fe25519_square(&t,&z2_20_0);
    /* 2^40 - 2^20 */ for (i = 1;i < 20;i++) { fe25519_square(&t,&t); }
    /* 2^40 - 2^0 */ fe25519_mul(&t,&t,&z2_20_0);

    /* 2^41 - 2^1 */ fe25519_square(&t,&t);
    /* 2^50 - 2^10 */ for (i = 1;i < 10;i++) { fe25519_square(&t,&t); }
    /* 2^50 - 2^0 */ fe25519_mul(&z2_50_0,&t,&z2_10_0);

    /* 2^51 - 2^1 */ fe25519_square(&t,&z2_50_0);
    /* 2^100 - 2^50 */ for (i = 1;i < 50;i++) { fe25519_square(&t,&t); }
    /* 2^100 - 2^0 */ fe25519_mul(&z2_100_0,&t,&z2_50_0);

    /* 2^101 - 2^1 */ fe25519_square(&t,&z2_100_0);
    /* 2^200 - 2^100 */ for (i = 1;i < 100;i++) { fe25519_square(&t,&t); }
    /* 2^200 - 2^0 */ fe25519_mul(&t,&t,&z2_100_0);

    /* 2^201 - 2^1 */ fe25519_square(&t,&t);
    /* 2^250 - 2^50 */ for (i = 1;i < 50;i++) { fe25519_square(&t,&t); }
    /* 2^250 - 2^0 */ fe25519_mul(&t,&t,&z2_50_0);

    /* 2^251 - 2^1 */ fe25519_square(&t,&t);
    /* 2^252 - 2^2 */ fe25519_square(&t,&t);
    /* 2^252 - 3 */ fe25519_mul(r,&t,x);
}

void fe25519_invsqrt(fe25519 *r, const fe25519 *x)
{
    fe25519 den2, den3, den4, den6, chk, t, t2;
    int b;

    fe25519_square(&den2, x);
    fe25519_mul(&den3, &den2, x);

    fe25519_square(&den4, &den2);
    fe25519_mul(&den6, &den2, &den4);
    fe25519_mul(&t, &den6, x); // r is now x^7

    fe25519_pow2523(&t, &t);
    fe25519_mul(&t, &t, &den3);

    fe25519_square(&chk, &t);
    fe25519_mul(&chk, &chk, x);

    fe25519_mul(&t2, &t, &fe25519_sqrtm1);
    b = 1 - fe25519_isone(&chk);

    fe25519_cmov(&t, &t2, b);

    *r = t;
}
