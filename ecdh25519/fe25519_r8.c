#include <stdio.h>
#include "fe25519.h"
#include "helper.h"

const fe25519 fe25519_zero = {{0}};
const fe25519 fe25519_one  = {{1}};
const fe25519 fe25519_two  = {{2}};

/* sqrt(-1) */
const fe25519 fe25519_sqrtm1 = {{0xB0, 0xA0, 0x0E, 0x4A, 0x27, 0x1B, 0xEE, 0xC4, 0x78, 0xE4, 0x2F, 0xAD, 0x06, 0x18, 0x43, 0x2F,
                                 0xA7, 0xD7, 0xFB, 0x3D, 0x99, 0x00, 0x4D, 0x2B, 0x0B, 0xDF, 0xC1, 0x4F, 0x80, 0x24, 0x83, 0x2B}};

/* -sqrt(-1) */
const fe25519 fe25519_msqrtm1 = {{0x3D, 0x5F, 0xF1, 0xB5, 0xD8, 0xE4, 0x11, 0x3B, 0x87, 0x1B, 0xD0, 0x52, 0xF9, 0xE7, 0xBC, 0xD0,
                                  0x58, 0x28, 0x4, 0xC2, 0x66, 0xFF, 0xB2, 0xD4, 0xF4, 0x20, 0x3E, 0xB0, 0x7F, 0xDB, 0x7C, 0x54}};

/* -1 */
const fe25519 fe25519_m1 = {{236, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
                             0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f}};

static void reduce_add_sub(fe25519 *r)
{
    uint32_t t;
    int i,rep;

    for(rep=0;rep<2;rep++)
    {
        t = r->v[31] >> 7;
        r->v[31] &= 127;
        t = times19(t);
        r->v[0] += t;
        for(i=0;i<31;i++)
        {
            t = r->v[i] >> 8;
            r->v[i+1] += t;
            r->v[i] &= 255;
        }
    }
}

static void reduce_mul(fe25519 *r)
{
    uint32_t t;
    int i,rep;

    for(rep=0;rep<3;rep++)
    {
        t = r->v[31] >> 7;
        r->v[31] &= 127;
        t = times19(t);
        r->v[0] += t;
        for(i=0;i<31;i++)
        {
            t = r->v[i] >> 8;
            r->v[i+1] += t;
            r->v[i] &= 255;
        }
    }
}

/* reduction modulo 2^255-19 */
void fe25519_freeze(fe25519 *r)
{
    int i;
    uint32_t m = equal(r->v[31],127);
    for(i=30;i>0;i--)
        m &= equal(r->v[i],255);
    m &= ge(r->v[0],237);

    m = -m;

    r->v[31] -= m&127;
    for(i=30;i>0;i--)
        r->v[i] -= m&255;
    r->v[0] -= m&237;
}

void fe25519_unpack(fe25519 *r, const unsigned char x[32])
{
    int i;
    for(i=0;i<32;i++) r->v[i] = x[i];
    r->v[31] &= 127;
}

/* Assumes input x being reduced below 2^255 */
void fe25519_pack(unsigned char r[32], const fe25519 *x)
{
    int i;
    fe25519 y = *x;
    fe25519_freeze(&y);
    for(i=0;i<32;i++)
        r[i] = y.v[i];
}



// not constant time impl
// void _fe25519_cmov(fe25519 *r, const fe25519 *x, unsigned char b)
// {
//     if(b) *r = *x;
// }

void fe25519_neg(fe25519 *r, const fe25519 *x)
{
    fe25519 t = fe25519_zero;
    fe25519_sub(r, &t, x);
}

void fe25519_add(fe25519 *r, const fe25519 *x, const fe25519 *y)
{
    int i;
    for(i=0;i<32;i++) r->v[i] = x->v[i] + y->v[i];
    reduce_add_sub(r);
}

// void fe25519_double(fe25519 *r, const fe25519 *x)
// {
//     int i;
//     for(i=0;i<32;i++) r->v[i] = 2*x->v[i];
//     reduce_add_sub(r);
// }

void fe25519_sub(fe25519 *r, const fe25519 *x, const fe25519 *y)
{
    int i;
    uint32_t t[32];
    t[0] = x->v[0] + 0x1da;
    t[31] = x->v[31] + 0xfe;
    for(i=1;i<31;i++) t[i] = x->v[i] + 0x1fe;
    for(i=0;i<32;i++) r->v[i] = t[i] - y->v[i];
    reduce_add_sub(r);
}

void fe25519_mul(fe25519 *r, const fe25519 *x, const fe25519 *y)
{
    int i,j;
    uint32_t t[63];
    for(i=0;i<63;i++)t[i] = 0;

    for(i=0;i<32;i++)
        for(j=0;j<32;j++)
            t[i+j] += x->v[i] * y->v[j];

    for(i=32;i<63;i++)
        r->v[i-32] = t[i-32] + times38(t[i]);
    r->v[31] = t[31]; /* result now in r[0]...r[31] */

    reduce_mul(r);
}

void fe25519_square(fe25519 *r, const fe25519 *x)
{
    fe25519_mul(r, x, x);
}

void fe25519_mul_small(fe25519 *r, const fe25519 *x, int32_t c)
{
    fe25519 k = fe25519_zero;                  /* c as a field element (3 bytes) */
    k.v[0] = (uint32_t)c & 0xff;
    k.v[1] = ((uint32_t)c >> 8) & 0xff;
    k.v[2] = ((uint32_t)c >> 16) & 0xff;
    fe25519_mul(r, x, &k);
}

