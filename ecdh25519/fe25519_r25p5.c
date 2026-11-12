#include <stdint.h>
#include <string.h> // for memcpy in pack()
#include "fe25519.h"

const fe25519 fe25519_zero    = {{0}};
const fe25519 fe25519_one     = {{1}};
const fe25519 fe25519_two     = {{2}};
/* sqrt(-1) */
const fe25519 fe25519_sqrtm1  = {{34513072, 25610706, 9377949, 3500415, 12389472, 33281959, 41962654, 31548777, 326685, 11406482}};
/* -sqrt(-1) */
const fe25519 fe25519_msqrtm1 = {{32595773, 7943725, 57730914, 30054016, 54719391, 272472, 25146209, 2005654, 66782178, 22147949}};
/* -1 */
const fe25519 fe25519_m1      = {{67108844, 33554431, 67108863, 33554431, 67108863, 33554431, 67108863, 33554431, 67108863, 33554431}};

static const int OFFSET[11] = {0, 26, 51, 77, 102, 128, 153, 179, 204, 230, 255};
#define WIDTH(i) (((i) & 1) ? 25 : 26)

static void carry(fe25519 *r, int64_t t[10])
{
    int pass, i;
    for (pass = 0; pass < 2; pass++) {
        for (i = 0; i < 10; i++) {
            int64_t c = t[i] >> WIDTH(i);
            t[i] -= c * ((int64_t)1 << WIDTH(i));
            if (i == 9) t[0] += 19 * c;
            else        t[i + 1] += c;
        }
    }
    for (i = 0; i < 10; i++) r->v[i] = (int32_t)t[i];
}

void fe25519_unpack(fe25519 *r, const unsigned char x[32])
{
    int i, k;
    for (i = 0; i < 10; i++) {
        uint64_t w = 0;
        int byte = OFFSET[i] >> 3;
        for (k = 0; k < 5 && byte + k < 32; k++)
            w |= (uint64_t)x[byte + k] << (8 * k);
        r->v[i] = (int32_t)((w >> (OFFSET[i] & 7)) & ((1u << WIDTH(i)) - 1));
    }
}

static void sub_p_if_ge(unsigned char v[33])
{
    unsigned char w[33], keep;
    int i, d, borrow = 0;
    for (i = 0; i < 33; i++) {
        int pi = (i == 0) ? 0xED : (i < 31) ? 0xFF : (i == 31) ? 0x7F : 0x00;
        d = v[i] - pi - borrow;
        borrow = (d >> 8) & 1;
        w[i] = (unsigned char)d;
    }
    keep = (unsigned char)(borrow - 1);          /* 0xFF if v >= p */
    for (i = 0; i < 33; i++) v[i] ^= keep & (v[i] ^ w[i]);
}

void fe25519_pack(unsigned char r[32], const fe25519 *x)
{
    int64_t t[10];
    unsigned char v[33] = {0};
    fe25519 h;
    int i, k;

    for (i = 0; i < 10; i++) t[i] = x->v[i];
    carry(&h, t);
    for (i = 0; i < 10; i++) t[i] = h.v[i];
    carry(&h, t);

    for (i = 0; i < 10; i++) {
        uint64_t acc = (uint64_t)(uint32_t)h.v[i] << (OFFSET[i] & 7);
        for (k = OFFSET[i] >> 3; k < 33; k++) {
            acc += v[k];
            v[k] = (unsigned char)acc;
            acc >>= 8;
        }
    }
    sub_p_if_ge(v);
    sub_p_if_ge(v);
    memcpy(r, v, 32);
}


// no carries
void fe25519_add(fe25519 *r, const fe25519 *x, const fe25519 *y)
{
    int i;
    for (i = 0; i < 10; i++) r->v[i] = x->v[i] + y->v[i];
}

void fe25519_sub(fe25519 *r, const fe25519 *x, const fe25519 *y)
{
    int i;
    for (i = 0; i < 10; i++) r->v[i] = x->v[i] - y->v[i];
}

void fe25519_neg(fe25519 *r, const fe25519 *x)
{
    for (int i = 0; i < 10; i++) r->v[i] = -x->v[i];
}

// multiplications -> 100 prods
void fe25519_mul(fe25519 *r, const fe25519 *x, const fe25519 *y)
{
    int64_t t[10] = {0};
    int64_t p;
    int i, j, k;
    for (i = 0; i < 10; i++) {
        for (j = 0; j < 10; j++) {
            p = (int64_t)x->v[i] * y->v[j];
            if (i & j & 1) p *= 2;
            k = i + j;
            if (k >= 10) { k -= 10; p *= 19; }
            t[k] += p;
        }
    }
    carry(r, t);
}

// squaring -> 55 prods
void fe25519_square(fe25519 *r, const fe25519 *x)
{
    int64_t t[10] = {0};
    int64_t p;
    int i, j, k;
    for (i = 0; i < 10; i++) {
        for (j = i; j < 10; j++) {
            p = (int64_t)x->v[i] * x->v[j];
            if (i != j) p *= 2;
            if (i & j & 1) p *= 2;
            k = i + j;
            if (k >= 10) { k -= 10; p *= 19; }
            t[k] += p;
        }
    }
    carry(r, t);
}

void fe25519_mul_small(fe25519 *r, const fe25519 *x, int32_t c)
{
    int i;
    int64_t t[10];
    for (i = 0; i < 10; i++) t[i] = (int64_t)x->v[i] * c;
    carry(r, t);
}
