#ifndef HELPER_H
#define HELPER_H

uint32_t equal(uint32_t a, uint32_t b);
uint32_t ge(uint32_t a, uint32_t b);
uint32_t times19(uint32_t a);
uint32_t times38(uint32_t a);

void nbyte_cmov(unsigned char *dst, const unsigned char *src, size_t n, unsigned char b);
void nbyte_cswap(unsigned char *c, unsigned char *d, size_t n, unsigned char b);

#endif

