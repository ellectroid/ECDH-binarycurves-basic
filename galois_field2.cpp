#include "galois_field2.h"

long gf2_degree_lsb(const unsigned char* in, unsigned long bytelen) {
    long i;
    for (i = (long)bytelen - 1; i >= 0; --i) {
        if (in[i]) {
            int bit;
            for (bit = 7; bit >= 0; --bit) {
                if (in[i] & (1 << bit))
                    return (long)(i * 8 + bit);
            }
        }
    }
    return -1;
}

void gf2_multiply_lsb(const unsigned char* in1, const unsigned char* in2, unsigned char* out, unsigned long bytelen) {
    unsigned long total_bits = bytelen * 8;
    unsigned long out_bytelen = bytelen * 2;
    unsigned long i, j;

    for (i = 0; i < out_bytelen; ++i)
        out[i] = 0;

    for (i = 0; i < total_bits; ++i) {
        unsigned char bit1 = (in1[i >> 3] >> (i & 7)) & 1;
        if (!bit1) continue;

        for (j = 0; j < total_bits; ++j) {
            unsigned char bit2 = (in2[j >> 3] >> (j & 7)) & 1;
            if (!bit2) continue;

            unsigned long k = i + j;
            out[k >> 3] ^= (1U << (k & 7));
        }
    }
}

void gf2_reduce_lsb(unsigned char* inout_reducible, unsigned long reducible_bytelen, const unsigned char* in_reducer, unsigned long reducer_bytelen) {
    long deg_r = gf2_degree_lsb(inout_reducible, reducible_bytelen);
    long deg_d = gf2_degree_lsb(in_reducer, reducer_bytelen);

    if (deg_r < 0 || deg_d < 0)
        return;

    while (deg_r >= deg_d) {
        unsigned long shift = (unsigned long)(deg_r - deg_d);
        unsigned long byte_shift = shift >> 3;
        unsigned long bit_shift = shift & 7;
        unsigned long i;

        for (i = 0; i < reducer_bytelen; ++i) {
            if (i + byte_shift < reducible_bytelen) {
                unsigned char lo = in_reducer[i] << bit_shift;
                inout_reducible[i + byte_shift] ^= lo;
            }
            if (bit_shift && (i + byte_shift + 1) < reducible_bytelen) {
                unsigned char hi = in_reducer[i] >> (8 - bit_shift);
                inout_reducible[i + byte_shift + 1] ^= hi;
            }
        }

        deg_r = gf2_degree_lsb(inout_reducible, reducible_bytelen);
    }
}

void gf2_lshift_lsb(unsigned char*       dst,
                       const unsigned char* src,
                       unsigned long        bytelen,
                       unsigned int         shift_bits)
{
    const unsigned int byte_shift = shift_bits >> 3;
    const unsigned int bit_shift  = shift_bits & 0x7;

    for (unsigned long i = 0; i < bytelen; ++i)
        dst[i] = 0;

    for (unsigned long i = 0; i + byte_shift < bytelen; ++i)
        dst[i + byte_shift] = src[i];

    if (bit_shift != 0) {
        for (long i = bytelen - 1; i > 0; --i) {
            dst[i] = (dst[i] << bit_shift) | (dst[i - 1] >> (8 - bit_shift));
        }
        dst[0] <<= bit_shift;
    }
}

void gf2_binary_inverse_lsb(
    const unsigned char* in,
    unsigned char*       out,
    unsigned long        bytelen,
    const unsigned char* modulus)
{
    alignas(8) unsigned char temp1[2 * GF2_VECTOR_MAX_BYTELEN] = {};
    alignas(8) unsigned char temp2[2 * GF2_VECTOR_MAX_BYTELEN] = {};
    alignas(8) unsigned char temp3[2 * GF2_VECTOR_MAX_BYTELEN] = {};
    alignas(8) unsigned char temp4[2 * GF2_VECTOR_MAX_BYTELEN] = {};
    alignas(8) unsigned char temp5[2 * GF2_VECTOR_MAX_BYTELEN] = {};

    long deg1, deg2, shift;
    unsigned long i;

    for (i = 0; i < bytelen; ++i) {
        temp1[i] = in[i];
        temp2[i] = modulus[i];
        temp3[i] = 0;
        temp4[i] = 0;
    }
    temp3[0] = 1;

    while (true) {
        deg1 = gf2_degree_lsb(temp1, bytelen);
        if (deg1 < 0 || deg1 == 0) break;

        deg2 = gf2_degree_lsb(temp2, bytelen);
        shift = deg1 - deg2;

        if (shift < 0) {
            for (i = 0; i < bytelen; ++i) {
                unsigned char t1 = temp1[i];
                unsigned char t2 = temp3[i];
                temp1[i] = temp2[i];
                temp2[i] = t1;
                temp3[i] = temp4[i];
                temp4[i] = t2;
            }
            shift = -shift;
        }

        gf2_lshift_lsb(temp5, temp2, bytelen, (unsigned)shift);
        for (i = 0; i < bytelen; ++i)
            temp1[i] ^= temp5[i];

        gf2_lshift_lsb(temp5, temp4, bytelen, (unsigned)shift);
        for (i = 0; i < bytelen; ++i)
            temp3[i] ^= temp5[i];
    }

    for (i = 0; i < bytelen; ++i)
        out[i] = temp3[i];

    if (gf2_degree_lsb(out, bytelen) >= gf2_degree_lsb(modulus, bytelen)) {
        gf2_reduce_lsb(out, bytelen, modulus, bytelen);
    }
}
