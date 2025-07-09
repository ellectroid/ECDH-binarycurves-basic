#ifndef GALOIS_FIELD2_H
#define GALOIS_FIELD2_H

#define GF2_VECTOR_MAX_BYTELEN (32UL)

long gf2_degree_lsb(const unsigned char* in, unsigned long bytelen);

void gf2_multiply_lsb(const unsigned char* in1,
                      const unsigned char* in2,
                      unsigned char* out,
                      unsigned long bytelen);

void gf2_reduce_lsb(unsigned char* inout_reducible,
                    unsigned long reducible_bytelen,
                    const unsigned char* in_reducer,
                    unsigned long reducer_bytelen);

void gf2_lshift_lsb(unsigned char*       dst,
                       const unsigned char* src,
                       unsigned long        bytelen,
                       unsigned int         shift_bits);

void gf2_binary_inverse_lsb(
    const unsigned char* in,
    unsigned char*       out,
    unsigned long        bytelen,
    const unsigned char* modulus);

#endif // GALOIS_FIELD2_H
