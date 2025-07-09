#ifndef ELLIPTIC_CURVE_H_
#define ELLIPTIC_CURVE_H_

#include "galois_field2.h"

typedef struct alignas(8){
    unsigned char a[GF2_VECTOR_MAX_BYTELEN];        // Curve coefficient a
    unsigned char b[GF2_VECTOR_MAX_BYTELEN];        // Curve coefficient b
    unsigned char xG[GF2_VECTOR_MAX_BYTELEN];       // x-coordinate of base point G
    unsigned char yG[GF2_VECTOR_MAX_BYTELEN];       // y-coordinate of base point G
    unsigned char modulus[GF2_VECTOR_MAX_BYTELEN];  // Field reduction polynomial f(x)
    unsigned char order[GF2_VECTOR_MAX_BYTELEN];    // Order of base point G
    unsigned char cofactor[GF2_VECTOR_MAX_BYTELEN]; // Curve cofactor (usually small: 2, 4)
    unsigned char curve_name_ascii[16];				// optional name
    unsigned long field_size_bytes;                 // Actual size of field element in bytes
    unsigned long binary_degree;					// 0 for prime fields, non-zero for binary
}EllipticCurve;

typedef struct alignas(8){
		unsigned char point_mem[2*GF2_VECTOR_MAX_BYTELEN];
}EllipticCurvePoint; //dynamically sized object, layout depends on curve's field_size_bytes

unsigned long elliptic_curve_get_maximum_vector_bytelen();

void elliptic_curve_binary_point_double(
    EllipticCurve* curve,
    EllipticCurvePoint* out,
    EllipticCurvePoint* in);
void elliptic_curve_binary_point_add(EllipticCurve* curve, EllipticCurvePoint* out, EllipticCurvePoint* in1, EllipticCurvePoint* in2);

void elliptic_curve_binary_point_multiply(EllipticCurve *curve, EllipticCurvePoint *out,
		EllipticCurvePoint *in, const unsigned char *exp,
		unsigned long bytelen);

int elliptic_curve_binary_point_on_curve(EllipticCurve *curve,
                                         EllipticCurvePoint *point);

unsigned char* elliptic_curve_point_get_coord_x(EllipticCurve *curve, EllipticCurvePoint *point);
unsigned char* elliptic_curve_point_get_coord_y(EllipticCurve *curve, EllipticCurvePoint *point);
unsigned long elliptic_curve_point_get_coord_one_bytelen(EllipticCurve *curve);
unsigned long elliptic_curve_point_get_coord_full_bytelen(EllipticCurve *curve);

#endif /* ELLIPTIC_CURVE_H_ */
