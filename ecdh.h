#ifndef ECDH_H_
#define ECDH_H_

#include "elliptic_curve.h"

void ecdh_generate_public_key(EllipticCurve *curve,
		unsigned char *in_private_key, unsigned char *out_public_key);
int ecdh_public_key_verify(EllipticCurve *curve, unsigned char *public_key);
void ecdh_generate_shared_secret(EllipticCurve *curve,
		unsigned char *in_private_key, unsigned char *in_public_key,
		unsigned char *out_shared_secret);

typedef struct
	alignas(8) {
		alignas(8) unsigned char private_key[GF2_VECTOR_MAX_BYTELEN];
		union {
			alignas(8) unsigned char public_key[2 * GF2_VECTOR_MAX_BYTELEN];
			EllipticCurvePoint public_key_point;
		};
		alignas(8) unsigned char shared_secret[GF2_VECTOR_MAX_BYTELEN];
	} ecdh_keygroup_t;

#endif /* ECDH_H_ */
