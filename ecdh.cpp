#include "ecdh.h"

void ecdh_generate_public_key(EllipticCurve *curve,
		unsigned char *in_private_key, unsigned char *out_public_key) {
	EllipticCurvePoint base_point = { 0 };
	unsigned char *base_x = elliptic_curve_point_get_coord_x(curve,
			&base_point);
	unsigned char *base_y = elliptic_curve_point_get_coord_y(curve,
			&base_point);
	for (unsigned long i = 0; i < curve->field_size_bytes; i++) {
		base_x[i] = curve->xG[i];
		base_y[i] = curve->yG[i];
	}
	elliptic_curve_binary_point_multiply(curve,
			(EllipticCurvePoint*) out_public_key, &base_point, in_private_key,
			curve->field_size_bytes);
}
int ecdh_public_key_verify(EllipticCurve *curve, unsigned char *public_key) {
	unsigned long len = curve->field_size_bytes;
	unsigned char *pk_x = elliptic_curve_point_get_coord_x(curve,
			(EllipticCurvePoint*) public_key);
	unsigned char *pk_y = elliptic_curve_point_get_coord_y(curve,
			(EllipticCurvePoint*) public_key);
	unsigned int is_inf = 1;
	for (unsigned long i = 0; i < len; i++) {
		if ((pk_x[i] != 0) || (pk_y[i] != 0)) {
			is_inf = 0;
			break;
		}
	}
	if (is_inf)
		return 0; //all zeroes (infinity point)

	//Check if the point is on the curve
	int is_on_curve = elliptic_curve_binary_point_on_curve(curve,
			(EllipticCurvePoint*) public_key);
	if (!is_on_curve)
		return 0;

	//Check if the point is in the correct subgroup
	EllipticCurvePoint multiplication_result = { 0 };
	unsigned char *mr_x = elliptic_curve_point_get_coord_x(curve,
			&multiplication_result);
	unsigned char *mr_y = elliptic_curve_point_get_coord_y(curve,
			&multiplication_result);
	elliptic_curve_binary_point_multiply(curve, &multiplication_result,
			(EllipticCurvePoint*) public_key, curve->order,
			curve->field_size_bytes);
	is_inf = 1;
	for (unsigned long i = 0; i < len; i++) {
		if (mr_x[i] != 0 || mr_y[i] != 0) {
			is_inf = 0;
			break;
		}
	}
	if (!is_inf)
		return 0; //multiplication by base order didn't yield zero
	return 1;
}
void ecdh_generate_shared_secret(EllipticCurve *curve,
		unsigned char *in_private_key, unsigned char *in_public_key,
		unsigned char *out_shared_secret) {
	EllipticCurvePoint *in_public_key_point =
			(EllipticCurvePoint*) in_public_key;
	EllipticCurvePoint output_shared_secret = { 0 };
	elliptic_curve_binary_point_multiply(curve, &output_shared_secret,
			in_public_key_point, in_private_key, curve->field_size_bytes);
	unsigned long len = curve->field_size_bytes;
	for (unsigned long i = 0; i < len; i++) {
		out_shared_secret[i] = output_shared_secret.point_mem[i];
	}
}
