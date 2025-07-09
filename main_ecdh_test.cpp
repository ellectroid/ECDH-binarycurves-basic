#include "ecdh.h"
#include <iostream>
#include <iomanip>

void curve_configure_k233(EllipticCurve *curve);

#define K233_VECTOR_MEMLEN (32UL)

int main() {
	std::cout
			<< "Compile-time maximum GF2 element capacity (affects temp buffer sizes),\n"
					"From galois_field2.h, GF2_VECTOR_MAX_BYTELEN: "
			<< +elliptic_curve_get_maximum_vector_bytelen() << ". " << std::endl;
	std::cout << "Have at least one leading zero byte." << std::endl;


	//Create a curve
	std::cout << "\n--- Testing ECDH with K-233 (sect233k1) ---\n";

	EllipticCurve k233_curve = { 0 };
	curve_configure_k233(&k233_curve);

	//Creating data structures to hold keys
	alignas(8) unsigned char k233_test_alex_private_key[K233_VECTOR_MEMLEN] = { 0 };
	alignas(8) unsigned char k233_test_alex_public_key[2 * K233_VECTOR_MEMLEN] = { 0 };
	alignas(8) unsigned char k233_test_alex_shared_secret[K233_VECTOR_MEMLEN] = { 0 };

	alignas(8) unsigned char k233_test_bethany_private_key[K233_VECTOR_MEMLEN] = { 0 };
	alignas(8) unsigned char k233_test_bethany_public_key[2 * K233_VECTOR_MEMLEN] = { 0 };
	alignas(8) unsigned char k233_test_bethany_shared_secret[K233_VECTOR_MEMLEN] = { 0 };

	alignas(8) unsigned char k233_test_shared_secret_expected_scalar[K233_VECTOR_MEMLEN] = { 0 };
	alignas(8) unsigned char k233_test_shared_secret_expected_value[2 * K233_VECTOR_MEMLEN] = { 0 };
	alignas(8) unsigned char k233_test_shared_secret_base_point[2 * K233_VECTOR_MEMLEN] = { 0 };
	unsigned char* k233_test_shared_secret_base_point_x = elliptic_curve_point_get_coord_x(&k233_curve, (EllipticCurvePoint*)k233_test_shared_secret_base_point);
	unsigned char* k233_test_shared_secret_base_point_y = elliptic_curve_point_get_coord_y(&k233_curve, (EllipticCurvePoint*)k233_test_shared_secret_base_point);

	k233_test_alex_private_key[0] = 2;
	k233_test_bethany_private_key[0] = 3;
	k233_test_shared_secret_expected_scalar[0] = 6; //2*3 = 3*2
	for (unsigned long i = 0; i < k233_curve.field_size_bytes; ++i){
		k233_test_shared_secret_base_point_x[i] = k233_curve.xG[i];
		k233_test_shared_secret_base_point_y[i] = k233_curve.yG[i];
	}

	ecdh_generate_public_key(&k233_curve, k233_test_alex_private_key, k233_test_alex_public_key);
	ecdh_generate_public_key(&k233_curve, k233_test_bethany_private_key,
			k233_test_bethany_public_key);
	//alex_public_key[0] ^= 0x01; //messing up the key
	std::cout << "Alex's public key validity test: " << +ecdh_public_key_verify(&k233_curve, k233_test_alex_public_key) << std::endl;
	std::cout << "Bethany's public key validity test: " << +ecdh_public_key_verify(&k233_curve, k233_test_bethany_public_key) << std::endl;

	std::cout << "Alex's public key is 2G." << std::endl;
	std::cout << "Bethany's public key is 3G." << std::endl;
	std::cout << "Shared secret must be 6G." << std::endl;

	ecdh_generate_shared_secret(&k233_curve, k233_test_alex_private_key,
			k233_test_bethany_public_key, k233_test_alex_shared_secret);
	ecdh_generate_shared_secret(&k233_curve, k233_test_bethany_private_key,
			k233_test_alex_public_key, k233_test_bethany_shared_secret);

	//shared secret expected value
	elliptic_curve_binary_point_multiply(&k233_curve, (EllipticCurvePoint *)k233_test_shared_secret_expected_value,
			(EllipticCurvePoint *)k233_test_shared_secret_base_point, k233_test_shared_secret_expected_scalar,
			k233_curve.field_size_bytes);

	std::cout << "Alex's shared secret:\nX: ";
	for (unsigned long i = 0; i < k233_curve.field_size_bytes; ++i)
		std::cout << std::hex << std::setw(2) << std::setfill('0')
				<< (int) k233_test_alex_shared_secret[i];
	std::cout << std::endl;
	std::cout << "Bethany's shared secret:\nX: ";
	for (unsigned long i = 0; i < k233_curve.field_size_bytes; ++i)
		std::cout << std::hex << std::setw(2) << std::setfill('0')
				<< (int) k233_test_bethany_shared_secret[i];
	std::cout << std::endl;
	std::cout << "Expected shared secret:\nX: ";
	for (unsigned long i = 0; i < k233_curve.field_size_bytes; ++i)
			std::cout << std::hex << std::setw(2) << std::setfill('0')
					<< (int) k233_test_shared_secret_expected_value[i];
	std::cout << std::endl;
}

void curve_configure_k233(EllipticCurve *curve) {
	curve->field_size_bytes = 30;
	curve->binary_degree = 233;

	const unsigned char name[] = "K-233";
	for (int i = 0; i < 6; ++i)
		curve->curve_name_ascii[i] = name[i];

	for (unsigned long i = 0; i < curve->field_size_bytes; ++i) {
		curve->a[i] = 0x00;
		curve->b[i] = (i == 0) ? 0x01 : 0x00;
	}

	for (unsigned long i = 0; i < curve->field_size_bytes; ++i)
		curve->modulus[i] = 0x00;
	curve->modulus[0] |= 0x01;
	curve->modulus[74 / 8] |= (1U << (74 % 8));
	curve->modulus[233 / 8] |= (1U << (233 % 8));

	// Base Point G.X
	unsigned char xG[30] = { 0x26, 0x61, 0xAD, 0xEF, 0x6E, 0x9D, 0x4C, 0x0A,
			0xF5, 0x6B, 0xC2, 0x19, 0xA4, 0x63, 0x95, 0x14, 0xF4, 0x2F, 0xF2,
			0x29, 0xF1, 0x1A, 0x73, 0x7E, 0x3A, 0x85, 0xBA, 0x32, 0x72, 0x01 };

	// Base Point G.Y
	unsigned char yG[30] = { 0xA3, 0xE6, 0xFA, 0x56, 0x10, 0xC1, 0xE0, 0x56,
			0x9B, 0xEB, 0x8A, 0xF1, 0x9B, 0xCD, 0xA8, 0x27, 0xC4, 0x67, 0x5A,
			0x55, 0x0F, 0xF7, 0xB7, 0x19, 0xE8, 0xEC, 0x7D, 0x53, 0xDB, 0x01 };

	for (unsigned long i = 0; i < curve->field_size_bytes; ++i) {
		curve->xG[i] = xG[i];
		curve->yG[i] = yG[i];
	}

	unsigned char order_n[30] = { 0xDF, 0xAB, 0x73, 0xF1, 0xD5, 0x1A, 0xFB,
			0x6E, 0xD4, 0xBC, 0x15, 0xB9, 0x5B, 0x9D, 0x06, 0x00, 0x00, 0x00,
			0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80,
			0x00 };

	for (unsigned long i = 0; i < curve->field_size_bytes; ++i)
		curve->order[i] = order_n[i];

	curve->cofactor[0] = 0x04;
}
