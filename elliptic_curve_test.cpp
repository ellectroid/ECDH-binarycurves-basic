#include <iostream>
#include <iomanip>
#include "elliptic_curve.h"


int test_elliptic_curves() {

    // Now test the NIST K-233 curve (sect233k1)
    std::cout << "\n--- Testing multiplying base point by base subgroup order. Expected result: (0;0) ---\n";

    std::cout << "\n--- Testing K-233 (sect233k1) ---\n";

	const unsigned long field_len = 30;
	const unsigned long y_offset = (field_len + 7UL) & (~7UL);

	EllipticCurve k233_curve = { 0 };
	EllipticCurvePoint point_A = { 0 };
	EllipticCurvePoint point_B = { 0 };
	EllipticCurvePoint point_C = { 0 };

	// Metadata
	k233_curve.field_size_bytes = field_len;
	k233_curve.binary_degree = 233;

	const unsigned char name[] = "K-233";
	for (int i = 0; i < 6; ++i)
		k233_curve.curve_name_ascii[i] = name[i];

	// Coefficients: a = 0, b = 1
	for (unsigned long i = 0; i < field_len; ++i) {
		k233_curve.a[i] = 0x00;
		k233_curve.b[i] = (i == 0) ? 0x01 : 0x00;
	}

	// Modulus: x^233 + x^74 + 1
	for (unsigned long i = 0; i < field_len; ++i)
		k233_curve.modulus[i] = 0x00;
	k233_curve.modulus[0] |= 0x01;
	k233_curve.modulus[74 / 8] |= (1U << (74 % 8));
	k233_curve.modulus[233 / 8] |= (1U << (233 % 8));

	// Base Point G.X
	unsigned char xG[30] = { 0x26, 0x61, 0xAD, 0xEF, 0x6E, 0x9D, 0x4C, 0x0A,
			0xF5, 0x6B, 0xC2, 0x19, 0xA4, 0x63, 0x95, 0x14, 0xF4, 0x2F, 0xF2,
			0x29, 0xF1, 0x1A, 0x73, 0x7E, 0x3A, 0x85, 0xBA, 0x32, 0x72, 0x01 };

	// Base Point G.Y
	unsigned char yG[30] = { 0xA3, 0xE6, 0xFA, 0x56, 0x10, 0xC1, 0xE0, 0x56,
			0x9B, 0xEB, 0x8A, 0xF1, 0x9B, 0xCD, 0xA8, 0x27, 0xC4, 0x67, 0x5A,
			0x55, 0x0F, 0xF7, 0xB7, 0x19, 0xE8, 0xEC, 0x7D, 0x53, 0xDB, 0x01 };

	for (unsigned long i = 0; i < field_len; ++i) {
		k233_curve.xG[i] = xG[i];
		k233_curve.yG[i] = yG[i];
		point_A.point_mem[i] = xG[i];
		point_A.point_mem[i + y_offset] = yG[i];
	}

	// Order of G (n)
	unsigned char order_n[30] = { 0xDF, 0xAB, 0x73, 0xF1, 0xD5, 0x1A, 0xFB,
			0x6E, 0xD4, 0xBC, 0x15, 0xB9, 0x5B, 0x9D, 0x06, 0x00, 0x00, 0x00,
			0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80,
			0x00 };



	for (unsigned long i = 0; i < field_len; ++i)
		k233_curve.order[i] = order_n[i];

	k233_curve.cofactor[0] = 0x04;

	unsigned char multiplication_exponent[GF2_VECTOR_MAX_BYTELEN] = { };
	multiplication_exponent[0] = 0x2;




	std::cout << "Computing for K-233...\n";
	elliptic_curve_binary_point_multiply(&k233_curve, &point_C, &point_A,
			order_n, field_len);
	std::cout << "Result:\nX: ";
	for (unsigned long i = 0; i < field_len; ++i)
		std::cout << std::hex << std::setw(2) << std::setfill('0')
				<< (int) point_C.point_mem[i];
	std::cout << "\nY: ";
	for (unsigned long i = 0; i < field_len; ++i)
		std::cout << std::hex << std::setw(2) << std::setfill('0')
				<< (int) point_C.point_mem[i + y_offset];
	std::cout << std::endl;


	// ------------------------------------------------------------
	// Now test the NIST K-163 curve (sect163k1)
	std::cout << "\n--- Testing K-163 (sect163k1) ---\n";

	const unsigned long field_len2 = 21;
	const unsigned long y_offset2  = (field_len2 + 7UL) & (~7UL);

	// curve & points
	EllipticCurve     k163_curve = { 0 };
	EllipticCurvePoint Pk         = { 0 };  // G
	EllipticCurvePoint Rk         = { 0 };  // order·G

	// metadata
	k163_curve.field_size_bytes = field_len2;
	k163_curve.binary_degree    = 163;
	const unsigned char name2[] = "K-163";
	for (int i = 0; i < 5; ++i)
	    k163_curve.curve_name_ascii[i] = name2[i];

	// a = 1, b = 1
	alignas(8) unsigned char a2[field_len2] = { 1 };
	alignas(8) unsigned char b2[field_len2] = { 1 };
	for (unsigned long i = 0; i < field_len2; ++i) {
	    k163_curve.a[i] = a2[i];
	    k163_curve.b[i] = b2[i];
	}

	// modulus polynomial x^163 + x^7 + x^6 + x^3 + 1
	// LSB‐first uint32 words: {0xC9,0,0,0,0,0x08}
	alignas(8) unsigned char poly2[field_len2] = {
	    0xC9,0x00,0x00,0x00,  // word0 = 0x000000C9
	    0x00,0x00,0x00,0x00,  // word1
	    0x00,0x00,0x00,0x00,  // word2
	    0x00,0x00,0x00,0x00,  // word3
	    0x00,0x00,0x00,0x00,  // word4
	    0x08                   // first byte of word5 = 0x00000008
	};
	for (unsigned long i = 0; i < field_len2; ++i)
	    k163_curve.modulus[i] = poly2[i];

	// base point G (LSB‐first bytes from given uint32 words)
	alignas(8) unsigned char xG2[field_len2] = {
	    0xE8,0xEE,0x94,0x5C,
	    0x5E,0x6D,0x4E,0xDE,
	    0x93,0xD7,0x07,0xAA,
	    0xAC,0x11,0xBC,0x7B,
	    0x53,0xC0,0x13,0xFE,
	    0x02
	};
	alignas(8) unsigned char yG2[field_len2] = {
	    0xD9,0xA3,0xDA,0xCC,
	    0x38,0xD5,0x36,0x05,
	    0x80,0x2E,0x1F,0x32,
	    0x58,0xFF,0x38,0x5D,
	    0xB0,0x0F,0x07,0x89,
	    0x02
	};
	for (unsigned long i = 0; i < field_len2; ++i) {
	    k163_curve.xG[i] = xG2[i];
	    k163_curve.yG[i] = yG2[i];
	    Pk.point_mem[i]         = xG2[i];
	    Pk.point_mem[i + y_offset2] = yG2[i];
	}

	// order of G (LSB‐first bytes from given uint32 words)
	alignas(8) unsigned char ord2[field_len2] = {
	    0xEF,0xA5,0xF8,0x99,
	    0x0D,0xCC,0xE0,0xA2,
	    0x08,0x01,0x02,0x00,
	    0x00,0x00,0x00,0x00,
	    0x00,0x00,0x00,0x00,
	    0x04
	};
	for (unsigned long i = 0; i < field_len2; ++i)
	    k163_curve.order[i] = ord2[i];

	k163_curve.cofactor[0] = 0x02;

	std::cout << "Computing for K-163...\n";
	elliptic_curve_binary_point_multiply(
	    &k163_curve,
	    &Rk,
	    &Pk,
		ord2,
	    field_len2
	);

	std::cout << "Result:\nX: ";
	for (unsigned long i = 0; i < field_len2; ++i)
	    std::cout << std::hex << std::setw(2) << std::setfill('0')
	              << (int)Rk.point_mem[i];
	std::cout << "\nY: ";
	for (unsigned long i = 0; i < field_len2; ++i)
	    std::cout << std::hex << std::setw(2) << std::setfill('0')
	              << (int)Rk.point_mem[i + y_offset2];
	std::cout << std::endl;
	// ------------------------------------------------------------






    // ------------------------------------------------------------
    // Now test the NIST B-163 curve (sect163r1)
    std::cout << "\n--- Testing B-163 (sect163r1) ---\n";

    const unsigned long field_len3 = 21;
    const unsigned long y_offset3  = (field_len3 + 7UL) & (~7UL);

    EllipticCurve     b163_curve = { 0 };
    EllipticCurvePoint Pb         = { 0 };  // G
    EllipticCurvePoint Rb         = { 0 };  // order·G

    // metadata
    b163_curve.field_size_bytes = field_len3;
    b163_curve.binary_degree    = 163;
    const unsigned char name3[] = "B-163";
    for (int i = 0; i < 5; ++i)
        b163_curve.curve_name_ascii[i] = name3[i];

    // Coefficient a = 1
    for (unsigned long i = 0; i < field_len3; ++i)
        b163_curve.a[i] = (i == 0 ? 0x01 : 0x00);

    // Coefficient b (LSB-first bytes of {0x4a3205fd,0x512f7874,0x1481eb10,0xb8c953ca,0x0a601907,0x00000002})
    alignas(8) unsigned char bcoef[field_len3] = {
        0xFD,0x05,0x32,0x4A,
        0x74,0x78,0x2F,0x51,
        0x10,0xEB,0x81,0x14,
        0xCA,0x53,0xC9,0xB8,
        0x07,0x19,0x60,0x0A,
        0x02
    };
    for (unsigned long i = 0; i < field_len3; ++i)
        b163_curve.b[i] = bcoef[i];

    // Modulus x^163 + x^7 + x^6 + x^3 + 1 (same as K-163)
    alignas(8) unsigned char poly3[field_len3] = {
        0xC9,0x00,0x00,0x00,
        0x00,0x00,0x00,0x00,
        0x00,0x00,0x00,0x00,
        0x00,0x00,0x00,0x00,
        0x00,0x00,0x00,0x00,
        0x08
    };
    for (unsigned long i = 0; i < field_len3; ++i)
        b163_curve.modulus[i] = poly3[i];

    // Base point G (LSB-first bytes of base_x, base_y)
    alignas(8) unsigned char xB3[field_len3] = {
        0x36,0x3E,0x34,0xE8,
        0x37,0x46,0x99,0xD4,
        0x68,0x11,0x99,0xA0,
        0x7E,0xD5,0xA2,0x86,
        0x62,0xA1,0xEB,0xF0,
        0x03
    };
    alignas(8) unsigned char yB3[field_len3] = {
        0xF1,0x24,0x73,0x79,
        0x0C,0x5C,0x1C,0xB1,
        0x45,0xD5,0xCD,0xA2,
        0x4F,0x09,0xA0,0x71,
        0x6C,0xBC,0x1F,0xD5,
        0x00
    };
    for (unsigned long i = 0; i < field_len3; ++i) {
        b163_curve.xG[i] = xB3[i];
        b163_curve.yG[i] = yB3[i];
        Pb.point_mem   [i]         = xB3[i];
        Pb.point_mem   [i + y_offset3] = yB3[i];
    }

    alignas(8) unsigned char ord3[field_len3] = {
        0x33,0x4C,0x23,0xA4,
        0x12,0x0C,0xE7,0x77,
        0xFE,0x92,0x02,0x00,
        0x00,0x00,0x00,0x00,
        0x00,0x00,0x00,0x00,
        0x04
    };
    for (unsigned long i = 0; i < field_len3; ++i)
        b163_curve.order[i] = ord3[i];

    b163_curve.cofactor[0] = 0x02;

	std::cout << "Computing for B-163...\n";
    elliptic_curve_binary_point_multiply(
        &b163_curve,
        &Rb,
        &Pb,
        ord3,
        field_len3
    );

	std::cout << "Result:\nX: ";
    for (unsigned long i = 0; i < field_len3; ++i)
        std::cout << std::hex << std::setw(2) << std::setfill('0')
                  << (int)Rb.point_mem[i];
    std::cout << "\nY: ";
    for (unsigned long i = 0; i < field_len3; ++i)
        std::cout << std::hex << std::setw(2) << std::setfill('0')
                  << (int)Rb.point_mem[i + y_offset3];
    std::cout << std::endl;
    // ------------------------------------------------------------


    // ------------------------------------------------------------
    // Now test the NIST K-283 curve (sect283k1)
    std::cout << "\n--- Testing K-283 (sect283k1) ---\n";

    const unsigned long field_len4 = 36;
    const unsigned long y_offset4  = (field_len4 + 7UL) & (~7UL);

    EllipticCurve     k283_curve = { 0 };
    EllipticCurvePoint P283       = { 0 };  // G
    EllipticCurvePoint R283       = { 0 };  // order·G

    // metadata
    k283_curve.field_size_bytes = field_len4;
    k283_curve.binary_degree    = 283;
    const unsigned char name4[] = "K-283";
    for (int i = 0; i < 5; ++i)
        k283_curve.curve_name_ascii[i] = name4[i];

    // a = 0
    for (unsigned long i = 0; i < field_len4; ++i)
        k283_curve.a[i] = 0x00;

    // b = 1
    alignas(8) unsigned char b283[field_len4] = {
        0x01,0x00,0x00,0x00,   0x00,0x00,0x00,0x00,
        0x00,0x00,0x00,0x00,   0x00,0x00,0x00,0x00,
        0x00,0x00,0x00,0x00,   0x00,0x00,0x00,0x00,
        0x00,0x00,0x00,0x00,   0x00,0x00,0x00,0x00,
        0x00,0x00,0x00,0x00
    };
    for (unsigned long i = 0; i < field_len4; ++i)
        k283_curve.b[i] = b283[i];

    // modulus x^283 + x^12 + x^7 + x^5 + 1
    // represented as LSB-first uint32 words {0x000010A1,...,0x08000000}
    alignas(8) unsigned char poly283[field_len4] = {
        0xA1,0x10,0x00,0x00,   0x00,0x00,0x00,0x00,
        0x00,0x00,0x00,0x00,   0x00,0x00,0x00,0x00,
        0x00,0x00,0x00,0x00,   0x00,0x00,0x00,0x00,
        0x00,0x00,0x00,0x00,   0x00,0x00,0x00,0x00,
        0x00,0x00,0x00,0x08
    };
    for (unsigned long i = 0; i < field_len4; ++i)
        k283_curve.modulus[i] = poly283[i];

    // base point G
    alignas(8) unsigned char x283[field_len4] = {
        0x36,0x28,0x49,0x58,   0x24,0xAC,0xC2,0xB0,
        0x13,0x69,0x87,0x16,   0x7A,0x56,0xC1,0x23,
        0x5F,0x26,0xCD,0x53,   0xE5,0x88,0xF1,0x62,
        0x81,0x3B,0x1A,0x3F,   0x88,0x44,0xCA,0x78,
        0x3F,0x21,0x03,0x05
    };
    alignas(8) unsigned char y283[field_len4] = {
        0x59,0x22,0xDD,0x77,   0x61,0x11,0x34,0x4E,
        0x36,0x62,0x59,0xE4,   0x98,0x46,0x18,0xE8,
        0xC0,0x45,0x7E,0xE8,   0x6F,0x42,0xE5,0x07,
        0x5D,0xF9,0x90,0x8D,   0x31,0x9E,0x1C,0x0F,
        0x38,0xDA,0xCC,0x01
    };
    for (unsigned long i = 0; i < field_len4; ++i) {
        k283_curve.xG[i] = x283[i];
        k283_curve.yG[i] = y283[i];
        P283.point_mem[i]         = x283[i];
        P283.point_mem[i + y_offset4] = y283[i];
    }

    // order of G
    alignas(8) unsigned char ord283[field_len4] = {
        0x61,0x3C,0x16,0x1E,   0x06,0x1E,0x45,0x94,
        0x7F,0xFF,0x5D,0x26,   0x77,0x75,0xD0,0x2E,
        0xAE,0xE9,0xFF,0xFF,   0xFF,0xFF,0xFF,0xFF,
        0xFF,0xFF,0xFF,0xFF,   0xFF,0xFF,0xFF,0xFF,
        0xFF,0xFF,0xFF,0x01
    };
    for (unsigned long i = 0; i < field_len4; ++i)
        k283_curve.order[i] = ord283[i];

    k283_curve.cofactor[0] = 0x04;

	std::cout << "Computing for K-283...\n";
    elliptic_curve_binary_point_multiply(
        &k283_curve,
        &R283,
        &P283,
		ord283,
        field_len4
    );

	std::cout << "Result:\nX: ";
    for (unsigned long i = 0; i < field_len4; ++i)
        std::cout << std::hex << std::setw(2) << std::setfill('0')
                  << (int)R283.point_mem[i];
    std::cout << "\nY: ";
    for (unsigned long i = 0; i < field_len4; ++i)
        std::cout << std::hex << std::setw(2) << std::setfill('0')
                  << (int)R283.point_mem[i + y_offset4];
    std::cout << std::endl;

    std::cout << "Testing if base point G is on the curve: " << +elliptic_curve_binary_point_on_curve(&k283_curve,
    		&P283) << std::endl;
    // ------------------------------------------------------------

	return 0;
}
