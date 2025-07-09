#include <iostream>
#include <iomanip>

#include "galois_field2.h"

void test_gf2_4_inverses()
{
    const unsigned long BYTES = 1;          // GF(2^4) fits in one byte
    const unsigned char MOD = 0x13;         // x^4 + x + 1 = 0b0001_0011

    unsigned char a, inv, prod;
    unsigned char tmp[2];  // ← updated to 2 bytes

    bool all_ok = true;
    for (int ai = 1; ai < 16; ++ai) {
        a = (unsigned char)ai;

        // compute inverse in GF(2^4)
        inv = 0;
        gf2_binary_inverse_lsb(&a, &inv, BYTES, &MOD);

        // multiply a * inv → tmp[2]
        tmp[0] = tmp[1] = 0;
        gf2_multiply_lsb(&a, &inv, tmp, BYTES);
        // reduce modulo P(x)
        gf2_reduce_lsb(tmp, 2, &MOD, BYTES);
        prod = tmp[0];

        if (prod != 1) {
            std::cerr << "FAIL: a=0x"
                      << std::hex << std::uppercase << std::setw(1)
                      << ai
                      << " inverse=0x"
                      << std::setw(1)
                      << (int)inv
                      << "  a*inv mod P = 0x"
                      << (int)prod
                      << std::dec << "\n";
            all_ok = false;
        }
    }

    if (all_ok) {
        std::cout << "GF(2^4) inverse test passed for all a=1..15\n";
    } else {
        std::cout << "GF(2^4) inverse test FAILED\n";
    }
}

// Test GF(2^12) inversion using brute-force comparison
void test_gf2_12_inverse()
{
    const unsigned long BYTES = 2;
    const unsigned char MOD[2] = { 0x09, 0x10 };  // x^12 + x^3 + 1

    bool all_ok = true;
    unsigned short inv_table[4096] = {0};

    // Step 1: Brute-force inverse table
    for (int a = 1; a < 4096; ++a) {
        for (int b = 1; b < 4096; ++b) {
            unsigned char A[2] = { (unsigned char)(a & 0xFF), (unsigned char)(a >> 8) };
            unsigned char B[2] = { (unsigned char)(b & 0xFF), (unsigned char)(b >> 8) };
            unsigned char tmp[4] = {0};

            gf2_multiply_lsb(A, B, tmp, BYTES);
            gf2_reduce_lsb(tmp, 4, MOD, BYTES);

            if (tmp[0] == 1 && tmp[1] == 0) {
                inv_table[a] = (unsigned short)b;
                break;
            }
        }
    }

    // Step 2: Compare to Euclidean inverse
    for (int a = 1; a < 4096; ++a) {
        unsigned char A[2] = { (unsigned char)(a & 0xFF), (unsigned char)(a >> 8) };
        unsigned char out[64] = {0};
        gf2_binary_inverse_lsb(A, out, BYTES, MOD);

        unsigned short result = (unsigned short)(out[0] | (out[1] << 8));
        if (result != inv_table[a]) {
            all_ok = false;
            std::cout << "FAIL: a = " << a
                      << " → inverse = 0x" << std::hex << std::setw(3) << std::setfill('0') << result
                      << ", expected = 0x" << std::setw(3) << inv_table[a]
                      << std::dec << "\n";
        }
    }

    std::cout << (all_ok ? "GF(2^12) inverse test PASSED ✅\n"
                         : "GF(2^12) inverse test FAILED ❌\n");
}

void test_gf2_24_inverse_with_bruteforce() {
    const unsigned long BYTES = 3;
    const unsigned char MOD[4] = { 0x87, 0x00, 0x00, 0x01 }; // x^24 + x^7 + x^2 + x + 1

    // Test value — can randomize or rotate through a few
    unsigned char a[3] = { 0x4D, 0x23, 0x9A };

    // Compute inverse via Euclidean method
    unsigned char inv[64] = {0};
    gf2_binary_inverse_lsb(a, inv, BYTES, MOD);

    // Brute-force expected inverse
    unsigned char expected[3] = {0};
    bool found = false;
    for (unsigned int i = 1; i < (1 << 24); ++i) {
        unsigned char b[3] = {
            static_cast<unsigned char>(i & 0xFF),
            static_cast<unsigned char>((i >> 8) & 0xFF),
            static_cast<unsigned char>((i >> 16) & 0xFF)
        };
        unsigned char prod[6] = {0};
        gf2_multiply_lsb(a, b, prod, BYTES);
        gf2_reduce_lsb(prod, 6, MOD, 4);
        if (prod[0] == 1 && prod[1] == 0 && prod[2] == 0) {
            std::copy(b, b + 3, expected);
            found = true;
            break;
        }
        if (i % 1000000 == 0) std::cout << "." << std::flush; // progress indicator
    }

    // Multiply a * inv and reduce
    unsigned char result[6] = {0};
    gf2_multiply_lsb(a, inv, result, BYTES);
    gf2_reduce_lsb(result, 6, MOD, 4);

    // Output
    auto print_hex = [](const unsigned char* x, int len) {
        for (int i = len - 1; i >= 0; --i)
            std::cout << std::hex << std::uppercase
                      << std::setw(2) << std::setfill('0') << (int)x[i];
    };

    std::cout << "\n\nInput (a):            0x"; print_hex(a, 3); std::cout << "\n";
    std::cout << "Expected inverse:     ";
    if (found) { std::cout << "0x"; print_hex(expected, 3); }
    else       { std::cout << "(not found in search range)"; }
    std::cout << "\n";
    std::cout << "Computed inverse:     0x"; print_hex(inv, 3); std::cout << "\n";
    std::cout << "Product a * inv mod f: 0x"; print_hex(result, 3); std::cout << "\n";

    if (result[0] == 1 && result[1] == 0 && result[2] == 0 &&
        (!found || std::equal(inv, inv + 3, expected))) {
        std::cout << "✅ PASS: inverse correct and multiplicative identity holds\n";
    } else {
        std::cout << "❌ FAIL: mismatch or incorrect inverse\n";
    }
}

static void print_hex_lsb(const unsigned char* data, unsigned long bytelen);
void test_gf2_233_inverse() {
    const unsigned long FIELD_SIZE = 32;  // 256-bit buffer
    const unsigned long TEMP_SIZE  = 64;  // for full 512-bit stack workspace

    // xG value from NIST K-233 base point, LSB-first padded to 32 bytes
    unsigned char xG[FIELD_SIZE] = {
        0x26, 0x61, 0xAD, 0xEF, 0x6E, 0x9D, 0x4C, 0x0A,
        0xF5, 0x6B, 0xC2, 0x19, 0xA4, 0x63, 0x95, 0x14,
        0xF4, 0x2F, 0xF2, 0x29, 0xF1, 0x1A, 0x73, 0x7E,
        0x3A, 0x85, 0xBA, 0x32, 0x72, 0x01, 0x00, 0x00
    };

    // Modulus: x^233 + x^74 + 1
    unsigned char modulus[FIELD_SIZE] = {0};
    modulus[0]  |= 0x01; // x^0
    modulus[9]  |= 0x04; // x^74 → bit 2 of byte 9
    modulus[29] |= 0x02; // x^233 → bit 1 of byte 29

    unsigned char inverse[TEMP_SIZE] = {0};
    gf2_binary_inverse_lsb(xG, inverse, FIELD_SIZE, modulus);

    std::cout << "Inverse of xG in GF(2^233): ";
    print_hex_lsb(inverse, FIELD_SIZE);
}

static void print_hex_lsb(const unsigned char* data, unsigned long bytelen) {
    std::cout << "0x";
    for (long i = (long)bytelen - 1; i >= 0; --i) {
        unsigned char byte = data[i];
        const char* hex = "0123456789ABCDEF";
        std::cout << hex[(byte >> 4) & 0xF];
        std::cout << hex[byte & 0xF];
    }
    std::cout << std::endl;
}
