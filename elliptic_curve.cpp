#include "elliptic_curve.h"
#include "galois_field2.h"

unsigned long elliptic_curve_get_maximum_vector_bytelen() {
	return GF2_VECTOR_MAX_BYTELEN;
}

void elliptic_curve_binary_point_add(EllipticCurve *curve,
		EllipticCurvePoint *out, EllipticCurvePoint *in1,
		EllipticCurvePoint *in2) {

	unsigned long len = curve->field_size_bytes;
	unsigned long dbl_len = 2 * len;
	unsigned long y_offset = (len + 7UL) & (~7UL);
	unsigned long zero_count = 0;
	unsigned long zero_count2 = 0;

	alignas(8) unsigned char x1[GF2_VECTOR_MAX_BYTELEN] = { };
	alignas(8) unsigned char y1[GF2_VECTOR_MAX_BYTELEN] = { };
	alignas(8) unsigned char x2[GF2_VECTOR_MAX_BYTELEN] = { }; // reused for x1 + x3
	alignas(8) unsigned char y2[GF2_VECTOR_MAX_BYTELEN] = { }; // reused for lambda inv

	alignas(8) unsigned char x3[GF2_VECTOR_MAX_BYTELEN] = { };
	alignas(8) unsigned char y3[GF2_VECTOR_MAX_BYTELEN] = { };

	alignas(8) unsigned char lambda[2 * GF2_VECTOR_MAX_BYTELEN] = { }; // holds lambda and lambda sq.
	alignas(8) unsigned char temp[2 * GF2_VECTOR_MAX_BYTELEN] = { }; // holds lambda numerator, denominator, intermediate

	for (unsigned long i = 0; i < len; ++i) {
		x1[i] = in1->point_mem[i];
		y1[i] = in1->point_mem[i + y_offset];
		x2[i] = in2->point_mem[i];
		y2[i] = in2->point_mem[i + y_offset];
	}

	for (unsigned long i = 0; i < len; ++i) {
		zero_count += (x1[i] != x2[i]);
		zero_count += (y1[i] != y2[i]);
	}
	if (zero_count == 0) {
		elliptic_curve_binary_point_double(curve, out, in1);
		return;
	}

	zero_count = 0;
	for (unsigned long i = 0; i < len; ++i) {
		zero_count += (x1[i] != 0);
		zero_count += (y1[i] != 0);
	}
	if (zero_count == 0) {
		//operand 1 is zero
		//copy operand 2 to output and return
		for (unsigned long i = 0; i < len; ++i) {
			out->point_mem[i] = in2->point_mem[i];
			out->point_mem[i + y_offset] = in2->point_mem[i + y_offset];
		}
		return;
	}
	zero_count = 0;
	for (unsigned long i = 0; i < len; ++i) {
		zero_count += (x2[i] != 0);
		zero_count += (y2[i] != 0);
	}
	if (zero_count == 0) {
		//operand 1 is zero
		//copy operand 2 to output and return
		for (unsigned long i = 0; i < len; ++i) {
			out->point_mem[i] = in1->point_mem[i];
			out->point_mem[i + y_offset] = in1->point_mem[i + y_offset];
		}
		return;
	}

	//Inverse-pair test (P + -P = inf
	zero_count = 0;
	zero_count2 = 0;
	for (unsigned long i = 0; i < len; ++i) {
		zero_count += (x1[i] != x2[i]);
	}
	for (unsigned long i = 0; i < len; ++i) {
		zero_count2 += (y2[i] != (x1[i] ^ y1[i]));
	}
	if((zero_count == 0) && (zero_count2 == 0)){
		for (unsigned long i = 0; i < len; ++i) {
					out->point_mem[i] = 0x00;
					out->point_mem[i + y_offset] = 0x00;
		}
		return;
	}



	for (unsigned long i = 0; i < len; ++i) {
		temp[i] = y1[i] ^ y2[i];                   // numerator
		temp[i + y_offset] = x1[i] ^ x2[i];        // denominator
	}

	gf2_binary_inverse_lsb(&temp[y_offset], y2, len, curve->modulus);
	gf2_multiply_lsb(temp, y2, lambda, len);
	gf2_reduce_lsb(lambda, dbl_len, curve->modulus, len);

	gf2_multiply_lsb(lambda, lambda, temp, len);
	gf2_reduce_lsb(temp, dbl_len, curve->modulus, len);

	for (unsigned long i = 0; i < len; ++i) {
		x3[i] = temp[i] ^ lambda[i] ^ x1[i] ^ x2[i] ^ curve->a[i];
		x2[i] = x1[i] ^ x3[i];

	}
	gf2_multiply_lsb(lambda, x2, temp, len);
	gf2_reduce_lsb(temp, dbl_len, curve->modulus, len);

	for (unsigned long i = 0; i < len; ++i)
		y3[i] = temp[i] ^ x3[i] ^ y1[i];

	for (unsigned long i = 0; i < len; ++i) {
		out->point_mem[i] = x3[i];
		out->point_mem[i + y_offset] = y3[i];
	}
}


void elliptic_curve_binary_point_double(EllipticCurve *curve,
		EllipticCurvePoint *out, EllipticCurvePoint *in) {
	unsigned long len = curve->field_size_bytes; //byte len of a gf(2) vector
	unsigned long y_offset = (len + 7UL) & (~7UL); //placement of y coordinate in curve object
	unsigned long dbl_len = 2 * len; //length of vector after multiplication
	unsigned long zero_cnt = 0; //helper

	alignas(8) unsigned char x1[GF2_VECTOR_MAX_BYTELEN] = { };
	alignas(8) unsigned char y1[GF2_VECTOR_MAX_BYTELEN] = { };
	alignas(8) unsigned char lambda[2 * GF2_VECTOR_MAX_BYTELEN] = { };
	alignas(8) unsigned char temp1[2 * GF2_VECTOR_MAX_BYTELEN] = { };
	alignas(8) unsigned char temp2[2 * GF2_VECTOR_MAX_BYTELEN] = { };

	for (unsigned long i = 0; i < len; ++i) {
		if (in->point_mem[i] == 0x00)
			zero_cnt++;
	}
	if (zero_cnt == len) {
		for (unsigned long i = 0; i < len; ++i) {
			out->point_mem[i] = x1[i];
			out->point_mem[i + y_offset] = y1[i];
		}
		return; //point at infinity
	}

	//Following formula from NIST SP 800-186
	for (unsigned long i = 0; i < len; ++i) {
		x1[i] = in->point_mem[i];
		y1[i] = in->point_mem[i + y_offset];
	}
	gf2_binary_inverse_lsb(x1, temp1, len, curve->modulus);
	gf2_multiply_lsb(y1, temp1, lambda, len);
	gf2_reduce_lsb(lambda, dbl_len, curve->modulus, len);
	for (unsigned long i = 0; i < len; ++i)
		lambda[i] ^= x1[i];
	gf2_multiply_lsb(x1, x1, temp1, len);
	gf2_reduce_lsb(temp1, dbl_len, curve->modulus, len);
	gf2_multiply_lsb(lambda, lambda, temp2, len);
	gf2_reduce_lsb(temp2, dbl_len, curve->modulus, len);
	for (unsigned long i = 0; i < len; ++i)
		x1[i] = temp2[i] ^ lambda[i] ^ curve->a[i];   //x1 contains coordinate X of the output

	for (unsigned long i = 0; i < len; ++i)
		temp1[i] = x1[i] ^ in->point_mem[i];
	gf2_multiply_lsb(lambda, temp1, temp2, len);
	gf2_reduce_lsb(temp2, dbl_len, curve->modulus, len);
	for (unsigned long i = 0; i < len; ++i)
		y1[i] = temp2[i] ^ x1[i] ^ in->point_mem[i + y_offset]; //y1 contains coordinate Y of the output

	for (unsigned long i = 0; i < len; ++i) {
		out->point_mem[i] = x1[i];
		out->point_mem[i + y_offset] = y1[i];
	}
}

void elliptic_curve_binary_point_multiply(EllipticCurve *curve, EllipticCurvePoint *out,
		EllipticCurvePoint *in, const unsigned char *exp,
		unsigned long bytelen) {
	unsigned long len = curve->field_size_bytes;
	unsigned long y_offset = (len + 7UL) & (~7UL);
	long total_bits = gf2_degree_lsb(exp, bytelen);

	alignas(8) EllipticCurvePoint tmp = { };
	alignas(8) EllipticCurvePoint dummy = { };

	for (unsigned long i = 0; i < 2 * GF2_VECTOR_MAX_BYTELEN; ++i) {
		tmp.point_mem[i] = 0x00;
		dummy.point_mem[i] = 0x00;
	}

	for (long i = (long) (total_bits); i >= 0; --i) {
		elliptic_curve_binary_point_double(curve, &tmp, &tmp);
		unsigned long byte_index = (unsigned long) i >> 3;
		unsigned long bit_index = (unsigned long) i & 7;
		unsigned char bit = (exp[byte_index] >> bit_index) & 1;
		if (bit) {
			elliptic_curve_binary_point_add(curve, &tmp, &tmp, in);
		} else {
			elliptic_curve_binary_point_add(curve, &tmp, &dummy, &tmp);
		}
	}

	for (unsigned long i = 0; i < len; ++i) {
		out->point_mem[i] = tmp.point_mem[i];
		out->point_mem[i + y_offset] = tmp.point_mem[i + y_offset];
	}
}

int elliptic_curve_binary_point_on_curve(EllipticCurve *curve,
                                         EllipticCurvePoint *point)
{
    unsigned long len      = curve->field_size_bytes;
    unsigned long y_offset = (len + 7UL) & (~7UL);
    unsigned long dbl_len  = 2 * len;
    unsigned char zero_flag = 0;

    alignas(8) unsigned char x[GF2_VECTOR_MAX_BYTELEN]     = {0};
    alignas(8) unsigned char y[GF2_VECTOR_MAX_BYTELEN]     = {0};
    alignas(8) unsigned char lhs[GF2_VECTOR_MAX_BYTELEN]   = {0};
    alignas(8) unsigned char rhs[GF2_VECTOR_MAX_BYTELEN]   = {0};
    alignas(8) unsigned char tmp[2 * GF2_VECTOR_MAX_BYTELEN] = {0};
    alignas(8) unsigned char x2[GF2_VECTOR_MAX_BYTELEN]    = {0};
    alignas(8) unsigned char x3[GF2_VECTOR_MAX_BYTELEN]    = {0};
    alignas(8) unsigned char ax2[GF2_VECTOR_MAX_BYTELEN]   = {0};

    for (unsigned long i = 0; i < len; i++) {
        x[i] = point->point_mem[i];
        y[i] = point->point_mem[i + y_offset];
    }

    for (unsigned long i = 0; i < len; i++)
        zero_flag |= x[i] | y[i];
    if (zero_flag == 0)
        return 1;

    gf2_multiply_lsb(y, y, tmp, len);
    gf2_reduce_lsb(tmp, dbl_len, curve->modulus, len);
    for (unsigned long i = 0; i < len; i++)
        lhs[i] = tmp[i];

    gf2_multiply_lsb(x, y, tmp, len);
    gf2_reduce_lsb(tmp, dbl_len, curve->modulus, len);
    for (unsigned long i = 0; i < len; i++)
        lhs[i] ^= tmp[i];

    gf2_multiply_lsb(x, x, tmp, len);
    gf2_reduce_lsb(tmp, dbl_len, curve->modulus, len);
    for (unsigned long i = 0; i < len; i++)
        x2[i] = tmp[i];

    gf2_multiply_lsb(x2, x, tmp, len);
    gf2_reduce_lsb(tmp, dbl_len, curve->modulus, len);
    for (unsigned long i = 0; i < len; i++)
        x3[i] = tmp[i];

    gf2_multiply_lsb(curve->a, x2, tmp, len);
    gf2_reduce_lsb(tmp, dbl_len, curve->modulus, len);
    for (unsigned long i = 0; i < len; i++)
        ax2[i] = tmp[i];

    for (unsigned long i = 0; i < len; i++)
        rhs[i] = x3[i] ^ ax2[i] ^ curve->b[i];

    unsigned char diff = 0;
    for (unsigned long i = 0; i < len; i++)
        diff |= lhs[i] ^ rhs[i];

    return (diff == 0);
}

unsigned char* elliptic_curve_point_get_coord_x(EllipticCurve *curve, EllipticCurvePoint *point){
	(void)curve;
	return &point->point_mem[0];
}
unsigned char* elliptic_curve_point_get_coord_y(EllipticCurve *curve, EllipticCurvePoint *point){
	return &point->point_mem[(curve->field_size_bytes + 7UL) & (~7UL)];
}
unsigned long elliptic_curve_point_get_coord_one_bytelen(EllipticCurve *curve){
	return (curve->field_size_bytes);
}
unsigned long elliptic_curve_point_get_coord_full_bytelen(EllipticCurve *curve){
	return (((curve->field_size_bytes + 7UL) & (~7UL)) + curve->field_size_bytes);
}
