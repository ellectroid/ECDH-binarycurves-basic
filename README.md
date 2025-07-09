# ECDH-binarycurves-basic
A small library for doing elliptic curve Diffie-Hellman key exchange using binary curves.

## WARNING
This library is <b>not a constant time library</b>. It only gives algebraically correct results!

## Description
A tiny recreational elliptic curve Diffie-Hellman key exchange library.   
Dependencies: none. Doesn't use standard libraries, doesn't use any runtime libraries (even no memset or memcpy), doesn't #include anything at all.  
It should be fully C-compatible C++ code.  
Should work on any platform for which one can compile C or C++.   
Min language requirements: C11/C++11 (only because of alignas)   
There is room for optimization of underlying algorithms for better use of CPU time   
This library can work as a reference   

Tested with curves K-163, B-163, K-233, K-283 (Binary curves only! Prime curves not supported!)   

Tested on desktop   
Tested on ARM Cortex-M7 STM32F746 @ 180MHz:
  - ~1.35KiB peak RAM usage with keys and calculations with 37 bytes allocated per GF(2) vector (e.g. with K-283)
  - ~1.15KiB peak RAM usage with keys and calculations with 32 bytes allocated per GF(2) vector (e.g. with K-233)
  - ~6.5KiB ROM use for debug build
  - ~2.9KiB ROM use for release build
  - Calculating a public key from a private key took around 25 seconds using curve K-233 (curves greater than K-233 seem impractical for it)
 
## Usage
The library doesn't provide curve parameters, but this repository includes library test code that contains parameters for K-163, B-163, K-233 and K-283.   
Make sure to adjust compile-time value GF2_VECTOR_MAX_BYTELEN in galois_field2.h to make sure it fits your GF(2) vector. Use one byte more than required to fit the entire vector.
By default it's 32, so it won't fit K-283 without adjustment.  
