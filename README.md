# ersbch

Reed-Solomon Erasure Coding (BCH) in C++

ersbch.cpp = C++ BCH based erasure code.

ersbchc.cpp+ersbch32.asm - C++ and assembly erasure code.
ersbch32.asm functions use C++ mangled names.

Visual Studio 2015 was used to compile and/or assemble code.

In the example code, there are defines that setup a data matrix
as a 20 row x 32768 column matrix, with the first 17 rows as data,
and the last 3 rows as parities, for up to 3 erasure correction.

The example code is based on "BCH view" Reed Solomon code as an
alternative to "original view Vandermonde matrix" code.
"Original view" and "BCH view" are actually two different codes,
using different encoding and decoding algorithms, but both are
called Reed Solomon.

Almost all Reed Solomon error or error+erasure correcting codes
are BCH view. They were developed much earlier and are much faster
at correcting errors.

BCH view is also faster as an erasure only code.
Let p = # parity rows, and e = # erasure rows for correction.

BCH view encodes p-1 rows via matrix multiply and XOR to encode
the remaining row, essentially correcting the last row to encode it.

BCH view corrects e-1 erasure rows via matrix multiply, and XOR to
correct the remaining erasure row. A single erasure only uses XOR.

The Wikipedia article describes the differences between "orignal view"
and "BCH view", as well as a descrption of the encoders. The article
describes error decoders, but not erasure only decoders.

[Wiki Reed Solomon](https://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction)

In the example code, the field is GF(2^8) based on

x^8 + x^4 + x^3 + x^2 + 1 (hex 11d) with primitive element x + 0 (hex 2).
 
The example code polynomial roots are 1,2,4, which translates
into a generator polynomial g(x) = (x-1)(x-2)(x-4). In a binary field, both
add and subtract are XOR, so the generator polynomial can be described as:

g(x) = (x+1)(x+2)(x+4) = x^3 + 7 x^2 + 14 x + 8  (the 14 is decimal == 0x0e)

A custom matrix class based on an array of pointers to rows is used.
This matrix class includes an mapped matrix option, where the mapped
matrix points to all or a sub-set of the rows of another matrix.

Main sets up the matrices used:
```
  mEnc =     3 x 20 parity   (for encoding) matrix
  mSyn =     3 x 20 syndrome (for decoding) matrix
  mDat = 20 x 32768 encoded data matrix
                    17 data rows, 3 parity rows
  mPar =  3 x 32768 mapped matrix = parity rows of mDat
```

The function Patterns tests all 1, 2, and 3 erasure patterns.
The inputs to Patterns are mSyn and mDat.
InitCombination initializes for NextCombination.
NextCombination generate the erasure indexes for each erasure pattern.

For a single erasure, the example code XOR's rows of data for correction.

For n erasures, where n == 2 or n == 3, the following steps are performed:

1. mDat's erased rows are filled with 0xAA (representing garbage).
2. mFix = mapped matrix = erased rows of mDat.
2. mLct = n x n locator matrix is generated based on erasure indexes.
3. mInv = inverse of mLct.
4. mInv is reduced in size by one row.
4. mCor = mInv x mSyn (using (n-1) x 20 of mSyn). This is the correcting matrix.
5. The columns of mCor corresponding to erasures are set = 0.
6. mFix = mCor x mDat. This corrects n-1 rows of mDat.
7. The remaining erased row of mDat is corrected using XOR.
