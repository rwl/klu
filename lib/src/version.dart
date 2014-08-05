/**
 * KLU: a sparse LU factorization algorithm.
 * Copyright (C) 2004-2009, Timothy A. Davis.
 * Copyright (C) 2011-2012, Richard W. Lincoln.
 * http://www.cise.ufl.edu/research/sparse/klu
 *
 * -------------------------------------------------------------------------
 *
 * KLU is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * KLU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this Module; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */

part of edu.ufl.cise.klu.tdouble;

const int KLU_OK = 0;
/** status > 0 is a warning, not an error */
const int KLU_SINGULAR = 1;
const int KLU_OUT_OF_MEMORY = -2;
const int KLU_INVALID = -3;
/** integer overflow has occured */
const int KLU_TOO_LARGE = -4;

/** enable diagnostic printing */
bool NPRINT = true;

const int INT_MAX = 0x7fffffff;

Float64List GET_I_POINTER(Float64List LU, Int32List Xip, int Xip_offset, Int32List Xi_offset, int k) {
  Xi_offset[0] = Xip[Xip_offset + k];
  return LU;
}

Float64List GET_POINTER(Float64List LU, Int32List Xip, int Xip_offset, Int32List Xlen, int Xlen_offset, Int32List Xi_offset, Int32List Xx_offset, int k, Int32List xlen) {
  int xp = Xip[Xip_offset + k];
  xlen[0] = Xlen[Xlen_offset + k];
  Xi_offset[0] = xp;
  Xx_offset[0] = xp + xlen[0];
  return LU;
}

bool SCALAR_IS_NAN(double x) {
  return x != x;
}

bool SCALAR_IS_ZERO(double x) {
  return x == 0.0;
}

bool SCALAR_IS_NONZERO(double x) {
  return x != 0.0;
}

bool SCALAR_IS_LTZERO(double x) {
  return x < 0.0;
}

/* scalar absolute value macro. If x is NaN, the result is NaN */
double SCALAR_ABS(double x) {
  return SCALAR_IS_LTZERO(x) ? -x : x;
}

void PRINTF(String format) {
  if (!NPRINT) {
    print(format);
  }
}

void PRINT_SCALAR(double a) {
  if (!NPRINT) {

    if (SCALAR_IS_NONZERO(a)) {
      PRINTF(" ($a)");
    } else {
      PRINTF(" (0)");
    }

  }
}

/* ---------------------------------------------------------------------- */
/* Real floating-point arithmetic */
/* ---------------------------------------------------------------------- */

/**
 * Returns [TRUE] if a complex number is in split form, [FALSE] if in packed
 * form.
 */
int SPLIT(double s) {
  return 1;
}

/**
 * @return True if a is NaN
 */
bool IS_NAN(double a) {
  return SCALAR_IS_NAN(a);
}

/**
 * @return True if a == 0
 */
bool IS_ZERO(double a) {
  return SCALAR_IS_ZERO(a);
}

/**
 * @return True if a != 0
 */
bool IS_NONZERO(double a) {
  return SCALAR_IS_NONZERO(a);
}

/**
 * c /= s
 */
double SCALE_DIV(double c, double s) {
  c /= s;
  return c;
}

/**
 * c -= a*b
 */
double MULT_SUB(double c, double a, double b) {
  c -= a * b;
  return c;
}

void PRINT_ENTRY(double a) {
  PRINT_SCALAR(a);
}

double ABS(double a) {
  return a.abs();
}

void CLEAR(Float64List A, int i) {
  A[i] = 0.0;
}

/* for flop counts */
const double MULTSUB_FLOPS = 2.0; // c -= a*b
const double DIV_FLOPS = 1.0; // c = a/b
const double ABS_FLOPS = 0.0; // c = abs(a)
const double ASSEMBLE_FLOPS = 1.0; // c += a
const double DECREMENT_FLOPS = 1.0; // c -= a
const double MULT_FLOPS = 1.0; // c = a*b
const double SCALE_FLOPS = 1.0; // c = a/s
