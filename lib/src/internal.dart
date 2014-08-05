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

/**
 * enable debugging and assertions
 */
bool NDEBUG = true;

void ASSERT(bool a) {
  if (!NDEBUG) {
    assert(a);
  }
}

void ASSERT_INT(int a) {
  ASSERT(a != 0);
}

/**
 * @return true if an integer (stored in double x) would overflow (or if
 * x is NaN)
 */
bool INT_OVERFLOW(double x) {
  return ((!(x * (1.0 + 1e-8) <= INT_MAX)) || SCALAR_IS_NAN(x));
}

final int TRUE = 1;
final int FALSE = 0;

num MAX(num a, num b) {
  return a > b ? a : b;
}

num MIN(num a, num b) {
  return a < b ? a : b;
}

/* FLIP is a "negation about -1", and is used to mark an integer i that is
 * normally non-negative.  FLIP (EMPTY) is EMPTY.  FLIP of a number > EMPTY
 * is negative, and FLIP of a number < EMTPY is positive.  FLIP (FLIP (i)) = i
 * for all integers i.  UNFLIP (i) is >= EMPTY. */
const int EMPTY = -1;
const double EMPTY_D = -1.0;

num FLIP(num i) {
  return -i - 2;
}

num UNFLIP(num i) {
  return (i < EMPTY) ? FLIP(i) : i;
}
