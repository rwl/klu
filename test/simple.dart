//package edu.ufl.cise.klu.test;

/*
import junit.framework.TestCase;
import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_numeric;
import edu.ufl.cise.klu.common.KLU_symbolic;

import static edu.ufl.cise.klu.tdouble.Dklu_defaults.klu_defaults;
import static edu.ufl.cise.klu.tdouble.Dklu_analyze.klu_analyze;
import static edu.ufl.cise.klu.tdouble.Dklu_factor.klu_factor;
import static edu.ufl.cise.klu.tdouble.Dklu_solve.klu_solve;
*/

//public class Dklu_simple extends TestCase {

static final double DELTA = 1e-09 ;

static int n = 5 ;
static List<int> Ap = [0, 2, 5, 9, 10, 12] ;
static List<int> Ai = [0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4] ;
static List<double> Ax = [2.0, 3.0, 3.0, -1.0, 4.0, 4.0, -3.0, 1.0, 2.0, 2.0, 6.0, 1.0] ;
static List<double> b = [8.0, 45.0, -3.0, 3.0, 19.] ;

/**
 * a simple KLU demo; solution is x = (1,2,3,4,5)
 */
main() {
  int i;
  KLU_symbolic Symbolic;
  KLU_numeric Numeric;
  KLU_common Common = new KLU_common();

  //Dklu_version.NPRINT = false ;
  //Dklu_internal.NDEBUG = false ;

  klu_defaults (Common);
  Symbolic = klu_analyze (n, Ap, Ai, Common);
  Numeric = klu_factor (Ap, Ai, Ax, Symbolic, Common);
  klu_solve (Symbolic, Numeric, 5, 1, b, 0, Common);

  for (i = 0 ; i < n ; i++) {
    System.out.printf("x [%d] = %g\n", i, b [i]) ;
    assertEquals(i + 1.0, b [i], DELTA) ;
  }
}
