/**
 * Read in a matrix and solve a linear system.
 */
import 'dart:io';
import 'dart:typed_data';
import 'package:unittest/unittest.dart';
import 'package:klu/klu.dart' as klu;
import 'package:csparse/test_util.dart';
import 'package:csparse/csparse.dart';

final String DIR = "matrix";
final String ARROW = "arrow";
final String IMPCOL_A = "impcol_a";
final String WEST0156 = "west0156";
//final String C1 = "1c";
//final String ARROW_C = "arrowc";
//final String C_TINA = "ctina";
//final String GD99_CC = "GD99_cc";
//final String ONE = "one";
//final String ONE_C = "onec";
//final String TWO = "two";
//final String W156 = "w156";

File get_file(String name) {
  return new File([Uri.base.toFilePath() + DIR, name].join('/'));
}

void REAL(Float64List X, int i, double v) {
  X[2 * i] = v;
}

void IMAG(Float64List X, int i, double v) {
  X[2 * i + 1] = v;
}

double MAX(double a, double b) {
  return (a) > (b) ? (a) : (b);
}

/**
 *
 * @param n A is n-by-n
 * @param Ap size n+1, column pointers
 * @param Ai size nz = Ap [n], row indices
 * @param Ax size nz, numerical values
 * @param isreal nonzero if A is real, 0 otherwise
 * @param B size n, right-hand-side
 * @param X size n, solution to Ax=b
 * @param R size n, residual r = b-A*x
 * @param lunz size 1, nnz(L+U+F)
 * @param rnorm size 1, norm(b-A*x,1) / norm(A,1)
 * @param Common default parameters and statistics
 * @return 1 if successful, 0 otherwise
 */
int backslash(int n, Int32List Ap, Int32List Ai, Float64List Ax, bool isreal, Float64List B, Float64List X, Float64List R, Int32List lunz, Float64List rnorm, klu.KLU_common Common) {
  double anorm = 0.0,
      asum;
  klu.KLU_symbolic Symbolic;
  klu.KLU_numeric Numeric;
  int i, j, p;

  if (Ap == null || Ai == null || Ax == null || B == null || X == null || B == null) return (0);

  /* ---------------------------------------------------------------------- */
  /* symbolic ordering and analysis */
  /* ---------------------------------------------------------------------- */

  Symbolic = klu.analyze(n, Ap, Ai, Common);
  if (Symbolic == null) return (0);

  if (isreal) {

    /* ------------------------------------------------------------------ */
    /* factorization */
    /* ------------------------------------------------------------------ */

    Numeric = klu.factor(Ap, Ai, Ax, Symbolic, Common);
    if (Numeric == null) {
      //klu_free_symbolic(Symbolic, Common);
      return (0);
    }

    /* ------------------------------------------------------------------ */
    /* statistics(not required to solve Ax=b) */
    /* ------------------------------------------------------------------ */

    klu.rgrowth(Ap, Ai, Ax, Symbolic, Numeric, Common);
    klu.condest(Ap, Ax, Symbolic, Numeric, Common);
    klu.rcond(Symbolic, Numeric, Common);
    klu.flops(Symbolic, Numeric, Common);
    lunz[0] = Numeric.lnz + Numeric.unz - n + (Numeric.Offp != null ? Numeric.Offp[n] : 0);

    /* ------------------------------------------------------------------ */
    /* solve Ax=b */
    /* ------------------------------------------------------------------ */

    for (i = 0; i < n; i++) {
      X[i] = B[i];
    }
    klu.solve(Symbolic, Numeric, n, 1, X, 0, Common);

    /* ------------------------------------------------------------------ */
    /* compute residual, rnorm = norm(b-Ax,1) / norm(A,1) */
    /* ------------------------------------------------------------------ */

    for (i = 0; i < n; i++) {
      R[i] = B[i];
    }
    for (j = 0; j < n; j++) {
      asum = 0.0;
      for (p = Ap[j]; p < Ap[j + 1]; p++) {
        /* R(i) -= A(i,j) * X(j) */
        R[Ai[p]] -= Ax[p] * X[j];
        asum += Ax[p].abs();
      }
      anorm = MAX(anorm, asum);
    }
    rnorm[0] = 0.0;
    for (i = 0; i < n; i++) {
      rnorm[0] = MAX(rnorm[0], R[i].abs());
    }

    /* ------------------------------------------------------------------ */
    /* free numeric factorization */
    /* ------------------------------------------------------------------ */

    //klu_free_numeric(Numeric, Common);
    Numeric = null;

  } else {
    throw new Error();

    /* ------------------------------------------------------------------ */
    /* statistics(not required to solve Ax=b) */
    /* ------------------------------------------------------------------ */

    //Numeric = klu.z_factor (Ap, Ai, Ax, Symbolic, Common);
    //if (Numeric == null)
    //{
    //	klu_free_symbolic (Symbolic, Common);
    //	return(0);
    //}
    //
    ///* ------------------------------------------------------------------ */
    ///* statistics */
    ///* ------------------------------------------------------------------ */
    //
    //klu_z_rgrowth (Ap, Ai, Ax, Symbolic, Numeric, Common);
    //klu_z_condest (Ap, Ax, Symbolic, Numeric, Common);
    //klu_z_rcond (Symbolic, Numeric, Common);
    //klu_z_flops (Symbolic, Numeric, Common);
    //lunz = Numeric.lnz + Numeric.unz - n +
    //	(Numeric.Offp != null ? Numeric.Offp [n] : 0);
    //
    ///* ------------------------------------------------------------------ */
    ///* solve Ax=b */
    ///* ------------------------------------------------------------------ */
    //
    //for(i = 0; i < 2*n; i++)
    //{
    //	X [i] = B [i];
    //}
    //klu_z_solve (Symbolic, Numeric, n, 1, X, Common);
    //
    ///* ------------------------------------------------------------------ */
    ///* compute residual, rnorm = norm(b-Ax,1) / norm(A,1) */
    ///* ------------------------------------------------------------------ */
    //
    //for(i = 0; i < 2*n; i++)
    //{
    //	R [i] = B [i];
    //}
    //for(j = 0; j < n; j++)
    //{
    //	asum = 0;
    //	for(p = Ap [j]; p < Ap [j+1]; p++)
    //	{
    //		/* R(i) -= A(i,j) * X(j) */
    //		i = Ai [p];
    //		REAL(R,i) -= REAL(Ax,p) * REAL(X,j) - IMAG(Ax,p) * IMAG(X,j);
    //		IMAG(R,i) -= IMAG(Ax,p) * REAL(X,j) + REAL(Ax,p) * IMAG(X,j);
    //		asum += CABS(Ax, p);
    //	}
    //	anorm = MAX(anorm, asum);
    //}
    //rnorm = 0;
    //for(i = 0; i < n; i++)
    //{
    //	rnorm = MAX (rnorm, CABS(R, i));
    //}
    //
    ///* ------------------------------------------------------------------ */
    ///* free numeric factorization */
    ///* ------------------------------------------------------------------ */
    //
    //klu_z_free_numeric (&Numeric, Common);
  }

  /* ---------------------------------------------------------------------- */
  /* free symbolic analysis, and residual */
  /* ---------------------------------------------------------------------- */

  //klu_free_symbolic (Symbolic, Common);
  Symbolic = null;

  return (1);
}

/**
 * Given a sparse matrix A, set up a right-hand-side and solve X = A\b.
 */
void demo(int n, Int32List Ap, Int32List Ai, Float64List Ax, bool isreal, Int32List lunz, Float64List rnorm, klu.KLU_common Common) {
  int i;
  Float64List B, X, R;

  stdout.write("KLU: ${klu.KLU_DATE}, version: ${klu.KLU_MAIN_VERSION}.${klu.KLU_SUB_VERSION}.${klu.KLU_SUBSUB_VERSION}\n");

  /* ---------------------------------------------------------------------- */
  /* set defaults */
  /* ---------------------------------------------------------------------- */

  klu.defaults(Common);

  /* ---------------------------------------------------------------------- */
  /* create a right-hand-side */
  /* ---------------------------------------------------------------------- */

  if (isreal) {
    /* B = 1 +(1:n)/n */
    B = new Float64List(n);
    X = new Float64List(n);
    R = new Float64List(n);
    if (B != null) {
      for (i = 0; i < n; i++) {
        B[i] = 1 + (i + 1) / n;
      }
    }
  } else {
    /* real(B) = 1 +(1:n)/n, imag(B) = (n:-1:1)/n */
    B = new Float64List(2 * n);
    X = new Float64List(2 * n);
    R = new Float64List(2 * n);
    if (B != null) {
      for (i = 0; i < n; i++) {
        REAL(B, i, 1 + (i + 1) / n);
        IMAG(B, i, (n - i) / n);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /* X = A\b using KLU and print statistics */
  /* ---------------------------------------------------------------------- */

  if (backslash(n, Ap, Ai, Ax, isreal, B, X, R, lunz, rnorm, Common) == 0) {
    stdout.write("KLU failed\n");
  } else {
    stdout.write("n $n nnz(A) ${Ap [n]} nnz(L+U+F) ${lunz[0]} resid ${rnorm[0]}\n" +
        "recip growth ${Common.rgrowth} condest ${Common.condest} rcond ${Common.rcond} flops ${Common.flops}\n");
  }

  /* ---------------------------------------------------------------------- */
  /* free the problem */
  /* ---------------------------------------------------------------------- */

  B = null;
  X = null;
  R = null;

  stdout.write("peak memory usage: ${Common.mempeak} bytes\n\n");
}

main() {
  /**
   * n 207 nnz(A) 572 nnz(L+U+F) 615 resid 6.98492e-10
   * recip growth 0.00957447 condest 4.35093e+07 rcond 4.5277e-05 flops 259
   */
  test('impcol_a', () {
    //klu.NPRINT = false ;
    //klu.NDEBUG = false ;

    klu.KLU_common Common = new klu.KLU_common();
    Int32List lunz = new Int32List(1);
    Float64List rnorm = new Float64List(1);

    final file = get_file(IMPCOL_A);
    Dproblem prob = get_problem(file, 0.0, 1);
    Dcs A = prob.A;
    demo(A.m, A.p, A.i, A.x, true, lunz, rnorm, Common);

    expect(207, equals(A.m));
    expect(207, equals(A.n));
    expect(572, equals(A.p[A.m]));
    expect(615, equals(lunz[0]));
    expect(6.98492e-10, closeTo(rnorm[0], 1e-14));
    expect(0.00957447, closeTo(Common.rgrowth, 1e-06));
    expect(4.35093e+07, closeTo(Common.condest, 1e+06)); // FIXME: improve assertion accuracy
    expect(4.5277e-05, closeTo(Common.rcond, 1e-08));
    expect(259, closeTo(Common.flops, 1e-03));
  });

  /**
   * n 100 nnz(A) 298 nnz(L+U+F) 298 resid 1.77636e-15
   * recip growth 0.0204082 condest 303 rcond 0.0204082 flops 297
   */
  test('arrow', () {
    klu.KLU_common Common = new klu.KLU_common();
    Int32List lunz = new Int32List(1);
    Float64List rnorm = new Float64List(1);

    final file = get_file(ARROW);
    Dproblem prob = get_problem(file, 0.0, 1);
    Dcs A = prob.A;
    demo(A.m, A.p, A.i, A.x, true, lunz, rnorm, Common);

    expect(100, equals(A.m));
    expect(100, equals(A.n));
    expect(298, equals(A.p[A.m]));
    expect(298, equals(lunz[0]));
    expect(1.77636e-15, closeTo(rnorm[0], 1e-18));
    expect(0.0204082, closeTo(Common.rgrowth, 1e-06));
    expect(303, closeTo(Common.condest, 1e-03));
    expect(0.0204082, closeTo(Common.rcond, 1e-06));
    expect(297, closeTo(Common.flops, 1e-03));
  });

  /**
   * n 156 nnz(A) 371 nnz(L+U+F) 406 resid 1.04858e+06
   * recip growth 0.0306751 condest 1.64225e+31 rcond 9.48528e-08 flops 188
   */
  test('west0156', () {
    klu.KLU_common Common = new klu.KLU_common();
    Int32List lunz = new Int32List(1);
    Float64List rnorm = new Float64List(1);

    final file = get_file(WEST0156);
    Dproblem prob = get_problem(file, 0.0, 1);
    Dcs A = prob.A;
    demo(A.m, A.p, A.i, A.x, true, lunz, rnorm, Common);

    expect(156, equals(A.m));
    expect(156, equals(A.n));
    expect(371, equals(A.p[A.m]));
    expect(406, equals(lunz[0]));
    expect(1.04858e+06, closeTo(rnorm[0], 1e+02));
    expect(0.0306751, closeTo(Common.rgrowth, 1e-06));
    expect(1.64225e+31, closeTo(Common.condest, 1e+26));
    expect(9.48528e-08, closeTo(Common.rcond, 1e-12));
    expect(188, closeTo(Common.flops, 1e-03));
  });
}
