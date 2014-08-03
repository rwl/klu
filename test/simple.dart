import 'dart:io';
import 'package:unittest/unittest.dart';
import 'package:klu/klu.dart' as klu;
import 'package:klu/common/common.dart';

const double DELTA = 1e-09 ;

final int n = 5 ;
final List<int> Ap = [0, 2, 5, 9, 10, 12] ;
final List<int> Ai = [0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4] ;
final List<double> Ax = [2.0, 3.0, 3.0, -1.0, 4.0, 4.0, -3.0, 1.0, 2.0, 2.0, 6.0, 1.0] ;
final List<double> b = [8.0, 45.0, -3.0, 3.0, 19.0] ;

main() {
  test('A simple KLU demo; solution is x = (1,2,3,4,5)', () {
    final Common = new KLU_common();

    //Dklu_version.NPRINT = false ;
    //Dklu_internal.NDEBUG = false ;

    klu.defaults (Common);
    final Symbolic = klu.analyze (n, Ap, Ai, Common);
    final Numeric = klu.factor (Ap, Ai, Ax, Symbolic, Common);
    klu.solve (Symbolic, Numeric, 5, 1, b, 0, Common);

    for (int i = 0 ; i < n ; i++) {
      stdout.write("x [$i] = ${b [i]}\n") ;
      expect(i + 1.0, closeTo(b [i], DELTA)) ;
    }
  });
}
