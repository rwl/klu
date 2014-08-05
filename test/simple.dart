import 'dart:io';
import 'dart:typed_data';
import 'package:unittest/unittest.dart';
import 'package:klu/klu.dart' as klu;

const double DELTA = 1e-09;

final int n = 5;
final Ap = new Int32List.fromList([0, 2, 5, 9, 10, 12]);
final Ai = new Int32List.fromList([0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4]);
final Ax = new Float64List.fromList([2.0, 3.0, 3.0, -1.0, 4.0, 4.0, -3.0, 1.0, 2.0, 2.0, 6.0, 1.0]);
final b = new Float64List.fromList([8.0, 45.0, -3.0, 3.0, 19.0]);

main() {
  test('A simple KLU demo; solution is x = (1,2,3,4,5)', () {
    final Common = new klu.KLU_common();

    //klu.NPRINT = false ;
    //klu.NDEBUG = false ;

    klu.defaults(Common);
    final Symbolic = klu.analyze(n, Ap, Ai, Common);
    final Numeric = klu.factor(Ap, Ai, Ax, Symbolic, Common);
    klu.solve(Symbolic, Numeric, 5, 1, b, 0, Common);

    for (int i = 0; i < n; i++) {
      stdout.write("x [$i] = ${b [i]}\n");
      expect(i + 1.0, closeTo(b[i], DELTA));
    }
  });
}
