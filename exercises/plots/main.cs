using System;
using static System.Math;
using System.Numerics;

class main {
    public static int Main() {
        for (double x = -3; x <= 3; x += 1.0 / 8) {
            System.Console.WriteLine($"{x} {sfuns.erf(x)}");
        }
        System.Console.WriteLine();
        System.Console.WriteLine();

        double dx = 1.0 / 64;
        for (double x = -5 + dx / 2; x <= 5; x += dx) {
            System.Console.WriteLine($"{x} {sfuns.gamma(x)}");
        }
        System.Console.WriteLine();
        System.Console.WriteLine();

        for (double x = dx; x <= 10; x += dx) {
            Console.WriteLine($"{x} {sfuns.lngamma(x).Real}");
        }
        System.Console.WriteLine();
        System.Console.WriteLine();

        double f = 1;
        for (int i = 1; i <= 4; i++) {
            f *= i;
            System.Console.WriteLine($"{i + 1} {f}");
        }
        System.Console.WriteLine();
        System.Console.WriteLine();

        double f1 = 1;
        for (int i = 1; i <= 9; i++) {
            f1 *= i;
            System.Console.WriteLine($"{i + 1} {Log(f1)}");
        }

        return 0;
    }
}
