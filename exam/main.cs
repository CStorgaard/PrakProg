using System;
using BerrutInterpolation;

class Program
{
    static void Main()
    {
        // 9 equispaced nodes in [0,2π]
        double[] xs = {
            0.0,
            Math.PI/4,
            Math.PI/2,
            3*Math.PI/4,
            Math.PI,
            5*Math.PI/4,
            3*Math.PI/2,
            7*Math.PI/4,
            2*Math.PI
        };
        
        double[] ys = new double[xs.Length];
        for (int i = 0; i < xs.Length; i++)
        ys[i] = xs[i] + 0.2*Math.Sin(5*xs[i]);

        int n = xs.Length;

        // Build both interpolators
        var interp1 = new BerrutInterpolator(xs, ys, linearReproduction: false);
        var interp2 = new BerrutInterpolator(xs, ys, linearReproduction: true);

        // Block 0: the raw data points
        // columns: x   y
        for (int i = 0; i < n; i++)
        {
            Console.WriteLine($"{xs[i]} {ys[i]}");
        }

        Console.WriteLine();  // blank line = new data block
        Console.WriteLine();

        // Block 1: the interpolant curves
        // columns: x   r1(x)   r2(x)
        int M = 200;
        double x0 = xs[0], x1 = xs[n - 1];
        for (int i = 0; i <= M; i++)
        {
            double t = x0 + (x1 - x0) * i / M;
            double r1 = interp1.Evaluate(t);
            double r2 = interp2.Evaluate(t);
            Console.WriteLine($"{t} {r1} {r2}");
        }
    }
}