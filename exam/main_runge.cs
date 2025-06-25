using System;
using BerrutInterpolation;

class Program
{
    static void Main()
    {
        // 11 equispaced nodes on [-1,1]
        int nNodes = 11;
        double[] xs = new double[nNodes];
        double[] ys = new double[nNodes];
        for (int i = 0; i < nNodes; i++)
        {
            xs[i] = -1.0 + 2.0*i/(nNodes-1);
            ys[i] = 1.0/(1.0 + 25.0*xs[i]*xs[i]);  // Function to interpolate
        }

        var r1 = new BerrutInterpolator(xs, ys, linearReproduction:false);
        var r2 = new BerrutInterpolator(xs, ys, linearReproduction:true);

        // === block 0: nodes ===
        for (int i = 0; i < nNodes; i++)
            Console.WriteLine($"{xs[i]} {ys[i]}");
        Console.WriteLine(); Console.WriteLine();

        // === block 1: interpolants + true f on a fine grid ===
        int M = 200;
        for (int k = 0; k <= M; k++)
        {
            double x = -1.0 + 2.0*k/M;
            double f = 1.0/(1.0 + 25.0*x*x);
            Console.WriteLine($"{x} {r1.Evaluate(x)} {r2.Evaluate(x)} {f}");
        }
    }
}