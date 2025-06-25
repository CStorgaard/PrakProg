using System;
using BerrutInterpolation;

class ProgramError
{
    static void Main()
    {
        // 9 equispaced nodes in [0,2π]
        double[] xs = {
            0.0, Math.PI/4, Math.PI/2, 3*Math.PI/4,
            Math.PI, 5*Math.PI/4, 3*Math.PI/2, 7*Math.PI/4, 2*Math.PI
        };
        double[] ys = new double[xs.Length];
        for(int i=0;i<xs.Length;i++)
            ys[i] = xs[i] + 0.2*Math.Sin(5*xs[i]);

        var r1 = new BerrutInterpolator(xs, ys, linearReproduction:false);
        var r2 = new BerrutInterpolator(xs, ys, linearReproduction:true);

        // header
        Console.WriteLine("# x    f(x)    r1(x)   r2(x)   err1    err2");

        int M = 200;
        double a = xs[0], b = xs[xs.Length-1];
        for(int k=0; k<=M; k++)
        {
            double x = a + (b-a)*k/M;
            double f = x + 0.2*Math.Sin(5*x);
            double v1 = r1.Evaluate(x);
            double v2 = r2.Evaluate(x);
            double e1 = Math.Abs(v1 - f);
            double e2 = Math.Abs(v2 - f);
            // eliminate near-zero roundoff so we get exact zeros at the nodes
            if (e1 < 1e-10) e1 = 0.0;
            if (e2 < 1e-10) e2 = 0.0;
            Console.WriteLine($"{x:F4} {f:F6} {v1:F6} {v2:F6} {e1:E3} {e2:E3}");
        }
    }
}
