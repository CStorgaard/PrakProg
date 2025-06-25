using System;
using System.Diagnostics;
using BerrutInterpolation;

class ProgramTime
{
    static void Main()
    {
        // different sizes to test O(n)
        int[] ns = {100, 200, 400, 800, 1600, 3200};
        const int M = 1000;   // number of eval calls per n
        var rnd = new Random(12345); // fixed seed for reproducibility

        // Header for Out_times.txt
        Console.WriteLine("# n  avg_eval_time_ns");

        foreach (int n in ns)
        {
            // 1) build equispaced nodes+values on [0,2π]
            double[] xs = new double[n];
            double[] ys = new double[n];
            for (int i = 0; i < n; i++)
            {
                xs[i] = 2*Math.PI * i/(n-1); // [0, 2π]
                ys[i] = xs[i] + 0.2*Math.Sin(5*xs[i]); // Function to interpolate
            }

            var interp = new BerrutInterpolator(xs, ys, linearReproduction: true); // Perform interpolation with linear reproduction

            // warm-up (JIT, cache, etc.)
            interp.Evaluate(xs[n/2]);

            // 2) time M evaluations at random points
            var sw = Stopwatch.StartNew();
            double sink = 0.0; 
            for (int k = 0; k < M; k++) 
            {
                double t = rnd.NextDouble() * 2*Math.PI;
                sink += interp.Evaluate(t);
            }
            sw.Stop();

            // average time in miliseconds
            double avgNs = sw.Elapsed.TotalMilliseconds * 1e6 / M;
            Console.WriteLine($"{n} {avgNs:F2}");
        }   
    }
}
