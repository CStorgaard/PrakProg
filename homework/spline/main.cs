using System;
using System.IO;
using static System.Math;

class main
{
    static void Main()
    {
        // Sample data points
        double[] xs = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        double[] ys = {Cos(0), Cos(1), Cos(2), Cos(3), Cos(4), Cos(5), Cos(6), Cos(7), Cos(8), Cos(9), Cos(10)};

        // Create a spline instance (using cubic spline in this example)
        Spline linearSpline = new Spline(xs, ys, SplineType.Linear);

        // Write original data points to "dataPoints.txt"
        using (StreamWriter writer = new StreamWriter("dataPoints.txt"))
        {
            writer.WriteLine("# x  y");
            for (int i = 0; i < xs.Length; i++)
            {
                writer.WriteLine($"{xs[i]} {ys[i]}");
            }
        }

        // Evaluate spline over a fine grid and write results to "splineResults.txt"
        int steps = 100;
        double xMin = xs[0];
        double xMax = xs[xs.Length - 1];
        double dx = (xMax - xMin) / (steps - 1);

        using (StreamWriter writer = new StreamWriter("linearSplineResults.txt"))
        {
            writer.WriteLine("# x  S(x)  S'(x)  Integral");
            for (int i = 0; i < steps; i++)
            {
                double xVal = xMin + i * dx;
                double splineVal = linearSpline.Evaluate(xVal);
                double splineDeriv = linearSpline.Derivative(xVal);
                double splineInt = linearSpline.Integral(xVal);
                writer.WriteLine($"{xVal} {splineVal} {splineDeriv} {splineInt}");
            }
        }
        
        Console.WriteLine("Data points written to dataPoints.txt");
        Console.WriteLine("Spline results written to linearSplineResults.txt");

        // Make a new quadratic spline and write results to "quadraticSplineResults.txt"
        Spline quadraticSpline = new Spline(xs, ys, SplineType.Quadratic);
        using (StreamWriter writer = new StreamWriter("quadraticSplineResults.txt"))
        {
            writer.WriteLine("# x  S(x)  S'(x)  Integral");
            for (int i = 0; i < steps; i++)
            {
                double xVal = xMin + i * dx;
                double splineVal = quadraticSpline.Evaluate(xVal);
                double splineDeriv = quadraticSpline.Derivative(xVal);
                double splineInt = quadraticSpline.Integral(xVal);
                writer.WriteLine($"{xVal} {splineVal} {splineDeriv} {splineInt}");
            }
        }
        Console.WriteLine("Quadratic spline results written to quadraticSplineResults.txt");

        // Make a cubic spline and write results to "cubicSplineResults.txt"
        Spline cubicSpline = new Spline(xs, ys, SplineType.Cubic);
        using (StreamWriter writer = new StreamWriter("cubicSplineResults.txt"))
        {
            writer.WriteLine("# x  S(x)  S'(x)  Integral");
            for (int i = 0; i < steps; i++)
            {
                double xVal = xMin + i * dx;
                double splineVal = cubicSpline.Evaluate(xVal);
                double splineDeriv = cubicSpline.Derivative(xVal);
                double splineInt = cubicSpline.Integral(xVal);
                writer.WriteLine($"{xVal} {splineVal} {splineDeriv} {splineInt}");
            }
        }
        Console.WriteLine("Cubic spline results written to cubicSplineResults.txt");
    }
}
