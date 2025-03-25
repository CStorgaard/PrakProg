using static System.Math;
using System;

public class IntegrationResult {
    public double Value { get; set; }
    public double Error { get; set; }
}

/// <summary>
/// Provides an adaptive integration method based on an open quadrature rule,
/// returning both the integral and an error estimate.
/// </summary>

public class Integrator {
    // Static counter for function evaluations
    public static long EvaluationCount { get; private set; } = 0;

    // Resets the evaluation counter.
    public static void ResetEvaluationCount() {
        EvaluationCount = 0;
    }

    // Helper method: increments the counter and evaluates f at x.
    private static double Evaluate(Func<double, double> f, double x) {
        EvaluationCount++;
        return f(x);
    }

    public static IntegrationResult Integrate(
        Func<double, double> f, 
        double a, 
        double b,
        double absTolerance = 0.001, 
        double relTolerance = 0.001, 
        double f2 = double.NaN, 
        double f3 = double.NaN,
        int depth = 0,
        int maxDepth = 1000)
    {
        double h = b - a;
        
        // If maximum recursion depth is reached, compute the quadrature and return.
        if (depth >= maxDepth)
        {
            double f1 = Evaluate(f, a + h / 6);
            double f4 = Evaluate(f, a + 5 * h / 6);
            if (double.IsNaN(f2))
            {
                f2 = Evaluate(f, a + 2 * h / 6);
                f3 = Evaluate(f, a + 4 * h / 6);
            }
            double Q = (2 * f1 + f2 + f3 + 2 * f4) / 6 * h;
            double q = (f1 + f2 + f3 + f4) / 4 * h;
            return new IntegrationResult { Value = Q, Error = Abs(Q - q) };
        }
        
        // On first call, compute f2 and f3 if not provided.
        if (double.IsNaN(f2))
        {
            f2 = Evaluate(f, a + 2 * h / 6);
            f3 = Evaluate(f, a + 4 * h / 6);
        }
        
        // Evaluate additional points.
        double f1_new = Evaluate(f, a + h / 6);
        double f4_new = Evaluate(f, a + 5 * h / 6);
        
        // Compute the higher-order and lower-order approximations.
        double Q_high = (2 * f1_new + f2 + f3 + 2 * f4_new) / 6 * h;
        double Q_low  = (f1_new + f2 + f3 + f4_new) / 4 * h;
        double error = Abs(Q_high - Q_low);
        
        // If within tolerance, return the result.
        if (error <= absTolerance + relTolerance * Abs(Q_high))
            return new IntegrationResult { Value = Q_high, Error = error };
        else
        {
            double mid = (a + b) / 2;
            IntegrationResult left = Integrate(f, a, mid, absTolerance / Sqrt(2), relTolerance, f1_new, f2, depth + 1, maxDepth);
            IntegrationResult right = Integrate(f, mid, b, absTolerance / Sqrt(2), relTolerance, f3, f4_new, depth + 1, maxDepth);
            return new IntegrationResult { Value = left.Value + right.Value, Error = left.Error + right.Error };
        }
    }


    /// <summary>
    /// Computes the integral using the Clenshaw-Curtis variable transformation.
    /// </summary>
    public static IntegrationResult CCIntegrate(
        Func<double, double> f, 
        double a, 
        double b, 
        double absTol = 0.001, 
        double relTol = 0.001)
    {
        // Define the transformed integrand g(θ)
        // x = (a+b)/2 + (b-a)/2*cosθ
        // dx/dθ = -(b-a)/2 sinθ, but reversing the limits gives a positive factor.
        Func<double, double> g = theta => f((a + b) / 2 + (b - a) / 2 * Cos(theta)) 
                                        * Sin(theta) * (b - a) / 2;
        // Now integrate g from 0 to π using the existing adaptive integrator.
        return Integrate(g, 0, PI, absTol, relTol);
    }
}
