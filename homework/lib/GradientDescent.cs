using System;
using static System.Math;

/// <summary>
/// Simple gradient descent with momentum.
/// </summary>
public static class GradientDescent
{
    /// <summary>Number of iterations performed in the last call.</summary>
    public static int Iterations { get; private set; }

    /// <summary>
    /// Minimize the objective φ starting from x0, using momentum.
    /// </summary>
    /// <param name="phi">Objective φ(p) (for logging).</param>
    /// <param name="gradFunc">Function to compute ∇φ(p).</param>
    /// <param name="x0">Initial parameter vector.</param>
    /// <param name="step">Learning rate (fixed).</param>
    /// <param name="maxIter">Maximum number of iterations.</param>
    /// <param name="tol">Stop when ‖grad‖ < tol.</param>
    /// <param name="beta">Momentum coefficient (0 ≤ β < 1).</param>
    /// <returns>Approximate minimizer vector.</returns>
    public static Vec Minimize(
        Func<Vec,double> phi,
        Func<Vec,Vec> gradFunc,
        Vec x0,
        double step    = 1e-1,
        int    maxIter = 20000,
        double tol     = 1e-8,
        double beta    = 0.9
    ){
        Vec x = x0.Copy();
        Vec v = new Vec(x.Length);

        Iterations = 0;
        for(int iter = 1; iter <= maxIter; iter++){
            Vec g    = gradFunc(x);
            double n = g.Norm();
            Iterations = iter;
            if(n < tol) break;
            v = v*beta + g*(1 - beta);
            x = x - v*step;
        }
        return x;
    }
}