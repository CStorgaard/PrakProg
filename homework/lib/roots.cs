using System;

namespace NewtonSolver
{
    public static class Solver
    {
        /// <summary>
        /// Computes the Jacobian matrix of the function f at the point x using finite differences.
        /// If the dx vector is not provided, it is set to δx[i] = Max(|x[i]|, 1) * 2^(-26).
        /// </summary>
        /// <param name="f">Function from Vec to Vec.</param>
        /// <param name="x">Point at which to compute the Jacobian.</param>
        /// <param name="fx">Optional: f(x). If null, it will be computed.</param>
        /// <param name="dx">Optional: step sizes for the finite differences.</param>
        /// <returns>A Mat representing the Jacobian (with dimensions m×n if f returns an m-dimensional vector and x is n-dimensional).</returns>
        public static Mat Jacobian(Func<Vec, Vec> f, Vec x, Vec fx = null, Vec dx = null)
        {
            int n = x.Length;
            // Evaluate f(x) if not provided.
            if (fx == null)
                fx = f(x);

            int m = fx.Length;
            // If dx is not provided, choose δx[i] = Max(|x[i]|, 1) * 2^(-26)
            if (dx == null)
            {
                dx = new Vec(n);
                for (int i = 0; i < n; i++)
                {
                    double factor = Math.Abs(x[i]) > 1.0 ? Math.Abs(x[i]) : 1.0;
                    dx[i] = factor * Math.Pow(2, -26);
                }
            }

            // Jacobian will be an m x n matrix.
            Mat J = new Mat(m, n);
            // Compute finite difference for each variable.
            for (int j = 0; j < n; j++)
            {
                // Save original value.
                double temp = x[j];
                // Perturb x[j] by dx[j]
                x[j] += dx[j];
                Vec fPerturbed = f(x);
                // Compute the j-th column: (f(x+δx) - f(x)) / δx[j]
                for (int i = 0; i < m; i++)
                {
                    J[i, j] = (fPerturbed[i] - fx[i]) / dx[j];
                }
                // Restore x[j]
                x[j] = temp;
            }
            return J;
        }

        /// <summary>
        /// Applies Newton's method to find a root of the function f.
        /// Uses a backtracking line-search strategy to ensure sufficient decrease.
        /// Iterations stop when ||f(x)|| is below the desired accuracy or when the step-size becomes too small.
        /// </summary>
        /// <param name="f">Function from Vec to Vec whose root is sought.</param>
        /// <param name="start">Initial guess.</param>
        /// <param name="acc">Accuracy goal: algorithm stops when ||f(x)|| &lt; acc.</param>
        /// <param name="dx">Optional: δx vector for finite difference estimation. If null, defaults are used.</param>
        /// <returns>The computed root as a Vec.</returns>
        public static Vec Newton(Func<Vec, Vec> f, Vec start, double acc = 1e-2, Vec dx = null)
        {
            // Create a copy of the starting vector.
            Vec x = start.Copy();
            Vec fx = f(x);
            const double lambdaMin = 1e-6;  // Minimum allowable step size factor.
            int maxIter = 100;             // Maximum iterations to avoid infinite loops.
            int iter = 0;

            while (fx.Norm() >= acc && iter < maxIter)
            {
                // Compute the Jacobian matrix at x.
                Mat J = Jacobian(f, x, fx, dx);

                // Solve for the Newton step Dx in J * Dx = -fx.
                // Our Mat class has SolveQR which assumes a square system.
                Vec Dx = J.SolveQR(-fx);

                // Backtracking line-search.
                double lambda = 1.0;
                Vec z;
                Vec fz;
                while (true)
                {
                    // Compute trial point: z = x + lambda * Dx.
                    z = x + lambda * Dx;
                    fz = f(z);
                    // Check for sufficient decrease: ||f(z)|| < (1 - lambda/2) * ||f(x)||
                    if (fz.Norm() < (1 - lambda / 2) * fx.Norm())
                        break;
                    // If step-size becomes too small, break out.
                    if (lambda < lambdaMin)
                        break;
                    lambda /= 2;
                }
                // Update x and fx.
                x = z;
                fx = fz;
                iter++;
            }
            return x;
        }
    }
}
