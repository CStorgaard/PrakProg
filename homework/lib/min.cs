using System;

public static class NewtonMinimizer
{
    // Public static property to store the number of iterations performed.
    public static int IterationCount { get; private set; }

    /// <summary>
    /// Minimizes the function φ using Newton's method with numerical derivatives,
    /// with an option to use central finite difference approximations, and backtracking line-search.
    /// </summary>
    /// <param name="phi">Objective function: maps a vector to a scalar.</param>
    /// <param name="x0">Initial guess vector.</param>
    /// <param name="tol">Tolerance for the gradient norm.</param>
    /// <param name="maxIter">Maximum number of Newton iterations.</param>
    /// <param name="central">If true, use central differences; otherwise, use forward differences.</param>
    /// <returns>The minimizer x found by the algorithm.</returns>
    public static Vec Minimize(Func<Vec, double> phi, Vec x0, double tol = 1e-3, int maxIter = 50, bool central = false)
    {
        // Reset the iteration counter.
        IterationCount = 0;
        
        // Make a copy of the initial guess.
        Vec x = x0.Copy();
        double phiVal = phi(x);

        while (IterationCount < maxIter)
        {
            Vec grad;
            if (central)
            {
                // Compute the numerical gradient using central differences.
                grad = CentralGradient(phi, x);
            }
            else
            {
                // Compute the numerical gradient using forward differences.
                grad = Gradient(phi, x);
            }
            if (grad.Norm() < tol)
                break;

            Mat H;
            if (central)
            {
                // Compute the numerical Hessian using central differences.
                H = CentralHessian(phi, x);
            }
            else
            {
                // Compute the numerical Hessian using forward differences.
                H = Hessian(phi, x);
            }

            // Regularize H by adding a small value to the diagonal.
            double reg = 1e-8;
            for (int i = 0; i < H.Rows; i++)
            {
                H[i, i] += reg;
            }

            // Solve the Newton system: H * dx = -grad.
            Vec dx = H.SolveQR(-grad);

            // Backtracking line-search.
            double lambda = 1.0;
            while (lambda > 1.0 / 128.0)
            {
                Vec trial = x + lambda * dx;
                if (phi(trial) < phiVal)
                    break;
                lambda /= 2.0;
            }

            // Update the iterate.
            x = x + lambda * dx;
            phiVal = phi(x);
            IterationCount++;
        }

        return x;
    }

    /// <summary>
    /// Computes the numerical gradient of φ at x using finite forward differences.
    /// f'(x) ≈ [f(x+δx)-f(x)]/δx
    /// </summary>
    public static Vec Gradient(Func<Vec, double> phi, Vec x)
    {
        int n = x.Length;
        Vec grad = new Vec(n);
        double phi_x = phi(x);
        double factor = Math.Pow(2, -26);
        
        for (int i = 0; i < n; i++)
        {
            double dxi = Math.Abs(x[i]) * factor;
            if (dxi == 0)
                dxi = factor;

            x[i] += dxi;
            double phiPerturbed = phi(x);
            grad[i] = (phiPerturbed - phi_x) / dxi;
            x[i] -= dxi;  // Restore original value.
        }
        return grad;
    }

    /// <summary>
    /// Computes the numerical Hessian of φ at x using finite forward differences.
    /// f''(x) ≈ [f'(x+δx)-f'(x)]/δx  (applied componentwise)
    /// </summary>
    public static Mat Hessian(Func<Vec, double> phi, Vec x)
    {
        int n = x.Length;
        Mat H = new Mat(n, n);
        Vec gphi = Gradient(phi, x);
        double factor = Math.Pow(2, -13);

        for (int j = 0; j < n; j++)
        {
            double dxj = Math.Abs(x[j]) * factor;
            if (dxj == 0)
                dxj = factor;

            x[j] += dxj;
            Vec gphiPerturbed = Gradient(phi, x);
            Vec dgphi = gphiPerturbed - gphi;
            for (int i = 0; i < n; i++)
            {
                H[i, j] = dgphi[i] / dxj;
            }
            x[j] -= dxj;  // Restore original value.
        }
        return H;
    }

    /// <summary>
    /// Computes the numerical gradient of φ at x using central finite differences.
    /// f'(x) ≈ [f(x+δx)-f(x-δx)]/[2δx]
    /// </summary>
    public static Vec CentralGradient(Func<Vec, double> phi, Vec x)
    {
        int n = x.Length;
        Vec grad = new Vec(n);
        double factor = Math.Pow(2, -13);
        for (int i = 0; i < n; i++)
        {
            double original = x[i];
            double delta = Math.Abs(original) * factor;
            if (delta == 0)
                delta = factor;

            x[i] = original + delta;
            double f_plus = phi(x);
            x[i] = original - delta;
            double f_minus = phi(x);
            grad[i] = (f_plus - f_minus) / (2 * delta);
            x[i] = original;
        }
        return grad;
    }

    /// <summary>
    /// Computes the numerical Hessian of φ at x using central finite differences.
    /// Diagonal: f''(x) ≈ [f(x+δx)-2f(x)+f(x-δx)]/δx²
    /// Off-diagonals: ∂²f/∂x_i∂x_j ≈ [f(x+δ_i+δ_j)-f(x+δ_i-δ_j)-f(x-δ_i+δ_j)+f(x-δ_i-δ_j)]/(4δ_iδ_j)
    /// </summary>
    public static Mat CentralHessian(Func<Vec, double> phi, Vec x)
    {
        int n = x.Length;
        Mat H = new Mat(n, n);
        double factor = Math.Pow(2, -13);
        double f0 = phi(x);

        // Diagonal elements.
        for (int i = 0; i < n; i++)
        {
            double original = x[i];
            double delta = Math.Abs(original) * factor;
            if (delta == 0)
                delta = factor;
            x[i] = original + delta;
            double f_plus = phi(x);
            x[i] = original - delta;
            double f_minus = phi(x);
            H[i, i] = (f_plus - 2 * f0 + f_minus) / (delta * delta);
            x[i] = original;
        }

        // Off-diagonals.
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                double orig_i = x[i], orig_j = x[j];
                double delta_i = Math.Abs(orig_i) * factor;
                if (delta_i == 0)
                    delta_i = factor;
                double delta_j = Math.Abs(orig_j) * factor;
                if (delta_j == 0)
                    delta_j = factor;

                // Evaluate f(x + δ_i e_i + δ_j e_j)
                x[i] = orig_i + delta_i;
                x[j] = orig_j + delta_j;
                double f_pp = phi(x);
                
                // Evaluate f(x + δ_i e_i - δ_j e_j)
                x[j] = orig_j - delta_j;
                double f_pm = phi(x);
                
                // Evaluate f(x - δ_i e_i + δ_j e_j)
                x[i] = orig_i - delta_i;
                x[j] = orig_j + delta_j;
                double f_mp = phi(x);
                
                // Evaluate f(x - δ_i e_i - δ_j e_j)
                x[j] = orig_j - delta_j;
                double f_mm = phi(x);
                
                double value = (f_pp - f_pm - f_mp + f_mm) / (4 * delta_i * delta_j);
                H[i, j] = value;
                H[j, i] = value;
                x[i] = orig_i;
                x[j] = orig_j;
            }
        }
        return H;
    }
}