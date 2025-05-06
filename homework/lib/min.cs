using System;
using static System.Math;

public static class NewtonMinimizer
{
    /// <summary>Number of Newton iterations performed in the last call.</summary>
    public static int IterationCount { get; private set; }

    /// <summary>
    /// Minimize the objective φ starting from x0.
    /// Uses finite‐difference gradient (forward or central) and numeric Hessian.
    /// </summary>
    /// <param name="phi">Objective function mapping a Vec to a double.</param>
    /// <param name="x0">Initial guess.</param>
    /// <param name="tol">Tolerance on gradient norm for convergence.</param>
    /// <param name="maxIter">Maximum Newton iterations.</param>
    /// <param name="central">If true, use central differences; otherwise forward differences.</param>
    /// <returns>Approximate minimizer x*.</returns>
    public static Vec Minimize(
        Func<Vec, double> phi,
        Vec x0,
        double tol = 1e-3,
        int maxIter = 50,
        bool central = false
    ){
        IterationCount = 0;
        Vec x = x0.Copy();
        double f0 = phi(x);

        while(IterationCount < maxIter){
            // 1) Build gradient via finite differences
            Vec grad = central
                ? CentralGradient(phi, x)
                : Gradient(phi, x);

            // Check convergence
            if(grad.Norm() < tol) break;

            // 2) Build Hessian numerically (central or forward)
            Mat H = central
                ? CentralHessian(phi, x)
                : Hessian(phi, x);

            // Regularize diagonal for stability
            for(int i = 0; i < H.Rows; i++)
                H[i,i] += 1e-8;

            // 3) Solve Newton step H * dx = -grad
            Vec dx = H.SolveQR(grad * -1.0);

            // 4) Backtracking line‐search
            double lambda = 1.0;
            while(lambda > 1.0/128.0){
                Vec xNew = x + dx * lambda;
                if(phi(xNew) < f0) break;
                lambda *= 0.5;
            }

            // 5) Update iterate
            x = x + dx * lambda;
            f0 = phi(x);
            IterationCount++;
        }

        return x;
    }

    /// <summary>Forward-difference gradient.</summary>
    public static Vec Gradient(Func<Vec, double> phi, Vec x){
        int n = x.Length;
        Vec grad = new Vec(n);
        double fx = phi(x);
        double eps = Pow(2, -26);
        for(int i=0;i<n;i++){
            double d = Abs(x[i])*eps;
            if(d == 0) d = eps;
            x[i] += d;
            double fpx = phi(x);
            grad[i] = (fpx - fx)/d;
            x[i] -= d;
        }
        return grad;
    }

    /// <summary>Forward-difference Hessian.</summary>
    public static Mat Hessian(Func<Vec, double> phi, Vec x){
        int n = x.Length;
        Mat H = new Mat(n, n);
        Vec g0 = Gradient(phi, x);
        double eps = Pow(2, -13);
        for(int j=0;j<n;j++){
            double d = Abs(x[j])*eps;
            if(d == 0) d = eps;
            x[j] += d;
            Vec gj = Gradient(phi, x);
            Vec diff = gj - g0;
            for(int i=0;i<n;i++) H[i,j] = diff[i]/d;
            x[j] -= d;
        }
        return H;
    }

    /// <summary>Central-difference gradient.</summary>
    public static Vec CentralGradient(Func<Vec, double> phi, Vec x){
        int n = x.Length;
        Vec grad = new Vec(n);
        double eps = Pow(2, -13);
        for(int i=0;i<n;i++){
            double orig = x[i];
            double d = Abs(orig)*eps;
            if(d == 0) d = eps;
            x[i] = orig + d;
            double fp = phi(x);
            x[i] = orig - d;
            double fm = phi(x);
            grad[i] = (fp - fm)/(2*d);
            x[i] = orig;
        }
        return grad;
    }

    /// <summary>Central-difference Hessian.</summary>
    public static Mat CentralHessian(Func<Vec, double> phi, Vec x){
        int n = x.Length;
        Mat H = new Mat(n, n);
        double f0 = phi(x);
        double eps = Pow(2, -13);

        // Diagonal terms
        for(int i=0;i<n;i++){
            double orig = x[i];
            double d = Abs(orig)*eps;
            if(d == 0) d = eps;
            x[i] = orig + d; double fp = phi(x);
            x[i] = orig - d; double fm = phi(x);
            H[i,i] = (fp - 2*f0 + fm)/(d*d);
            x[i] = orig;
        }
        // Off-diagonal terms
        for(int i = 0; i < n; i++){
            for(int j = i+1; j < n; j++){
                double oi = x[i], oj = x[j];
                double di = Abs(oi)*eps; if(di == 0) di = eps;
                double dj = Abs(oj)*eps; if(dj == 0) dj = eps;

                x[i] = oi + di;  x[j] = oj + dj;  double fpp = phi(x);
                x[i] = oi + di;  x[j] = oj - dj;  double fpm = phi(x);
                x[i] = oi - di;  x[j] = oj + dj;  double fmp = phi(x);
                x[i] = oi - di;  x[j] = oj - dj;  double fmm = phi(x);

                double hij = (fpp - fpm - fmp + fmm) / (4 * di * dj);
                H[i,j] = hij;
                H[j,i] = hij;

                x[i] = oi;
                x[j] = oj;
            }
        }
        return H;
    }
}