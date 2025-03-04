// Spline.cs
using System;

public enum SplineType
{
    Linear,
    Quadratic,
    Cubic
}

public class Spline
{
    // Node arrays (copied from input)
    private double[] x;   // Knot positions
    private double[] y;   // Function values at knots
    private SplineType type;

    // Coefficient arrays for each interval (i=0..n-2)
    // Linear: b[i] = slope on interval i.
    // Quadratic: b[i] = derivative estimate at x[i], c[i] = quadratic coefficient.
    // Cubic: b[i], c[i], d[i] such that:
    //   S_i(x) = y[i] + b[i]*(x-x[i]) + c[i]*(x-x[i])^2 + d[i]*(x-x[i])^3.
    private double[] b;
    private double[] c;
    private double[] d;  // Only used for cubic

    // Constructor: copies input arrays and builds the spline of the specified type.
    public Spline(double[] xs, double[] ys, SplineType splineType)
    {
        if (xs.Length != ys.Length)
            throw new ArgumentException("x and y must have the same length.");
        if (xs.Length < 2)
            throw new ArgumentException("At least two data points are required.");

        x = (double[])xs.Clone();
        y = (double[])ys.Clone();
        type = splineType;

        switch (type)
        {
            case SplineType.Linear:
                BuildLinearCoefficients();
                break;
            case SplineType.Quadratic:
                BuildQuadraticCoefficients();
                break;
            case SplineType.Cubic:
                BuildCubicCoefficients();
                break;
            default:
                throw new ArgumentException("Unknown spline type.");
        }
    }

    // =====================================================
    // Build routines for each spline type
    // =====================================================

    // Linear spline: precompute the slope on each interval.
    private void BuildLinearCoefficients()
    {
        int n = x.Length;
        b = new double[n - 1];
        for (int i = 0; i < n - 1; i++)
        {
            double dx = x[i + 1] - x[i];
            b[i] = (y[i + 1] - y[i]) / dx;
        }
    }

    // Quadratic spline:
    // For each interval [x[i], x[i+1]], define
    //   S_i(x) = y[i] + d_i*(x - x[i]) + a_i*(x - x[i])^2.
    // We estimate the derivative d_i at x[i]:
    //   - For i = 0, use forward difference.
    //   - For interior nodes, use central difference.
    // Then set:
    //   a_i = (y[i+1] - y[i] - d_i*(x[i+1]-x[i])) / (x[i+1]-x[i])^2.
    private void BuildQuadraticCoefficients()
    {
        int n = x.Length;
        int m = n - 1;  // number of intervals
        b = new double[m];
        c = new double[m];

        // Compute derivative estimates at the knots.
        double[] dEst = new double[n];
        dEst[0] = (y[1] - y[0]) / (x[1] - x[0]);
        for (int i = 1; i < n - 1; i++)
        {
            dEst[i] = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1]);
        }
        dEst[n - 1] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);

        // For each interval, store derivative and quadratic coefficient.
        for (int i = 0; i < m; i++)
        {
            double h = x[i + 1] - x[i];
            b[i] = dEst[i];
            c[i] = (y[i + 1] - y[i] - b[i] * h) / (h * h);
        }
    }

    // Cubic spline (Natural spline):
    // Solve for second derivatives m[i] with m[0] = m[n-1] = 0.
    // Then for each interval:
    //   b[i] = (y[i+1] - y[i])/h_i - h_i*(2*m[i] + m[i+1])/6,
    //   c[i] = m[i]/2,
    //   d[i] = (m[i+1] - m[i])/(6*h_i).
    private void BuildCubicCoefficients()
    {
        int n = x.Length;
        int m_intervals = n - 1;
        double[] h = new double[m_intervals];
        for (int i = 0; i < m_intervals; i++)
        {
            h[i] = x[i + 1] - x[i];
        }

        // Set up tridiagonal system for second derivatives m[i].
        double[] mu = new double[n];      // subdiagonal factors (mu[0] unused)
        double[] lambda = new double[n];  // superdiagonal factors (lambda[n-1] unused)
        double[] d_sys = new double[n];     // right-hand side, d_sys[0]=d_sys[n-1]=0
        mu[0] = 0;
        lambda[n - 1] = 0;
        d_sys[0] = 0;
        d_sys[n - 1] = 0;
        for (int i = 1; i < n - 1; i++)
        {
            double hi = h[i - 1];
            double hi1 = h[i];
            double denom = hi + hi1;
            mu[i] = hi / denom;
            lambda[i] = hi1 / denom;
            d_sys[i] = 6 * ((y[i + 1] - y[i]) / hi1 - (y[i] - y[i - 1]) / hi) / denom;
        }

        // Solve the tridiagonal system using the Thomas algorithm.
        double[] m_arr = new double[n];
        m_arr[0] = 0;
        m_arr[n - 1] = 0;
        double[] cprime = new double[n];
        double[] dprime = new double[n];
        cprime[1] = lambda[1] / 2.0;  // denominator for i=1 is 2
        dprime[1] = d_sys[1] / 2.0;
        for (int i = 2; i < n - 1; i++)
        {
            double denom = 2 - mu[i] * cprime[i - 1];
            cprime[i] = lambda[i] / denom;
            dprime[i] = (d_sys[i] - mu[i] * dprime[i - 1]) / denom;
        }
        for (int i = n - 2; i >= 1; i--)
        {
            m_arr[i] = dprime[i] - cprime[i] * m_arr[i + 1];
        }

        // Allocate coefficient arrays for each interval.
        b = new double[m_intervals];
        c = new double[m_intervals];
        d = new double[m_intervals];
        for (int i = 0; i < m_intervals; i++)
        {
            double hi = h[i];
            b[i] = (y[i + 1] - y[i]) / hi - hi * (2 * m_arr[i] + m_arr[i + 1]) / 6.0;
            c[i] = m_arr[i] / 2.0;
            d[i] = (m_arr[i + 1] - m_arr[i]) / (6.0 * hi);
        }
    }

    // =====================================================
    // Public methods: Evaluate, Derivative, Integral
    // =====================================================

    public double Evaluate(double z)
    {
        switch (type)
        {
            case SplineType.Linear:
                return EvaluateLinear(z);
            case SplineType.Quadratic:
                return EvaluateQuadratic(z);
            case SplineType.Cubic:
                return EvaluateCubic(z);
            default:
                throw new InvalidOperationException("Unknown spline type.");
        }
    }

    public double Derivative(double z)
    {
        switch (type)
        {
            case SplineType.Linear:
                return DerivativeLinear(z);
            case SplineType.Quadratic:
                return DerivativeQuadratic(z);
            case SplineType.Cubic:
                return DerivativeCubic(z);
            default:
                throw new InvalidOperationException("Unknown spline type.");
        }
    }

    public double Integral(double z)
    {
        switch (type)
        {
            case SplineType.Linear:
                return IntegralLinear(z);
            case SplineType.Quadratic:
                return IntegralQuadratic(z);
            case SplineType.Cubic:
                return IntegralCubic(z);
            default:
                throw new InvalidOperationException("Unknown spline type.");
        }
    }

    // =====================================================
    // Linear methods
    // =====================================================
    private double EvaluateLinear(double z)
    {
        int i = FindInterval(z);
        return y[i] + b[i] * (z - x[i]);
    }

    private double DerivativeLinear(double z)
    {
        int i = FindInterval(z);
        return b[i];
    }

    private double IntegralLinear(double z)
    {
        int i = FindInterval(z);
        double area = 0.0;
        for (int k = 0; k < i; k++)
        {
            double dx = x[k + 1] - x[k];
            area += 0.5 * (y[k] + y[k + 1]) * dx;
        }
        double dx_last = z - x[i];
        double y_z = y[i] + b[i] * dx_last;
        area += 0.5 * (y[i] + y_z) * dx_last;
        return area;
    }

    // =====================================================
    // Quadratic methods
    // =====================================================
    private double EvaluateQuadratic(double z)
    {
        int i = FindInterval(z);
        double t = z - x[i];
        return y[i] + b[i] * t + c[i] * t * t;
    }

    private double DerivativeQuadratic(double z)
    {
        int i = FindInterval(z);
        double t = z - x[i];
        return b[i] + 2 * c[i] * t;
    }

    private double IntegralQuadratic(double z)
    {
        int i = FindInterval(z);
        double area = 0.0;
        // Sum complete intervals
        for (int k = 0; k < i; k++)
        {
            double h = x[k + 1] - x[k];
            area += y[k] * h + 0.5 * b[k] * h * h + (1.0 / 3.0) * c[k] * h * h * h;
        }
        // Partial interval [x[i], z]
        double t = z - x[i];
        area += y[i] * t + 0.5 * b[i] * t * t + (1.0 / 3.0) * c[i] * t * t * t;
        return area;
    }

    // =====================================================
    // Cubic methods
    // =====================================================
    private double EvaluateCubic(double z)
    {
        int i = FindInterval(z);
        double t = z - x[i];
        return y[i] + b[i] * t + c[i] * t * t + d[i] * t * t * t;
    }

    private double DerivativeCubic(double z)
    {
        int i = FindInterval(z);
        double t = z - x[i];
        return b[i] + 2 * c[i] * t + 3 * d[i] * t * t;
    }

    private double IntegralCubic(double z)
    {
        int i = FindInterval(z);
        double area = 0.0;
        // Sum complete intervals
        for (int k = 0; k < i; k++)
        {
            double h = x[k + 1] - x[k];
            area += y[k] * h + 0.5 * b[k] * h * h + (1.0 / 3.0) * c[k] * h * h * h + 0.25 * d[k] * h * h * h * h;
        }
        // Partial interval [x[i], z]
        double t = z - x[i];
        area += y[i] * t + 0.5 * b[i] * t * t + (1.0 / 3.0) * c[i] * t * t * t + 0.25 * d[i] * t * t * t * t;
        return area;
    }

    // =====================================================
    // Helper: binary search to find the interval index i such that x[i] <= z <= x[i+1]
    // =====================================================
    private int FindInterval(double z)
    {
        int n = x.Length;
        if (z <= x[0])
            return 0;
        if (z >= x[n - 1])
            return n - 2;
        int i = 0, j = n - 1;
        while (j - i > 1)
        {
            int mid = (i + j) / 2;
            if (z >= x[mid])
                i = mid;
            else
                j = mid;
        }
        return i;
    }
}