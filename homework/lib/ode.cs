using System;
using System.Collections.Generic;
using static System.Math;

/// <summary>
/// Holds the result of one Runge–Kutta step.
/// </summary>
public class RKResult
{
    public Vec YNext { get; }
    public Vec ErrorEstimate { get; }

    public RKResult(Vec yNext, Vec errorEstimate)
    {
        YNext = yNext;
        ErrorEstimate = errorEstimate;
    }
}

/// <summary>
/// Holds the result of an ODE integration: lists of x values and corresponding Vec solutions.
/// </summary>
public class ODEResult
{
    public List<double> XList { get; }
    public List<Vec> YList { get; }

    public ODEResult(List<double> xList, List<Vec> yList)
    {
        XList = xList;
        YList = yList;
    }
}

/// <summary>
/// Provides ODE solving methods using an adaptive Runge–Kutta method.
/// </summary>
public static class ODESolver
{
    /// <summary>
    /// Computes one integration step of an embedded Runge–Kutta method (orders 1 and 2)
    /// to solve the ODE dy/dx = f(x, y) using both Euler (first order) and midpoint (second order) methods.
    /// </summary>
    public static RKResult RKStep12(
        Func<double, Vec, Vec> f, // The derivative function f(x, y)
        double x,                 // The current x value
        Vec y,                    // The current state vector y(x)
        double h                  // The step size to be taken
    )
    {
        if (f == null)
            throw new ArgumentNullException("f", "The derivative function cannot be null.");

        Vec k0 = f(x, y);                     // Compute slope using Euler's method (first order estimate)
        Vec k1 = f(x + h / 2, y + k0 * (h / 2)); // Compute slope at the midpoint (second order estimate)

        Vec yNext = y + k1 * h;                 // Estimate y(x+h) using the midpoint method
        Vec errorEstimate = (k1 - k0) * h;      // Estimate the error as the difference between midpoint and Euler estimates

        return new RKResult(yNext, errorEstimate); // Return the result as an RKResult object
    }

    /// <summary>
    /// Integrates the ODE dy/dx = F(x, y) over a given interval using an adaptive step-size method
    /// based on the embedded Runge–Kutta (orders 1 and 2) method. Returns lists of x values and y vectors.
    /// </summary>
    public static ODEResult driver(
        Func<double, Vec, Vec> F,            // The ODE function F(x, y)
        (double start, double end) interval,  // The integration interval as (start, end)
        Vec yStart,                          // The initial state vector at the start of the interval
        double h = 0.125,                    // Initial step size
        double acc = 0.01,                   // Absolute accuracy goal
        double eps = 0.01                    // Relative accuracy goal
    )
    {
        var (a, b) = interval;               // Unpack the integration interval into a and b
        double x = a;                      // Initialize x to the start value

        Vec y = yStart.Copy();               // Copy the initial vector using your Copy() method

        List<double> xList = new List<double>(); // List to store x values
        xList.Add(x);                        // Add the starting x value

        List<Vec> yList = new List<Vec>();   // List to store state vectors
        yList.Add(y);                        // Add the initial state vector

        while (true)
        {
            if (x >= b)                    // If we've reached or passed the end of the interval,
                return new ODEResult(xList, yList); // return the result

            if (x + h > b)                 // If the next step overshoots the end,
                h = b - x;                 // adjust h so the final step lands exactly at b

            RKResult result = RKStep12(F, x, y, h); // Compute one Runge–Kutta step and error estimate

            // Calculate the tolerance based on absolute and relative accuracy goals,
            // scaled by the square root of the step's fraction of the total interval.
            double tol = (acc + eps * result.YNext.Norm()) * Math.Sqrt(h / (b - a));
            double err = result.ErrorEstimate.Norm();  // Compute the norm of the error estimate

            if (err <= tol)                // If the error is within tolerance,
            {
                x += h;                  // update x by the step size,
                y = result.YNext;        // accept the new state vector,
                xList.Add(x);            // record x,
                yList.Add(y);            // record y.
            }
            // Readjust the step size:
            // Scale h by (tol/err)^0.25, apply a safety factor of 0.95, and cap the increase to a factor of 2.
            h *= Math.Min(Math.Pow(tol / err, 0.25) * 0.95, 2);
        }
    }
}

public static class Interpolator
{
    // Helper: Binary search to find the index i such that x[i] <= z < x[i+1].
    private static int BinSearch(List<double> x, double z)
    {
        int low = 0, high = x.Count - 1;
        while (high - low > 1)
        {
            int mid = (low + high) / 2;
            if (x[mid] > z)
                high = mid;
            else
                low = mid;
        }
        return low;
    }

    /// <summary>
    /// Returns a cubic Hermite interpolant based on the given data points.
    /// This function first computes approximated derivatives at the nodes using finite differences.
    /// Then, for any z in [x[0], x[n-1]], it performs cubic Hermite interpolation.
    /// </summary>
    public static Func<double, Vec> MakeCubicInterpolant(List<double> x, List<Vec> y)
    {
        int n = x.Count;
        if (n < 2)
            throw new ArgumentException("At least two data points are required.", nameof(x));

        // Compute derivative approximations for each node.
        // We assume that Vec supports subtraction, division by a scalar, and scalar multiplication.
        List<Vec> d = new List<Vec>(n);

        // Forward difference for the first point.
        double h0 = x[1] - x[0];
        d.Add((y[1] - y[0]) / h0);

        // Central differences for interior points.
        for (int i = 1; i < n - 1; i++)
        {
            double h = x[i + 1] - x[i - 1];
            d.Add((y[i + 1] - y[i - 1]) / h);
        }

        // Backward difference for the last point.
        double hLast = x[n - 1] - x[n - 2];
        d.Add((y[n - 1] - y[n - 2]) / hLast);

        // Return the interpolant function.
        return delegate (double z)
        {
            // Clamp z to the interpolation range.
            if (z <= x[0])
                return y[0];
            if (z >= x[n - 1])
                return y[n - 1];

            // Find the interval index i so that x[i] <= z < x[i+1].
            int i = BinSearch(x, z);
            double hInterval = x[i + 1] - x[i];
            double t = (z - x[i]) / hInterval;

            // Cubic Hermite basis functions.
            double H00 = 2 * t * t * t - 3 * t * t + 1;
            double H10 = t * t * t - 2 * t * t + t;
            double H01 = -2 * t * t * t + 3 * t * t;
            double H11 = t * t * t - t * t;

            // Return the interpolated vector:
            // y[i]*H00 + h*d[i]*H10 + y[i+1]*H01 + h*d[i+1]*H11.
            return y[i] * H00 + d[i] * (hInterval * H10) + y[i + 1] * H01 + d[i + 1] * (hInterval * H11);
        };
    }

    /// <summary>
    /// Solves the ODE IVP using the driver and returns an interpolant for the solution.
    /// </summary>
    /// <returns>A function that interpolates the solution, so that for any z in the interval, you get y(z).</returns>
    public static Func<double, Vec> MakeOdeIvpInterpolant(
        Func<double, Vec, Vec> f,
        (double start, double end) interval,
        Vec y,
        double acc = 0.01,
        double eps = 0.01,
        double hstart = 0.01)
    {
        // Assume ODESolver.driver returns a result containing List<double> XList and List<Vec> YList.
        ODEResult result = ODESolver.driver(f, interval, y, hstart, acc, eps);
        return MakeCubicInterpolant(result.XList, result.YList);
    }
}
