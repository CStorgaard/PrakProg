using System;

/// <summary>
/// Provides methods for Monte Carlo integration.
/// </summary>
public static class MonteCarlo
{
    /// <summary>
    /// Approximates the integral of a function over a multidimensional rectangular region using the Monte Carlo method.
    /// </summary>
    /// <param name="f">The function to integrate. It accepts a point (array of doubles) and returns a double.</param>
    /// <param name="a">An array representing the lower bounds of the integration region.</param>
    /// <param name="b">An array representing the upper bounds of the integration region.</param>
    /// <param name="N">The number of random samples to generate.</param>
    /// <returns>A tuple where the first element is the estimated integral and the second element is the error estimate.</returns>
    public static (double estimatedIntegral, double errorEstimate) PlainMC(
        Func<double[], double> f, double[] a, double[] b, int N)
    {
        // Ensure that the lower and upper bounds have the same number of dimensions.
        int dim = a.Length;
        if (b.Length != dim)
            throw new ArgumentException("Lower and upper bounds must have the same number of dimensions.");

        // Calculate the volume of the integration region.
        double volume = 1;
        for (int i = 0; i < dim; i++)
        {
            volume *= b[i] - a[i];
        }

        double sum = 0;
        double sumOfSquares = 0;
        double[] x = new double[dim];
        Random rnd = new Random();

        // Perform Monte Carlo sampling.
        for (int i = 0; i < N; i++)
        {
            // Generate a random point within the integration region.
            for (int k = 0; k < dim; k++)
            {
                x[k] = a[k] + rnd.NextDouble() * (b[k] - a[k]);
            }
            double fx = f(x);
            sum += fx;
            sumOfSquares += fx * fx;
        }

        // Compute the mean and standard deviation.
        double mean = sum / N;
        double sigma = Math.Sqrt(sumOfSquares / N - mean * mean);

        // Scale the mean and sigma by the volume of the region.
        double estimatedIntegral = mean * volume;
        double errorEstimate = sigma * volume / Math.Sqrt(N);

        return (estimatedIntegral, errorEstimate);
    }

    /// <summary>
    /// Approximates the integral of a function over a multidimensional rectangular region using a quasi-random (low-discrepancy) sequence.
    /// This implementation uses the Van der Corput sequence for each dimension (with different prime bases) to form a Halton sequence.
    /// The error is estimated by computing two independent quasi-random estimates (using different index offsets) and comparing them.
    /// </summary>
    /// <param name="f">The function to integrate. It accepts a point (array of doubles) and returns a double.</param>
    /// <param name="a">An array representing the lower bounds of the integration region.</param>
    /// <param name="b">An array representing the upper bounds of the integration region.</param>
    /// <param name="N">The number of quasi-random samples per sequence.</param>
    /// <returns>A tuple where the first element is the estimated integral and the second element is an error estimate.</returns>
    public static (double estimatedIntegral, double errorEstimate) QMC(
        Func<double[], double> f, double[] a, double[] b, int N)
    {
        int dim = a.Length;
        if (b.Length != dim)
            throw new ArgumentException("Lower and upper bounds must have the same number of dimensions.");

        // Calculate the volume of the integration region.
        double volume = 1;
        for (int i = 0; i < dim; i++)
        {
            volume *= b[i] - a[i];
        }

        // We'll compute two independent quasi-random estimates.
        double sum1 = 0;
        double sum2 = 0;
        double[] x = new double[dim];

        // For multi-dimensional quasi-MC we use a Halton sequence.
        // We use a different prime base for each dimension.
        int[] bases = new int[] { 2, 3, 5, 7, 11, 13, 17, 19 };

        // First quasi-random sequence: use indices (i + 1)
        for (int i = 0; i < N; i++)
        {
            for (int k = 0; k < dim; k++)
            {
                double vdc = VanDerCorput(i + 1, bases[k]);
                x[k] = a[k] + vdc * (b[k] - a[k]);
            }
            double fx = f(x);
            sum1 += fx;
        }
        double estimate1 = (sum1 / N) * volume;

        // Second quasi-random sequence: use a different offset (i + N + 1)
        for (int i = 0; i < N; i++)
        {
            for (int k = 0; k < dim; k++)
            {
                double vdc = VanDerCorput(i + N + 1, bases[k]);
                x[k] = a[k] + vdc * (b[k] - a[k]);
            }
            double fx = f(x);
            sum2 += fx;
        }
        double estimate2 = (sum2 / N) * volume;

        // Combine the two estimates.
        double estimate = (estimate1 + estimate2) / 2.0;
        // Use the absolute difference as an error estimate.
        double errorEstimate = Math.Abs(estimate1 - estimate2);
        return (estimate, errorEstimate);
    }

    /// <summary>
    /// Computes the Van der Corput sequence value for index n in the given base.
    /// </summary>
    /// <param name="n">The index (should be non-negative, typically n+1 is used for the first sample).</param>
    /// <param name="baseValue">The base for the sequence (typically a prime number).</param>
    /// <returns>A quasi-random number in the interval [0,1).</returns>
    private static double VanDerCorput(int n, int baseValue)
    {
        double vdc = 0;
        double denominator = baseValue;
        while (n > 0)
        {
            int remainder = n % baseValue;
            vdc += remainder / denominator;
            n /= baseValue;
            denominator *= baseValue;
        }
        return vdc;
    }
}
