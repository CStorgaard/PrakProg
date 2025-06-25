using System;

namespace BerrutInterpolation
{
    /// <summary>
    /// Berrut’s rational function interpolator.
    /// </summary>
    public class BerrutInterpolator
    {
        private readonly double[] x; // interpolation nodes
        private readonly double[] y; // function values at the nodes
        private readonly double[] w; // weights
        private readonly int n; // number of nodes - 1

        /// <summary>
        /// Constructs a Berrut interpolator.
        /// </summary>
        /// <param name="nodes">Distinct interpolation nodes, length n+1.</param>
        /// <param name="values">Function values at the nodes; same length as nodes.</param>
        /// <param name="linearReproduction">
        /// If true, uses the “second” form (reproduces linear polynomials).
        /// If false, uses the “first” form (reproduces constants).
        /// </param>
        public BerrutInterpolator(double[] nodes, double[] values, bool linearReproduction = false)
        {
            if (nodes == null || values == null)
                throw new ArgumentNullException("Input arrays cannot be null.");
            if (nodes.Length != values.Length || nodes.Length < 2)
                throw new ArgumentException("nodes and values must have the same length >= 2.");
            n = nodes.Length - 1;

            x = new double[n + 1];
            y = new double[n + 1];
            Array.Copy(nodes, x, n + 1);
            Array.Copy(values, y, n + 1);

            w = new double[n + 1];
            ComputeWeights(linearReproduction);
        }

        /// <summary>
        /// Precompute weights.
        /// </summary>
        private void ComputeWeights(bool linearReproduction)
        {
            // First form: w_j = (-1)^j
            // Second form: w_0 = 1/2, w_n = (-1)^n/2, w_j = (-1)^j for 1 <= j <= n-1
            for (int j = 0; j <= n; j++)
            {
                double sign = (j % 2 == 0) ? +1.0 : -1.0;
                if (linearReproduction)
                {
                    if (j == 0)
                        w[j] = 0.5;
                    else if (j == n)
                        w[j] = sign * 0.5;
                    else
                        w[j] = sign;
                }
                else
                {
                    w[j] = sign;
                }
            }
        }

        /// <summary>
        /// Evaluate the interpolant at z.
        /// </summary>
        public double Evaluate(double z)
        {
            // If z exactly matches a node, return the known value.
            for (int j = 0; j <= n; j++)
            {
                if (z == x[j])
                    return y[j];
            }

            // Otherwise form the barycentric sums.
            double num = 0.0;
            double den = 0.0;
            for (int j = 0; j <= n; j++)
            {
                double diff = z - x[j];
                double temp = w[j] / diff;
                num += temp * y[j];
                den += temp;
            }

            return num / den;
        }
    }
}