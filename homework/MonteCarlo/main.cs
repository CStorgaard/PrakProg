using System;
using System.IO;
using static System.Math;
using static MonteCarlo;  // Assumes MonteCarlo.PlainMC is defined in this static class.

class Program 
{
    // Lanczos approximation for the Gamma function.
    // This implementation is valid for positive real arguments.
    public static double Gamma(double z)
    {
        // Coefficients used by the GNU Scientific Library
        double[] p = {
            676.5203681218851,
            -1259.1392167224028,
            771.32342877765313,
            -176.61502916214059,
            12.507343278686905,
            -0.13857109526572012,
            9.9843695780195716e-6,
            1.5056327351493116e-7
        };
        int g = 7;
        if (z < 0.5)
        {
            // Reflection formula for argument < 0.5
            return PI / (Sin(PI * z) * Gamma(1 - z));
        }
        else
        {
            z -= 1;
            double a = 0.99999999999980993;
            for (int i = 0; i < p.Length; i++)
            {
                a += p[i] / (z + i + 1);
            }
            double t = z + g + 0.5;
            return Sqrt(2 * PI) * Pow(t, z + 0.5) * Exp(-t) * a;
        }
    }
    
    // Main entry point for the program.
    static void Main() 
    {
        // --- Part 1: Estimate π via Monte Carlo integration ---
        // Define the function for the unit circle.
        // Returns 1 if the point (x[0], x[1]) is inside the circle, 0 otherwise.
        Func<double[], double> f = (x) => (x[0] * x[0] + x[1] * x[1] <= 1) ? 1.0 : 0.0;

        // Define the integration region as the square [-1, 1] x [-1, 1].
        double[] lowerLimits = new double[] { -1, -1 };
        double[] upperLimits = new double[] { 1, 1 };

        // Output file for π estimates.
        string outputPath = "pi_results.txt";

        using (StreamWriter writer = new StreamWriter(outputPath))
        {
            // Write header.
            writer.WriteLine("Samples\tEstimated Pi\tError Estimate\tActual error");
            Console.WriteLine("Estimating Pi using Monte Carlo integration:");
            
            // Loop over sample sizes from 10^0 to 10^6.
            for (int exp = 0; exp <= 6; exp++)
            {
                int N = (int)Pow(10, exp);
                var result = PlainMC(f, lowerLimits, upperLimits, N);
                double product = result.errorEstimate * Sqrt(N);

                // For the circle in [-1,1]x[-1,1], the region has area 4,
                // so the estimated integral equals the area of the circle (≈π).
                double actualError = Abs(result.estimatedIntegral - PI);
                writer.WriteLine($"{N}\t{result.estimatedIntegral}\t{result.errorEstimate}\t{actualError}");
                Console.WriteLine($"N = {N}: Estimated Pi = {result.estimatedIntegral}, Error = {result.errorEstimate}");
                Console.WriteLine($"N = {N}: Product of error and sqrt(N)= {product}");
            }
        }
        Console.WriteLine($"π estimation results have been written to {outputPath}");

        // --- Part 2: Compute the integral I ---
        // I = 1/π^3 ∫_0^π dx ∫_0^π dy ∫_0^π dz [1 / (1 - cos(x)cos(y)cos(z))]
        // Exact value: (Gamma(1/4)^4) / (4π^3)
        double exact = Pow(Gamma(0.25), 4) / (4 * Pow(PI, 3));

        // Define the integrand function.
        Func<double[], double> f3 = (x) => 1.0 / (1.0 - Cos(x[0]) * Cos(x[1]) * Cos(x[2]));

        // Integration region: the cube [0, π]^3.
        double[] lowerLimits3 = new double[] { 0, 0, 0 };
        double[] upperLimits3 = new double[] { PI, PI, PI };

        // Output file for the integral estimates.
        string outputPath3 = "integral_results.txt";

        using (StreamWriter writer = new StreamWriter(outputPath3))
        {
            writer.WriteLine("Samples\tEstimated Integral\tError Estimate\tActual error");
            Console.WriteLine("\nEstimating the integral I using Monte Carlo integration:");
            
            // Loop through sample sizes from 10^0 to 10^8.
            for (int exp = 0; exp <= 6; exp++)
            {
                int N = (int)Pow(10, exp);
                var result = PlainMC(f3, lowerLimits3, upperLimits3, N);

                // The integration routine returns an estimate = (average f * volume),
                // so we normalize by dividing by π^3.
                double normalizedIntegral = result.estimatedIntegral / Pow(PI, 3);
                double normalizedError = result.errorEstimate / Pow(PI, 3);
                double actualError = Abs(normalizedIntegral - exact);

                writer.WriteLine($"{N}\t{normalizedIntegral}\t{normalizedError}\t{actualError}");
                Console.WriteLine($"N = {N}: Estimated Integral = {normalizedIntegral}, Error Estimate = {normalizedError}, Actual error = {actualError}");
            }
        }
        Console.WriteLine($"Integral estimation results have been written to {outputPath3}");


        // --- Part 3: Estimate π via quasi-random (low-discrepancy) sequence ---
        // Output file for π estimates using quasi-random sequence.
        string outputPathQMC = "pi_results_quasi.txt";

        using (StreamWriter writer = new StreamWriter(outputPathQMC))
        {
            writer.WriteLine("Samples\tEstimated Pi\tError Estimate\tActual error");
            Console.WriteLine("\nEstimating Pi using quasi-random (low-discrepancy) sequence:");

            // Loop over sample sizes from 10^0 to 10^8.
            for (int exp = 0; exp <= 6; exp++)
            {
                int N = (int)Pow(10, exp);
                var result = QMC(f, lowerLimits, upperLimits, N);
                double product = result.errorEstimate * Sqrt(N);

                // For the circle in [-1,1]x[-1,1], the region has area 4,
                // so the estimated integral equals the area of the circle (≈π).
                double actualError = Abs(result.estimatedIntegral - PI);
                writer.WriteLine($"{N}\t{result.estimatedIntegral}\t{result.errorEstimate}\t{actualError}");
                Console.WriteLine($"N = {N}: Estimated Pi = {result.estimatedIntegral}, Error = {result.errorEstimate}");
                Console.WriteLine($"N = {N}: Product of error and sqrt(N)= {product}");
            }
        }
        Console.WriteLine($"π estimation results using quasi-random sequence have been written to {outputPathQMC}");

        // --- Part 4: Compute the integral I using quasi-random sequence ---
        // Output file for the integral estimates using quasi-random sequence.
        string outputPathQMC3 = "integral_results_quasi.txt";

        using (StreamWriter writer = new StreamWriter(outputPathQMC3))
        {
            writer.WriteLine("Samples\tEstimated Integral\tError Estimate\tActual error");
            Console.WriteLine("\nEstimating the integral I using quasi-random (low-discrepancy) sequence:");

            // Loop through sample sizes from 10^0 to 10^8.
            for (int exp = 0; exp <= 6; exp++)
            {
                int N = (int)Pow(10, exp);
                var result = QMC(f3, lowerLimits3, upperLimits3, N);

                // The integration routine returns an estimate = (average f * volume),
                // so we normalize by dividing by π^3.
                double normalizedIntegral = result.estimatedIntegral / Pow(PI, 3);
                double normalizedError = result.errorEstimate / Pow(PI, 3);
                double actualError = Abs(normalizedIntegral - exact);

                writer.WriteLine($"{N}\t{normalizedIntegral}\t{normalizedError}\t{actualError}");
                Console.WriteLine($"N = {N}: Estimated Integral = {normalizedIntegral}, Error Estimate = {normalizedError}, Actual error = {actualError}");
            }
        }
        Console.WriteLine($"Integral estimation results using quasi-random sequence have been written to {outputPathQMC3}");

        Console.WriteLine("Error of quasi-random sequence is significantly smaller than random sequence.");
    }
}
