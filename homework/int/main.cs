using System;
using System.IO;
using static System.Math;

class Program {

    /// <summary>
    /// Computes the error function erf(z) along with an estimated error.
    /// 
    /// Piecewise definition:
    ///   - For z < 0:    erf(z) = -erf(-z)
    ///   - For 0 ≤ z ≤ 1: erf(z) = (2/√π) ∫₀ᶻ e^(-x²) dx
    ///   - For z > 1:   erf(z) = 1 - (2/√π) ∫₀¹ [e^(-(z+(1-t)/t)²)/(t²)] dt
    /// 
    /// Note: IntegrationResult is now defined only in your integrator.dll.
    /// </summary>
    static IntegrationResult ErfWithError(double z)
    {
        if (z < 0) {
            IntegrationResult res = ErfWithError(-z);
            return new IntegrationResult { Value = -res.Value, Error = res.Error };
        }
        else if (z <= 1) {
            Func<double, double> integrand = x => Exp(-x * x);
            IntegrationResult res = Integrator.Integrate(integrand, 0, z);
            double factor = 2 / Sqrt(PI);
            return new IntegrationResult { Value = factor * res.Value, Error = factor * res.Error };
        }
        else {
            Func<double, double> integrand = t => Exp(-Pow(z + ((1-t)/t), 2)) / (t * t);
            IntegrationResult res = Integrator.Integrate(integrand, 0, 1);
            double factor = 2 / Sqrt(PI);
            return new IntegrationResult { Value = 1 - factor * res.Value, Error = factor * res.Error };
        }
    }

    static void Main()
    {

        // Reset counter before integration
        Integrator.ResetEvaluationCount();

        // Define the function to integrate: f(x) = sqrt(x)
        Func<double, double> f = x => Sqrt(x);
        var result = Integrator.Integrate(f, 0, 1);
        Console.WriteLine("Integral of sqrt(x) [0:1] is approximately: " + result.Value +
                          " with estimated error: " + result.Error);
        Console.WriteLine("Function evaluations: " + Integrator.EvaluationCount);

        // Reset counter before integration
        Integrator.ResetEvaluationCount();

        // Define the function to integrate: f(x) = 1/sqrt(x)
        Func<double, double> f1 = x1 => 1 / Sqrt(x1);
        var result1 = Integrator.Integrate(f1, 0, 1);
        Console.WriteLine("Integral of 1/sqrt(x) [0:1] is approximately: " + result1.Value +
                          " with estimated error: " + result1.Error);
        Console.WriteLine("Function evaluations: " + Integrator.EvaluationCount);

        // Reset counter before integration
        Integrator.ResetEvaluationCount();

        // Define the function to integrate: f(x) = 4*sqrt(1-x^2)
        Func<double, double> f2 = x2 => 4 * Sqrt(1 - x2 * x2);
        var result2 = Integrator.Integrate(f2, 0, 1);
        Console.WriteLine("Integral of 4*sqrt(1-x^2) [0:1] is approximately: " + result2.Value +
                          " with estimated error: " + result2.Error);
        Console.WriteLine("Function evaluations: " + Integrator.EvaluationCount);

        // Reset counter before integration
        Integrator.ResetEvaluationCount();

        // Define the function to integrate: f(x) = ln(x)/sqrt(x)
        Func<double, double> f3 = x3 => Log(x3) / Sqrt(x3);
        var result3 = Integrator.Integrate(f3, 0, 1);
        Console.WriteLine("Integral of ln(x)/sqrt(x) [0:1] is approximately: " + result3.Value +
                          " with estimated error: " + result3.Error);
        Console.WriteLine("Function evaluations: " + Integrator.EvaluationCount);

        // Compute and write the error function values for z from -3 to 3 into "erf.txt"
        using (StreamWriter writer = new StreamWriter("erf.txt"))
        {
            double step = 0.1;
            for (double z = -3.0; z <= 3.0; z += step)
            {
                IntegrationResult erfRes = ErfWithError(z);
                // Write z, erf(z), and the estimated error (tab-separated)
                writer.WriteLine($"{z}\t{erfRes.Value}\t{erfRes.Error}");
            }
        }
        Console.WriteLine("The error function values (with estimated errors) have been written to 'erf.txt'.");

        // Reset counter before integration
        Integrator.ResetEvaluationCount();

         // Define the function to integrate: f(x) = 1/sqrt(x)
        Func<double, double> f4 = x4 => 1 / Sqrt(x4);

        // Use the Clenshaw-Curtis transformation to compute the integral of f4 from 0 to 1
        var result4 = Integrator.CCIntegrate(f4, 0, 1);

        // Output the result and estimated error 
        Console.WriteLine("Integral of 1/sqrt(x) [0:1] using Clenshaw-Curtis transformation is approximately: " + result4.Value +
                          " with estimated error: " + result4.Error);
        Console.WriteLine("Function evaluations: " + Integrator.EvaluationCount);

        // Reset counter before integration
        Integrator.ResetEvaluationCount();

        // Define the function to integrate: f(x) = ln(x) / sqrt(x)
        Func<double, double> f5 = x5 => Log(x5) / Sqrt(x5);

        // Use the Clenshaw-Curtis transformation to compute the integral of f5 from 0 to 1
        var result5 = Integrator.CCIntegrate(f5, 0, 1);

        // Output the result and estimated error
        Console.WriteLine("Integral of ln(x) / sqrt(x) [0:1] using Clenshaw-Curtis transformation is approximately: " + result5.Value +
                          " with estimated error: " + result5.Error);
        Console.WriteLine("Function evaluations: " + Integrator.EvaluationCount);

        // Output that scipy uses 231 function evaluations for the same integral
        Console.WriteLine("Scipy uses 231 function evaluations for the same integrals.");
    }  
}