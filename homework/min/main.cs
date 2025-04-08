using System;
using System.Collections.Generic;
using System.IO;

public class Program
{
    // Rosenbrock: f(x,y) = (1-x)² + 100*(y-x²)²
    public static double Rosenbrock(Vec v)
    {
        double x = v[0];
        double y = v[1];
        return Math.Pow(1 - x, 2) + 100 * Math.Pow(y - x * x, 2);
    }

    // Himmelblau: f(x,y) = (x²+y-11)² + (x+y²-7)²
    public static double Himmelblau(Vec v)
    {
        double x = v[0];
        double y = v[1];
        return Math.Pow(x * x + y - 11, 2) + Math.Pow(x + y * y - 7, 2);
    }

    // Breit–Wigner function: F(E|m,Γ,A) = A/[(E-m)² + Γ²/4]
    public static double BreitWigner(double E, double m, double Gamma, double A)
    {
        return A / (Math.Pow(E - m, 2) + Math.Pow(Gamma, 2) / 4.0);
    }

    // Deviation function: D(m,Γ,A)= Σi[((F(Ei|m,Γ,A)-σi)/Δσi)²].
    public static double Deviation(Vec p, double[] energy, double[] signal, double[] error)
    {
        double m = p[0];
        double Gamma = p[1];
        double A = p[2];
        double sum = 0.0;
        for (int i = 0; i < energy.Length; i++)
        {
            double diff = BreitWigner(energy[i], m, Gamma, A) - signal[i];
            sum += Math.Pow(diff / error[i], 2);
        }
        return sum;
    }

    public static void Main(string[] args)
    {
        // First, minimize the Rosenbrock function using forward differences.
        Vec x0 = new Vec(new double[] { -1, 1 });
        Vec solRosenbrock = NewtonMinimizer.Minimize(Rosenbrock, x0, 1e-3, 100, false);
        Console.WriteLine("Rosenbrock minimum:");
        Console.WriteLine("  x = " + solRosenbrock[0]);
        Console.WriteLine("  y = " + solRosenbrock[1]);
        Console.WriteLine("  Iterations: " + NewtonMinimizer.IterationCount);
        Console.WriteLine();

        // Next, minimize the Rosenbrock function with central differences.
        x0 = new Vec(new double[] { -1, 1 });
        solRosenbrock = NewtonMinimizer.Minimize(Rosenbrock, x0, 1e-3, 100, true);
        Console.WriteLine("Rosenbrock minimum (central differences):");
        Console.WriteLine("  x = " + solRosenbrock[0]);
        Console.WriteLine("  y = " + solRosenbrock[1]);
        Console.WriteLine("  Iterations: " + NewtonMinimizer.IterationCount);
        Console.WriteLine();

        // Next, minimize the Himmelblau function using forward differences.
        x0 = new Vec(new double[] { 2, 2 });
        Vec solHimmelblau = NewtonMinimizer.Minimize(Himmelblau, x0, 1e-3, 100, false);
        Console.WriteLine("Himmelblau minimum:");
        Console.WriteLine("  x = " + solHimmelblau[0]);
        Console.WriteLine("  y = " + solHimmelblau[1]);
        Console.WriteLine("  Iterations: " + NewtonMinimizer.IterationCount);
        Console.WriteLine();

        // Now, minimize the Himmelblau function with central differences.
        x0 = new Vec(new double[] { 2, 2 });
        solHimmelblau = NewtonMinimizer.Minimize(Himmelblau, x0, 1e-3, 100, true);
        Console.WriteLine("Himmelblau minimum (central differences):");
        Console.WriteLine("  x = " + solHimmelblau[0]);
        Console.WriteLine("  y = " + solHimmelblau[1]);
        Console.WriteLine("  Iterations: " + NewtonMinimizer.IterationCount);
        Console.WriteLine();

        // Load experimental data from a file.
        string dataFile = "higgs.data.txt";
        if (!File.Exists(dataFile))
        {
            Console.WriteLine("Data file \"" + dataFile + "\" not found.");
            return;
        }

        var energyList = new List<double>();
        var signalList = new List<double>();
        var errorList = new List<double>();
        string[] lines = File.ReadAllLines(dataFile);
        char[] separators = new char[] { ' ', '\t' };
        foreach (string line in lines)
        {
            if (string.IsNullOrWhiteSpace(line) || line.TrimStart().StartsWith("#"))
                continue;
            string[] words = line.Split(separators, StringSplitOptions.RemoveEmptyEntries);
            if (words.Length < 3)
                continue;
            energyList.Add(double.Parse(words[0]));
            signalList.Add(double.Parse(words[1]));
            errorList.Add(double.Parse(words[2]));
        }
        double[] energy = energyList.ToArray();
        double[] signal = signalList.ToArray();
        double[] error = errorList.ToArray();

        // Define the deviation function for Breit–Wigner fitting.
        // p[0] = mass m, p[1] = width Γ, p[2] = scale factor A.
        Func<Vec, double> D = p => Deviation(p, energy, signal, error);

        // Set an initial guess: mass ~125 GeV, width ~6 GeV, scale factor ~43.
        Vec initialGuess = new Vec(new double[] { 125.0, 6.0, 43.0 });

        // Perform the minimization/fitting (using default forward differences).
        Vec bestFit = NewtonMinimizer.Minimize(D, initialGuess, 1e-3, 50, false);
        double m_fit = bestFit[0];
        double Gamma_fit = bestFit[1];
        double A_fit = bestFit[2];
        int iterations = NewtonMinimizer.IterationCount;

        // Print the fit results.
        Console.WriteLine("Best-fit parameters for Breit-Wigner:");
        Console.WriteLine("  Mass m         = " + m_fit + " GeV/c²");
        Console.WriteLine("  Width Γ        = " + Gamma_fit + " GeV/c² (experimental upper limit)");
        Console.WriteLine("  Scale factor A = " + A_fit);
        Console.WriteLine("  Iterations     = " + iterations);

        // Generate a file for plotting: "plot.dat"
        // First section: experimental data; second section: the fitted curve.
        string outFile = "plot.dat";
        using (StreamWriter writer = new StreamWriter(outFile))
        {
            writer.WriteLine("# Experimental data: Energy [GeV]   Signal   Uncertainty");
            for (int i = 0; i < energy.Length; i++)
            {
                writer.WriteLine("{0} {1} {2}", energy[i], signal[i], error[i]);
            }
            writer.WriteLine(); // blank line to separate datasets
            writer.WriteLine("# Fitted Breit-Wigner curve: Energy [GeV]   Fitted Signal");
            double Emin = energy[0];
            double Emax = energy[energy.Length - 1];
            int Npoints = 200;
            double dE = (Emax - Emin) / (Npoints - 1);
            for (int i = 0; i < Npoints; i++)
            {
                double E = Emin + i * dE;
                double fitValue = BreitWigner(E, m_fit, Gamma_fit, A_fit);
                writer.WriteLine("{0} {1}", E, fitValue);
            }
        }
        Console.WriteLine("Plot data written to " + outFile);
    }
}
    