using NewtonSolver;
using System;
using System.IO;
using System.Collections.Generic;
using static ODESolver;
using static System.Math;
using static System.Console;

public class Program
{
    public static Vec Rosenbrock(Vec v)
    {
        // Extract x and y.
        double x = v[0];
        double y = v[1];
        // Return f(x) as a vector.
        return new Vec(new double[]
        {
            (1-x)*(1-x),
            100*(y-x*x)*(y-x*x)
        });
    }

public static Vec Himmelblau(Vec v)
{
    // Extract x and y.
    double x = v[0];
    double y = v[1];
    // Return f(x) as a vector.
    return new Vec(new double[]
    {
        (x * x + y - 11) * (x * x + y - 11),
        (x + y * y - 7) * (x + y * y - 7)
    });
}

// Nested static class containing the hydrogen shooting methods.
public static class HydrogenShootingProgram
{
    // ------------------------------------------------------------------
    // The ODE: f''(r) = -2*(E + 1/r)*f(r)
    // Written as a first-order system:
    //   y[0] = f,  y[1] = f'
    //   y[0]' = y[1]
    //   y[1]' = -2*(E + 1/r)*y[0]
    // ------------------------------------------------------------------
    public static Vec SchrodingerODE(double r, Vec y, double E)
    {
        Vec dy = new Vec(2);
        dy[0] = y[1];
        dy[1] = -2 * (E + 1 / r) * y[0];
        return dy;
    }

    // For a given energy E, integrate from rmin to rmax using our ODE solver.
    // The initial conditions are taken from the known small-r behavior:
    //   f(rmin) = rmin - rmin^2   and   f'(rmin) = 1 - 2*rmin.
    // Returns M(E) = f(rmax).
    public static double M(double E, double rmin, double rmax, double acc, double eps, double hstart)
    {
        Vec y0 = new Vec(new double[] { rmin - rmin * rmin, 1 - 2 * rmin });
        Func<double, Vec, Vec> F = (r, y) => SchrodingerODE(r, y, E);
        var result = ODESolver.driver(F, (rmin, rmax), y0, hstart, acc, eps);
        Vec yEnd = result.YList[result.YList.Count - 1];
        return yEnd[0];
    }

    // Wraps the scalar function M(E) as a function from 1D Vec to 1D Vec.
    // Given X = [E], returns [ M(E) ].
    public static Vec MEFunction(Vec X, double rmin, double rmax, double acc, double eps, double hstart)
    {
        double E = X[0];
        double mValue = M(E, rmin, rmax, acc, eps, hstart);
        return new Vec(new double[] { mValue });
    }
}

    public static void Main(string[] args)
    {
        // Parameters for the ODE integration and shooting method:
        double rmin = 1e-3;      // starting radius (small but > 0)
        double rmax = 8;         // ending radius (should be large enough)
        double acc = 1e-6;       // absolute accuracy for the ODE integrator
        double eps = 1e-6;       // relative accuracy for the ODE integrator
        double hstart = 0.001;   // initial step size for ODE integration

        // Starting guess for the optimization problems.
        Vec start = new Vec(new double[] { 0.1, 0.1 });
        // Solve f(x)=0 using Newton's method.
        Vec solution = Solver.Newton(Rosenbrock, start, 1e-6);
        Vec solution1 = Solver.Newton(Himmelblau, start, 1e-6);
        // Print the solutions.
        Console.WriteLine("Solution to Rosenbrock:");
        for (int i = 0; i < solution.Length; i++)
        {
            Console.WriteLine($"x[{i}] = {solution[i]}");
        }
        Console.WriteLine("Solution to Himmelblau:");
        for (int i = 0; i < solution1.Length; i++)
        {
            Console.WriteLine($"x[{i}] = {solution1[i]}");
        }

        // Use the Newton solver from roots.cs to find the energy E0 such that M(E0)=0.
        // Wrap the function into one that maps a 1D Vec -> 1D Vec.
        // For the hydrogen ground state the exact energy is -0.5.
        Vec Einitial = new Vec(new double[] { -0.5 });
        Func<Vec, Vec> F_E = X => HydrogenShootingProgram.MEFunction(X, rmin, rmax, acc, eps, hstart);
        Vec Eroot = Solver.Newton(F_E, Einitial, 1e-6);
        double E0 = Eroot[0];

        // Write the computed energy to a text file.
        using (StreamWriter energyWriter = new StreamWriter("EnergyOutput.txt"))
        {
            energyWriter.WriteLine("Computed bound state energy E0 = " + E0);
        }

        // With the computed E0, integrate the ODE once more to obtain the wavefunction.
        Func<double, Vec, Vec> F0 = (r, y) => HydrogenShootingProgram.SchrodingerODE(r, y, E0);
        Vec y0 = new Vec(new double[] { rmin - rmin * rmin, 1 - 2 * rmin });
        var sol = ODESolver.driver(F0, (rmin, rmax), y0, hstart, acc, eps);
        List<double> rValues = sol.XList;
        List<Vec> yValues = sol.YList;  // y[0] is f(r)

        // Write the wavefunction to a text file.
        // The exact ground state (up to normalization) is f_exact(r) = r * exp(-r)
        using (StreamWriter waveWriter = new StreamWriter("WavefunctionOutput.txt"))
        {
            waveWriter.WriteLine("r\tf_numeric\tf_exact");
            for (int i = 0; i < rValues.Count; i++)
            {
                double r = rValues[i];
                double f_numeric = yValues[i][0];
                double f_exact = r * Math.Exp(-r);
                waveWriter.WriteLine($"{r:G5}\t{f_numeric:G5}\t{f_exact:G5}");
            }
        }

        // Perform a convergence study with respect to rmax and rmin.
        using (StreamWriter convWriter = new StreamWriter("ConvergenceStudy.txt"))
        {
            convWriter.WriteLine("Convergence study with respect to rmax:");
            double[] rmax_values = new double[] {4, 6, 8, 10, 12, 14, 16 };
            foreach (double rm in rmax_values)
            {
                Vec Eguess = new Vec(new double[] { -0.5 });
                Func<Vec, Vec> F_E_rm = X => HydrogenShootingProgram.MEFunction(X, rmin, rm, acc, eps, hstart);
                Vec E_root_rm = Solver.Newton(F_E_rm, Eguess, 1e-6);
                convWriter.WriteLine($"rmax = {rm}, E0 = {E_root_rm[0]}");
            }

            convWriter.WriteLine("\nConvergence study with respect to rmin:");
            double[] rmin_values = new double[] {1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 0.00005, 0.00001 };
            foreach (double rmin_val in rmin_values)
            {
                Vec Eguess = new Vec(new double[] { -0.5 });
                Func<Vec, Vec> F_E_rmin = X => HydrogenShootingProgram.MEFunction(X, rmin_val, rmax, acc, eps, hstart);
                Vec E_root_rmin = Solver.Newton(F_E_rmin, Eguess, 1e-6);
                convWriter.WriteLine("rmin = {0:0.00000}, E0 = {1:0.#################}", rmin_val, E_root_rmin[0]);
            }

            convWriter.WriteLine("\nConvergence study with respect to acc:");
            double[] acc_values = new double[] { 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7 };
            foreach (double acc_val in acc_values)
            {
                Vec Eguess = new Vec(new double[] { -0.5 });
                Func<Vec, Vec> F_E_acc = X => HydrogenShootingProgram.MEFunction(X, rmin, rmax, acc_val, eps, hstart);
                Vec E_root_acc = Solver.Newton(F_E_acc, Eguess, 1e-6);
                convWriter.WriteLine("acc = {0:0.0000000}, E0 = {1:0.#################}", acc_val, E_root_acc[0]);
            }
            convWriter.WriteLine("\nConvergence study with respect to eps:");
            double[] eps_values = new double[] { 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7 };
            foreach (double eps_val in eps_values)
            {
                Vec Eguess = new Vec(new double[] { -0.5 });
                Func<Vec, Vec> F_E_eps = X => HydrogenShootingProgram.MEFunction(X, rmin, rmax, acc, eps_val, hstart);
                Vec E_root_eps = Solver.Newton(F_E_eps, Eguess, 1e-6);
                convWriter.WriteLine("eps = {0:0.0000000}, E0 = {1:0.#################}", eps_val, E_root_eps[0]);
            }
        }
    }
}
