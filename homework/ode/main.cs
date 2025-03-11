using System;
using System.IO;
using System.Collections.Generic;
using static System.Math;

class Program
    {
    static void Main(string[] args)
    {
        // Choose simulation based on command-line arguments.
        // "friction" runs the damped pendulum,
        // "planet" or "planet_orbit_gr" runs the three-scenario planet orbit simulation,
        // otherwise, the harmonic oscillator is simulated.
        if (args.Length > 0)
        {
            string sim = args[0].ToLower();
            if (sim == "friction")
                RunFrictionHarmonicOscillator();
            else if (sim == "planet" || sim == "planet_orbit_gr")
                RunPlanetOrbitGRMultiple();
            else
                RunHarmonicOscillator();
        }
        else
        {
            RunHarmonicOscillator();
        }
    }


    static void RunHarmonicOscillator()
    {
        // Define the ODE for the harmonic oscillator: y'' = -y.
        // Convert it to a first-order system: let u = y and v = y', so that u' = v and v' = -u.
        Func<double, Vec, Vec> harmonicOscillator = (x, y) =>
        {
            Vec dydx = new Vec(2);
            dydx[0] = y[1];    // u' = v
            dydx[1] = -y[0];   // v' = -u
            return dydx;
        };

        // Initial conditions: y(0)=1, y'(0)=0.
        Vec yStart = new Vec(new double[] { 1.0, 0.0 });
        var interval = (start: 0.0, end: 10.0);

        ODEResult result = ODESolver.driver(harmonicOscillator, interval, yStart);
        List<double> xList = result.XList;
        List<Vec> yList = result.YList;

        using (StreamWriter sw = new StreamWriter("harmonic.txt"))
        {
            sw.WriteLine("x\tu\tv");
            for (int i = 0; i < xList.Count; i++)
                sw.WriteLine($"{xList[i]}\t{yList[i][0]}\t{yList[i][1]}");
        }
        Console.WriteLine("Harmonic oscillator data written to harmonic.txt");
    }

    static void RunFrictionHarmonicOscillator()
    {
        // Define the ODE for a damped pendulum:
        // theta' = y[1] and theta'' = -b*y[1] - c*sin(theta).
        Func<double, Vec, Vec> frictionHarmonicOscillator = (x, y) =>
        {   
            double b = 0.25;
            double c = 5.0;
            Vec dydx = new Vec(2);
            dydx[0] = y[1];                      
            dydx[1] = -b * y[1] - c * Sin(y[0]);
            return dydx;
        };

        // Initial conditions: theta(0)=PI/4, theta'(0)=0.
        Vec yStart = new Vec(new double[] { PI / 4, 0.0 });
        var interval = (start: 0.0, end: 10.0);

        ODEResult result = ODESolver.driver(frictionHarmonicOscillator, interval, yStart);
        List<double> xList = result.XList;
        List<Vec> yList = result.YList;

        using (StreamWriter sw = new StreamWriter("friction.txt"))
        {
            sw.WriteLine("x\t theta\t theta'");
            for (int i = 0; i < xList.Count; i++)
                sw.WriteLine($"{xList[i]}\t{yList[i][0]}\t{yList[i][1]}");
        }
        Console.WriteLine("Friction harmonic oscillator data written to friction.txt");
    }

    // Run three scenarios of planet orbit with GR and output them to a single file.
    static void RunPlanetOrbitGRMultiple()
    {
        var interval = (start: 0.0, end: 10.0);
        double hstart = 0.01, acc = 0.01, eps = 0.01;

        // Helper lambda: returns a planet orbit ODE function for a given gamma.
        Func<double, Func<double, Vec, Vec>> makePlanetOrbitGR = (gamma) =>
        {
            return (phi, y) =>
            {
                Vec dydphi = new Vec(2);
                dydphi[0] = y[1];
                dydphi[1] = 1 - y[0] + gamma * y[0] * y[0];
                return dydphi;
            };
        };

        // Scenario 1: gamma = 0, u(0)=1, u'(0)=0.
        double gamma1 = 0.0;
        Vec yStart1 = new Vec(new double[] { 1.0, 0.0 });
        Func<double, Vec> ip1 = Interpolator.MakeOdeIvpInterpolant(makePlanetOrbitGR(gamma1), interval, yStart1, hstart, acc, eps);

        // Scenario 2: gamma = 0, u(0)=1, u'(0)=-0.5.
        double gamma2 = 0.0;
        Vec yStart2 = new Vec(new double[] { 1.0, -0.5 });
        Func<double, Vec> ip2 = Interpolator.MakeOdeIvpInterpolant(makePlanetOrbitGR(gamma2), interval, yStart2, hstart, acc, eps);

        // Scenario 3: gamma = 0.01, u(0)=1, u'(0)=-0.5.
        double gamma3 = 0.05;
        Vec yStart3 = new Vec(new double[] { 1.0, -0.5 });
        Func<double, Vec> ip3 = Interpolator.MakeOdeIvpInterpolant(makePlanetOrbitGR(gamma3), interval, yStart3, hstart, acc, eps);

        // Write combined data to one file.
        using (StreamWriter sw = new StreamWriter("planet_orbit_gr.txt"))
        {
            // Header: phi, then u and u' for scenario 1, scenario 2, and scenario 3.
            sw.WriteLine("phi\t u1\t u1'\t u2\t u2'\t u3\t u3'");
            for (double phi = interval.start; phi <= interval.end; phi += 0.1)
            {
                Vec sol1 = ip1(phi);
                Vec sol2 = ip2(phi);
                Vec sol3 = ip3(phi);
                sw.WriteLine($"{phi}\t{sol1[0]}\t{sol1[1]}\t{sol2[0]}\t{sol2[1]}\t{sol3[0]}\t{sol3[1]}");
            }
        }
        Console.WriteLine("Multiple planet orbit GR scenarios written to planet_orbit_gr.txt");
    }
}