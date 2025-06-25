using System;
using static System.Math;
using System.IO;

class main
{
    static void Main()
    {
        // 1) Define your data (time, activity, uncertainty):
        double[] tData  = { 1,  2,  3,  4,  6,   9, 10, 13, 15 };
        double[] yData  = {117,100, 88, 72, 53,  29.5, 25.2, 15.2, 11.1 };
        double[] dyData = {  6,  5,  4,  4,  4,   3,  3,  2,  2 };

        int n = tData.Length;

        // 2) Transform y -> ln(y) and compute uncertainties d(ln(y)) = dy / y:
        double[] lnYData   = new double[n];
        double[] dlnYData  = new double[n];
        for (int i = 0; i < n; i++)
        {
            lnYData[i]  = Log(yData[i]);
            dlnYData[i] = dyData[i] / yData[i];
        }

        // 3) Create Vec objects for t, ln(y), and d(ln(y)):
        Vec x    = new Vec(tData);
        Vec lnY  = new Vec(lnYData);
        Vec dlnY = new Vec(dlnYData);

        // 4) Define the basis functions for the fit:
        //    ln(y) = c0 + c1*t  =>  c0=ln(a), c1=-lambda
        Func<double,double>[] fs = {
            (double t) => 1.0,      // c0 part
            (double t) => t        // c1 * t part (we will interpret c1 as negative lambda)
        };

        // 5) Perform the weighted least squares fit:
        //    LsFit(fs, x, lnY, dlnY) returns (c, Cov)
        //    c is the best-fit coefficients, Cov is the covariance matrix
        var (c, Cov) = Mat.LsFit(fs, x, lnY, dlnY);

        // 6) Extract the parameters:
        double c0       = c[0];         // c0 = ln(a)
        double c1       = c[1];         // c1 = -lambda
        double a        = Exp(c0);      // a = exp(c0)
        double lambda   = -c1;          // lambda = -c1

        // Find the uncertainty in the half-life:
        double dlambda = Sqrt(Cov[1,1]);
        double dHalfLife = Log(2) * dlambda / lambda / lambda;

        // 7) Print the results:
        Console.WriteLine("Fitted parameters:");
        Console.WriteLine($"  ln(a)   = {c0:F3}");
        Console.WriteLine($"  a       = {a:F3}");
        Console.WriteLine($"  lambda  = {lambda:F3}");
        Console.WriteLine($"  half-life = {Log(2)/lambda:F3} days");
        Console.WriteLine();
        Console.WriteLine("Covariance matrix (Cov):");
        Cov.Print();
        Console.WriteLine();
        Console.WriteLine($"Calculated uncertainty in half-life: {dHalfLife:F3} days");
        Console.WriteLine($"Modern uncertainty in half-life: 0.0014 days");



        // -------------------------------
        // 1) Write out data points to a file (for error bars in Gnuplot).
        //    We'll plot in normal y-scale, so just store (t, y, dy).
        // -------------------------------
        using (var dataFile = new StreamWriter("dataPoints.txt"))
        {
            for (int i = 0; i < n; i++)
            {
                dataFile.WriteLine($"{tData[i]} {yData[i]} {dyData[i]}");
            }
        }

        // -------------------------------
        // 2) Write out a smooth set of fit points, 
        //    i.e. t from min to max in small steps, 
        //    yFit(t) = exp(c0 + c1 * t).
        // -------------------------------
        using (var fitFile = new StreamWriter("fitPoints.txt"))
        {
            double tMin = tData[0];
            double tMax = tData[n-1];
            int steps   = 100;   // Number of points in the smooth curve
            double dt    = (tMax - tMin)/(steps-1);
            for(int i = 0; i < steps; i++)
            {
                double T = tMin + i*dt;
                double lnFit = c0 + c1*T;
                double yFit  = Exp(lnFit);
                fitFile.WriteLine($"{T} {yFit}");
            }
        }

        using (var envelopeFile = new StreamWriter("fitEnvelope.txt"))
            {
            double tMin = tData[0];
            double tMax = tData[n - 1];
            int    steps = 100;
            double dt    = (tMax - tMin) / (steps - 1);

            for (int i = 0; i < steps; i++)
            {
                double T     = tMin + i * dt;
                double lnFit = c0 + c1 * T;            // best-fit in log-space
                double bestY = Math.Exp(lnFit);        

                // Variance *in ln(y)* = Cov[0,0] + 2*T*Cov[0,1] + (T^2)*Cov[1,1]
                double varLn   = Cov[0,0]
                            + 2.0 * T * Cov[0,1]
                            + (T * T) * Cov[1,1];
                double sigmaLn = Math.Sqrt(varLn);

                // Now shift *in ln-space*, then exponentiate
                double yPlus   = Math.Exp(lnFit + sigmaLn);
                double yMinus  = Math.Exp(lnFit - sigmaLn);

                // Columns: t, best-fit y, upper band, lower band
                envelopeFile.WriteLine(
                    $"{T} {bestY} {yPlus} {yMinus}"
                );
            }
        }
    }
}
