using System;

class Program
{
    /// <summary>
    /// Builds the Hamiltonian for the s-wave hydrogen problem with given rmax, dr,
    /// diagonalizes it, and returns the lowest eigenvalue (ground-state energy).
    /// </summary>
    static double ComputeGroundState(double rmax, double dr)
    {
        int npoints = (int)(rmax / dr) - 1;
        if (npoints < 1) return double.NaN;

        // Construct the radial grid
        Vec r = new Vec(npoints);
        for (int i = 0; i < npoints; i++)
            r[i] = dr * (i + 1);

        // Build the Hamiltonian matrix
        Mat H = new Mat(npoints, npoints);
        for (int i = 0; i < npoints - 1; i++)
        {
            H[i, i]     = -2 * (-0.5 / (dr * dr));
            H[i, i + 1] =  1 * (-0.5 / (dr * dr));
            H[i + 1, i] =  1 * (-0.5 / (dr * dr));
        }
        H[npoints - 1, npoints - 1] = -2 * (-0.5 / (dr * dr));

        // Add the -1/r potential
        for (int i = 0; i < npoints; i++)
            H[i, i] += -1.0 / r[i];

        // Diagonalize via Jacobi
        (Vec e, Mat V) = Mat.cyclic(H);
        (Vec eSorted, Mat _) = Mat.Sort(e, V);

        return eSorted[0];
    }

    static void Main(string[] args)
    {
        // Default parameters
        double rmax = 10.0;
        double dr   = 0.3;
        string study = ""; // possible values: "dr", "rmax", "wavefunction"

        // Parse command-line arguments
        foreach (var arg in args)
        {
            if (arg.StartsWith("-rmax:"))
                rmax = double.Parse(arg.Substring("-rmax:".Length));
            else if (arg.StartsWith("-dr:"))
                dr = double.Parse(arg.Substring("-dr:".Length));
            else if (arg.StartsWith("-study:"))
                study = arg.Substring("-study:".Length).ToLower();
        }

        // 1) Δr study
        if (study == "dr")
        {
            Console.WriteLine("# dr   epsilon0");
            for (double drVar = 0.05; drVar <= 0.5; drVar += 0.01)
            {
                double e0 = ComputeGroundState(rmax, drVar);
                Console.WriteLine($"{drVar} {e0}");
            }
            return;
        }

        // 2) rmax study
        if (study == "rmax")
        {
            Console.WriteLine("# rmax   epsilon0");
            for (double rmaxVar = 5.0; rmaxVar <= 20.0; rmaxVar += 1.0)
            {
                double e0 = ComputeGroundState(rmaxVar, dr);
                Console.WriteLine($"{rmaxVar} {e0}");
            }
            return;
        }

        // 3) Wavefunction study: numeric + analytic
        if (study == "wavefunction")
        {
            int npoints = (int)(rmax / dr) - 1;
            if (npoints < 1)
            {
                Console.WriteLine("Invalid grid: increase rmax or decrease dr.");
                return;
            }

            // Build the grid
            Vec r = new Vec(npoints);
            for (int i = 0; i < npoints; i++)
                r[i] = dr * (i + 1);

            // Build Hamiltonian
            Mat H = new Mat(npoints, npoints);
            for (int i = 0; i < npoints - 1; i++)
            {
                H[i, i]     = -2 * (-0.5 / (dr * dr));
                H[i, i + 1] =  1 * (-0.5 / (dr * dr));
                H[i + 1, i] =  1 * (-0.5 / (dr * dr));
            }
            H[npoints - 1, npoints - 1] = -2 * (-0.5 / (dr * dr));
            for (int i = 0; i < npoints; i++)
                H[i, i] += -1.0 / r[i];

            // Diagonalize
            (Vec e, Mat V) = Mat.cyclic(H);
            (Vec eSorted, Mat VSorted) = Mat.Sort(e, V);

            // We print 3 states: n=1,2,3
            int statesToPlot = 2;
            double normConst = 1.0 / Math.Sqrt(dr); // factor for reduced radial wavefunctions

            // Print a header line for gnuplot
            Console.WriteLine("# i   r[i]    f0(r[i])  f0_analytic   f1(r[i])  f1_analytic");

            for (int i = 0; i < npoints; i++)
            {
                // Numeric wavefunctions
                double f0_num = normConst * VSorted[i, 0]; // n=1
                double f1_num = normConst * VSorted[i, 1]; // n=2
                // double f2_num = normConst * VSorted[i, 2]; // n=3

                // Analytic wavefunctions
                // n=1
                double f0_ana = 2.0 * r[i] * Math.Exp(-r[i]);
                // n=2
                double f1_ana = (1.0 / Math.Sqrt(2.0)) * (1 - r[i] / 2.0) * r[i] * Math.Exp(-r[i] / 2.0);
                // n=3 (the 3s state)
                // double f2_ana = (2.0 / (27.0 * Math.Sqrt(3.0)))
                // * r[i]
                // * (27.0 - 18.0*r[i] + 2.0*r[i]*r[i])
                // * Math.Exp(-r[i] / 3.0);


                Console.WriteLine($"{i} {r[i]} {f0_num} {f0_ana} {f1_num} {f1_ana}");
            }

            // Optionally print eigenvalues as comments
            Console.WriteLine("# Eigenvalues:");
            for (int k = 0; k < statesToPlot; k++)
            {
                Console.WriteLine($"# e[{k}] = {eSorted[k]}");
            }
            return;
        }

        // 4) Exercise A
        Console.WriteLine("Exercise A");
        Console.WriteLine("Checking Jacobi diagonalization");

        Mat symMatrix = new Mat(2, 2);
        symMatrix.SetElement(0, 0, 1.0);
        symMatrix.SetElement(0, 1, 2.0);
        symMatrix.SetElement(1, 0, 2.0);
        symMatrix.SetElement(1, 1, 4.0);

        Console.WriteLine("Input Matrix:");
        symMatrix.Print();

        (Vec eigenvals, Mat eigenvecs) = Mat.cyclic(symMatrix);

        Console.WriteLine("Eigenvalues:");
        eigenvals.Print();

        Console.WriteLine("Normalized eigenvectors:");
        eigenvecs.Print();

        Mat diagMatrix = eigenvecs.Transpose() * symMatrix * eigenvecs;
        Console.WriteLine("Diagonal matrix:");
        diagMatrix.Print();

        Mat originalMatrix = eigenvecs * diagMatrix * eigenvecs.Transpose();
        Console.WriteLine("Reconstructed matrix:");
        originalMatrix.Print();

        Mat identityMatrix = eigenvecs * eigenvecs.Transpose();
        Console.WriteLine("V*V^T:");
        identityMatrix.Print();

        identityMatrix = eigenvecs.Transpose() * eigenvecs;
        Console.WriteLine("V^T*V:");
        identityMatrix.Print();
    }
}