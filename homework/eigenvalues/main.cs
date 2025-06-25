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
        double diag = -2.0 * (-0.5 / (dr * dr));
        double off  =  1.0 * (-0.5 / (dr * dr));

        for (int i = 0; i < npoints - 1; i++)
        {
            H[i, i]     = diag;
            H[i, i + 1] = off;
            H[i + 1, i] = off;
        }
        H[npoints - 1, npoints - 1] = diag;

        // Add the -1/r potential
        for (int i = 0; i < npoints; i++)
            H[i, i] += -1.0 / r[i];

        // Diagonalize via Jacobi
        (Vec e, Mat V) = Mat.cyclic(H);
        (Vec eSorted, Mat VSorted) = Mat.Sort(e, V);

        return eSorted[0];
    }

    static void Main(string[] args)
    {
        // Default parameters
        double rmax  = 15.0;
        double dr    = 0.2;
        string study = ""; // "dr", "rmax", or "wavefunction"

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

        // 3) Wavefunction study: numeric + analytic f₀, f₁, f₂
        if (study == "wavefunction")
        {
            int npoints = (int)(rmax / dr) - 1;
            if (npoints < 1)
            {
                Console.WriteLine("Invalid grid: increase rmax or decrease dr.");
                return;
            }

            // Build the radial grid
            Vec r = new Vec(npoints);
            for (int i = 0; i < npoints; i++)
                r[i] = dr * (i + 1);

            // Build Hamiltonian (same as in ComputeGroundState)
            Mat H = new Mat(npoints, npoints);
            double diag = -2.0 * (-0.5 / (dr * dr));
            double off  =  1.0 * (-0.5 / (dr * dr));
            for (int i = 0; i < npoints - 1; i++)
            {
                H[i, i]     = diag;
                H[i, i + 1] = off;
                H[i + 1, i] = off;
            }
            H[npoints - 1, npoints - 1] = diag;
            for (int i = 0; i < npoints; i++)
                H[i, i] += -1.0 / r[i];

            // Diagonalize
            (Vec e, Mat V) = Mat.cyclic(H);
            (Vec eSorted, Mat VSorted) = Mat.Sort(e, V);

            // Prepare storage
            int    N         = npoints;
            double normConst = 1.0 / Math.Sqrt(dr);

            var f0_num = new double[N];
            var f1_num = new double[N];
            var f2_num = new double[N];
            var f0_ana = new double[N];
            var f1_ana = new double[N];
            var f2_ana = new double[N];

            // Detect sign flips so numeric and analytic agree at small r
            bool flip0 = VSorted[0, 0] < 0.0;
            bool flip1 = VSorted[0, 1] < 0.0;
            bool flip2 = VSorted[0, 2] < 0.0;

            // Accumulate for discrete normalization of analytic
            double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;
            for (int i = 0; i < N; i++)
            {
                // numeric reduced wavefunctions
                f0_num[i] = (flip0 ? -1.0 : +1.0) * (VSorted[i, 0] * normConst);
                f1_num[i] = (flip1 ? -1.0 : +1.0) * (VSorted[i, 1] * normConst);
                f2_num[i] = (flip2 ? -1.0 : +1.0) * (VSorted[i, 2] * normConst);

                // continuous analytic forms
                // 1s: f₀ ∝ r e^{-r}
                f0_ana[i] = 2.0 * r[i] * Math.Exp(-r[i]);
                // 2s: f₁ ∝ (1 - r/2) r e^{-r/2}
                f1_ana[i] = (1.0 / Math.Sqrt(2.0)) 
                             * (1.0 - r[i] / 2.0) 
                             * r[i] 
                             * Math.Exp(-r[i] / 2.0);
                // 3s: f₂ ∝ (27 - 18r + 2r²) r e^{-r/3}
                double c3 = 1.0 / (81.0 * Math.Sqrt(3.0));
                f2_ana[i] = c3
                             * r[i]
                             * (27.0 - 18.0 * r[i] + 2.0 * r[i] * r[i])
                             * Math.Exp(-r[i] / 3.0);

                // accumulate ∑ f_ana²
                sum0 += f0_ana[i] * f0_ana[i];
                sum1 += f1_ana[i] * f1_ana[i];
                sum2 += f2_ana[i] * f2_ana[i];
            }

            // scale analytic so that Σ f_ana² Δr = 1
            double scale0 = 1.0 / Math.Sqrt(sum0 * dr);
            double scale1 = 1.0 / Math.Sqrt(sum1 * dr);
            double scale2 = 1.0 / Math.Sqrt(sum2 * dr);

            // Print header for Gnuplot
            Console.WriteLine("# i   r     f0_num    f0_ana    f1_num    f1_ana    f2_num    f2_ana");
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine(
                    $"{i,3} {r[i],6:F3}  " +
                    $"{f0_num[i],8:F5}  " +
                    $"{f0_ana[i] * scale0,8:F5}  " +
                    $"{f1_num[i],8:F5}  " +
                    $"{f1_ana[i] * scale1,8:F5}  " +
                    $"{f2_num[i],8:F5}  " +
                    $"{f2_ana[i] * scale2,8:F5}"
                );
            }

            // Print the first three eigenvalues
            Console.WriteLine("# Eigenvalues:");
            for (int k = 0; k < 3; k++)
                Console.WriteLine($"# e[{k}] = {eSorted[k]:F8}");

            return;
        }

        // 4) Exercise A: test Jacobi diagonalization
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

        Mat diagMatrix      = eigenvecs.Transpose() * symMatrix * eigenvecs;
        Console.WriteLine("Diagonal matrix:");
        diagMatrix.Print();

        Mat originalMatrix  = eigenvecs * diagMatrix * eigenvecs.Transpose();
        Console.WriteLine("Reconstructed matrix:");
        originalMatrix.Print();

        Mat identityMatrix  = eigenvecs * eigenvecs.Transpose();
        Console.WriteLine("V*V^T:");
        identityMatrix.Print();

        identityMatrix      = eigenvecs.Transpose() * eigenvecs;
        Console.WriteLine("V^T*V:");
        identityMatrix.Print();
    }
}
