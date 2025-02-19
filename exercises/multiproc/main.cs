using System;
using System.Linq;            // Needed for .Sum()
using System.Threading;
using System.Threading.Tasks;

class main {
    // Data structure for the manual threading version.
    public class datum {
        public int start, stop;
        public double sum;
    }

    // Thread function for the manual threading version.
    public static void harm(object obj) {
        datum d = (datum)obj;
        d.sum = 0;
        // Compute partial harmonic sum for the assigned range.
        for (int i = d.start + 1; i <= d.stop; i++) {
            d.sum += 1.0 / i;
        }
    }

    // Compute the harmonic sum using a flawed Parallel.For version
    // that directly updates a shared variable (leading to race conditions).
    public static void HarmonicParallel(int N) {
        double sum = 0;
        Parallel.For(1, N + 1, (int i) => sum += 1.0 / i);
        Console.WriteLine("Flawed Parallel.For computed sum = " + sum);
    }

    // Compute the harmonic sum using ThreadLocal to safely accumulate
    // partial results from each thread, then summing them up.
    public static void HarmonicThreadLocal(int N) {
        var sum = new ThreadLocal<double>(() => 0, trackAllValues: true);
        Parallel.For(1, N + 1, (int i) => sum.Value += 1.0 / i);
        double totalSum = sum.Values.Sum();
        Console.WriteLine("ThreadLocal Parallel.For computed sum = " + totalSum);
    }

    public static int Main(string[] argv) {
        // Default parameters.
        int nterms = 100000000; // 1e8 terms by default.
        int nthreads = 1;       // 1 thread by default.
        bool useParallel = false;      // Flawed parallel mode.
        bool useThreadLocal = false;   // Safe thread-local parallel mode.

        // Process command-line arguments.
        for (int i = 0; i < argv.Length; i++) {
            string arg = argv[i];
            if (arg == "-threads" && i + 1 < argv.Length) {
                nthreads = int.Parse(argv[++i]);
            } else if (arg == "-terms" && i + 1 < argv.Length) {
                nterms = int.Parse(argv[++i]);
            } else if (arg == "-parallel") {
                useParallel = true;
            } else if (arg == "-threadlocal") {
                useThreadLocal = true;
            }
        }

        // If -threadlocal is specified, run the safe ThreadLocal version.
        if (useThreadLocal) {
            Console.WriteLine("Computing harmonic sum using Parallel.For with ThreadLocal:");
            HarmonicThreadLocal(nterms);
        }
        // Else if -parallel is specified, run the flawed Parallel.For version.
        else if (useParallel) {
            Console.WriteLine("Computing harmonic sum using flawed Parallel.For (shared variable race condition):");
            HarmonicParallel(nterms);
        }
        // Otherwise, use manual threading.
        else {
            Console.WriteLine("Computing harmonic sum using manual threading:");
            Thread[] threads = new Thread[nthreads];
            datum[] data = new datum[nthreads];

            // Initialize the data objects and assign each its portion of the work.
            for (int i = 0; i < nthreads; i++) {
                data[i] = new datum();
                data[i].start = i * nterms / nthreads;
                data[i].stop = (i + 1) * nterms / nthreads;
            }

            // Create and start each thread.
            for (int i = 0; i < nthreads; i++) {
                threads[i] = new Thread(harm);
                threads[i].Start(data[i]);
            }

            // Wait for all threads to finish.
            foreach (Thread thread in threads)
                thread.Join();

            // Combine the partial sums.
            double sum = 0;
            foreach (datum d in data)
                sum += d.sum;

            Console.WriteLine("Manual threading computed sum = " + sum);
        }

        return 0;
    }
}
