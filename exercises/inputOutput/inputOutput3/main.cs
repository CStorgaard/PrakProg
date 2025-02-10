using static System.Math;   // import static members of System.Math
using static System.Console; // import static members of System.Console

class main {
    static int Main(string[] args) {
        string infile = null, outfile = null; // input/output file names

        // Parse command-line arguments
        foreach (var arg in args) {
            var words = arg.Split(':');
            if (words[0] == "-input")
                infile = words[1]; // input file name
            if (words[0] == "-output")
                outfile = words[1]; // output file name
        }

        // Check if input/output file names are given
        if (infile == null || outfile == null) {
            Error.WriteLine("wrong filename argument");
            return 1;
        }

        // Open input and output file streams
        var instream = new System.IO.StreamReader(infile);
        var outstream = new System.IO.StreamWriter(outfile, append: false);

        // Write header to output file
        outstream.WriteLine("x \t sin(x) \t cos(x)");

        // Read input file line by line
        for (string line = instream.ReadLine(); line != null; line = instream.ReadLine()) {
            string trimmed = line.Trim();

            // Check if the line starts with "-numbers:"
            if (trimmed.StartsWith("-numbers:")) {
                // Extract the comma-separated numbers after "-numbers:"
                string numbersPart = trimmed.Substring("-numbers:".Length);
                var numberStrings = numbersPart.Split(',');

                foreach (var ns in numberStrings) {
                    double x;
                    if (double.TryParse(ns.Trim(), out x)) {
                        outstream.WriteLine($"{x:F1} \t {Sin(x):F3} \t {Cos(x):F3}");
                    } else {
                        Error.WriteLine($"Warning: Could not parse '{ns}' as a number. Skipping.");
                    }
                }
            }
            else {
                // Try to parse the line as a single number
                double x;
                if (double.TryParse(trimmed, out x)) {
                    outstream.WriteLine($"{x:F1} \t {Sin(x):F3} \t {Cos(x):F3}");
                } else {
                    Error.WriteLine($"Warning: Could not parse '{line}' as a number. Skipping.");
                }
            }
        }

        // Close file streams
        instream.Close();
        outstream.Close();

        return 0; // return success
    }
}
