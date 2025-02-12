using static System.Math;
using static System.Console;

public class Program {
public static void Main(string[] args) {
    foreach(var arg in args) {
        var words = arg.Split(':');   // Split each argument on the colon.
        if (words[0] == "-numbers") {  // Check if the option is "-numbers".
            var numbers = words[1].Split(',');  // Split the list of numbers by commas.
            WriteLine("x \t sin(x) \t cos(x)");  // Print the header.
            foreach(var number in numbers) {
                double x = double.Parse(number); // Convert the string to a double.
                WriteLine($"{x:F1} \t {Sin(x):F3} \t {Cos(x):F3}");  // Print the number, its sine, and its cosine.
                }
            }
        }
    }
}