using System;
using doublecompare; // Reference the external library

class Program
{
    static void Main()
    {
        // Find max integer
        int maxInt = 1;
        while (maxInt + 1 > maxInt)
        {
            maxInt++;
        }
        
        int minInt = -1;
        while (minInt - 1 < minInt)
        {
            minInt--;
        }
        
        // Find min integer
        double epsilonDouble = 1.0;
        while (1.0 + epsilonDouble != 1.0)
        {
            epsilonDouble /= 2.0;
        }
        epsilonDouble *= 2.0;
        
        // Find the lowest value that when summed with 1 changes its value
        float epsilonFloat = 1.0F;
        while ((float)(1.0F + epsilonFloat) != 1.0F)
        {
            epsilonFloat /= 2.0F;
        }
        epsilonFloat *= 2.0F;
        
        // Print results
        Console.WriteLine("My max int = {0}", maxInt);
        Console.WriteLine("int.MaxValue = {0}", int.MaxValue);
        Console.WriteLine("Match: {0}", maxInt == int.MaxValue);
        
        Console.WriteLine("My min int = {0}", minInt);
        Console.WriteLine("int.MinValue = {0}", int.MinValue);
        Console.WriteLine("Match: {0}", minInt == int.MinValue);
        
        Console.WriteLine("Machine epsilon for double: {0}", epsilonDouble);
        Console.WriteLine("Machine epsilon for float: {0}", epsilonFloat);
        
        Console.WriteLine("2^-52: {0}", Math.Pow(2, -52));
        Console.WriteLine("2^-23: {0}", Math.Pow(2, -23));

        // Make a smaller value than machine epsilon and see how it affects arithmetic
        double epsilon = Math.Pow(2, -52);
        double tiny = epsilon / 2;
        double a = 1 + tiny + tiny;
        double b = tiny + tiny + 1;
        Console.WriteLine($"a==b ? {a == b}");
        Console.WriteLine($"a>1  ? {a > 1}");
        Console.WriteLine($"b>1  ? {b > 1}");

        // Compare two floating point numbers
        double d1 = 0.1 + 0.1 + 0.1 + 0.1;
        double d2 = 4 * 0.1;

        // Try comparing d1 and d2
        Console.WriteLine($"d1==d2 ? {DoubleCompare.Approx(d1, d2)}");
    }
}