using System.Numerics;
using System;

public class Program{
    public static void Main(string[] args) {
        Complex a = Complex.Sqrt(-1);
        Console.WriteLine($"Square root of -1 is {a}");
        Complex b = Complex.Log(Complex.ImaginaryOne);
        Console.WriteLine($"Natural logarithm of i is {b}");
        Complex c = Complex.Sqrt(Complex.ImaginaryOne);
        Console.WriteLine($"Square root of i is {c}");
        Complex d = Complex.Pow(Complex.ImaginaryOne, Complex.ImaginaryOne);
        Console.WriteLine($"i to the power of i is {d}");
        Complex e = Complex.Sinh(Complex.ImaginaryOne);
        Console.WriteLine($"Hyperbolic sine of i is {e}");
        Complex f = Complex.Cosh(Complex.ImaginaryOne);
        Console.WriteLine($"Hyperbolic cosine of i is {f}");
    }
}