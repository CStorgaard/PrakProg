using System;

static class program{

    static void Main(){
    // Find square root of 2
    double res0 = Math.Sqrt(2);
    Console.WriteLine($"Square root of 2 is {res0}");

    // Find 2 to the 1/5 power
    double res1 = Math.Pow(2.0, 0.2);
    Console.WriteLine($"2 to the one-fifth power is {res1}");

    // Find e^pi
    double res2 = Math.Pow(Math.E, Math.PI);
    Console.WriteLine($"e to the pi'th power is {res2}");

    // Find pi^e
    double res3 = Math.Pow(Math.PI, Math.E);
    Console.WriteLine($"pi to the e'th power is {res3}");

    for (int i=1; i<= 10; i++){
        double res4 = sfuns.fgamma(i);
        Console.WriteLine($"Value of gamma function with input {i} is {res4}");
    }
for (int i=1; i<= 10; i++){
        double res4 = sfuns.lnfgamma(i);
        Console.WriteLine($"Value of log gamma function with input {i} is {res4}");
    }
    }
}