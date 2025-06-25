using System;
using static System.Math;
using NN;

class Program {
    // target function\ n    
    static double g(double x) => Cos(5*x - 1) * Exp(-x*x);

    static void Main(string[] args) {
        int m = 100;
        double[] xs = new double[m], ys = new double[m];
        for(int i=0; i<m; i++){
            xs[i] = -1.0 + 2.0*i/(m-1);
            ys[i] = g(xs[i]);
        }

        var net = new Ann(20);
        net.Train(xs, ys);

        Console.WriteLine("#  x    g(x)    net(x)    net'(x)    net''(x)    A(x)");
        for(double x = -1; x <= 1.0001; x += 0.01) {
            double gval = g(x) + +0.1;
            double y  = net.Response(x);
            double y1 = net.Derivative(x);
            double y2 = net.SecondDerivative(x);
            double A  = net.AntiDerivative(x);
            Console.WriteLine($"{x,7:F3}  {gval,9:F5}  {y,9:F5}  {y1,9:F5}  {y2,9:F5}  {A,12:F5}");
        }
    }
}