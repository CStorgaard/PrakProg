using System;
using static System.Math;

namespace NN {
  public class Ann {
    int n;
    Vec p;

    // Activation function and its derivatives and integral
        double f(double z)    => z * Exp(-0.5*z*z);
    double df(double z)   => (1 - z*z) * Exp(-0.5*z*z);
    double d2f(double z)  => (z*z*z - 3*z) * Exp(-0.5*z*z);
    double F(double z)    => -Exp(-0.5*z*z);



    public Ann(int hiddenNeurons) {
      n = hiddenNeurons;
      p = new Vec(3*n + 1);
      var rnd = new Random();
      for(int i=0; i<p.Length; i++) p[i] = (rnd.NextDouble()-0.5)*0.1;
    }

    // Forward pass
    public double Response(double x) {
      double sum = 0;
      for(int j=0; j<n; j++){
        double z = p[j]*x + p[n+j];
        sum += p[2*n+j] * f(z);
      }
      return sum + p[3*n];
    }

    // 1st derivative wrt x
    public double Derivative(double x) {
      double sum = 0;
      for(int j=0; j<n; j++){
        double w1 = p[j], b1 = p[n+j], w2 = p[2*n+j];
        double z = w1*x + b1;
        sum += w2 * df(z) * w1;
      }
      return sum;
    }

    // 2nd derivative wrt x
    public double SecondDerivative(double x) {
      double sum = 0;
      for(int j=0; j<n; j++){
        double w1 = p[j], b1 = p[n+j], w2 = p[2*n+j];
        double z = w1*x + b1;
        sum += w2 * d2f(z) * w1*w1;
      }
      return sum;
    }

    // Indefinite integral (constant of integration omitted)
    public double AntiDerivative(double x) {
      double sum = 0;
      for(int j=0; j<n; j++){
        double w1 = p[j], b1 = p[n+j], w2 = p[2*n+j];
        double z = w1*x + b1;
        sum += w2 * (F(z)/w1);
      }
      sum += p[3*n] * x;
      return sum;
    }

    // Analytic gradient for training
    Vec AnalyticGradient(double[] xs, double[] ys, Vec pTry) {
      var old = p; p = pTry;
      int m = xs.Length;
      Vec g = new Vec(p.Length);
      int oW1=0, oB1=n, oW2=2*n, oB2=3*n;
      for(int i=0;i<m;i++){
        double x = xs[i], y = ys[i];
        double[] z = new double[n], a = new double[n];
        for(int j=0; j<n; j++){ z[j] = p[oW1+j]*x + p[oB1+j]; a[j] = f(z[j]); }
        double ypred = p[oB2];
        for(int j=0; j<n; j++) ypred += p[oW2+j]*a[j];
        double r = ypred - y;
        for(int j=0; j<n; j++) g[oW2+j] += r * a[j];
        g[oB2] += r;
        for(int j=0; j<n; j++){
          double w2 = p[oW2+j];
          g[oW1+j] += r * w2 * df(z[j]) * xs[i];
          g[oB1+j] += r * w2 * df(z[j]);
        }
      }
      p = old;
      return g;
    }

    // Training with gradient descent
    public void Train(double[] xs, double[] ys,
                      double tol=1e-8, int maxIter=50000) {
      if(xs.Length != ys.Length) throw new ArgumentException();
      Func<Vec,double> phi = param => {
        var old = p; p = param;
        double s = 0;
        for(int i=0; i<xs.Length; i++){
          double d = Response(xs[i]) - ys[i]; s += d*d;
        }
        p = old; return 0.5*s;
      };
      Func<Vec,Vec> gphi = param => AnalyticGradient(xs, ys, param);
      var p0 = p.Copy();

      p = GradientDescent.Minimize(
        phi, gphi, p0,
        step:    1e-1,
        maxIter: maxIter,
        tol:     tol,
        beta:    0.9
      );

      Console.WriteLine(
        $"Training done in {GradientDescent.Iterations} iterations, SSE={2*phi(p):E}"
      );
    }
  }
}