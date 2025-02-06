using System;

namespace doublecompare
{
    public static class DoubleCompare
    {
        public static bool Approx(double a, double b, double acc = 1e-9, double eps = 1e-9)
        {
            if (Math.Abs(b - a) <= acc) return true;
            if (Math.Abs(b - a) <= Math.Max(Math.Abs(a), Math.Abs(b)) * eps) return true;
            return false;
        }
    }
}
