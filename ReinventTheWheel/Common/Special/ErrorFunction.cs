using System;

namespace Common.Special
{
    public static class ErrorFunction
    {
        public static double Erf(double x)
        {
            double a = 8.0 * (Math.PI - 3) / (3 * Math.PI * (4 - Math.PI));
            double v1 = 4.0 / Math.PI + a * x * x;
            double v2 = 1 + a * x * x;
            double v3 = -x * x * v1 / v2;
            double v4 = Math.Sqrt(1.0 - Math.Exp(v3));
            return v4;
        }

        public static double Erfc(double x)
        {
            return 1.0 - Erf(x);
        }
    }
}