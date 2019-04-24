using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Common.Special
{
    public static class Functions
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

        public static double GammaLower05(double x)
        {
            return Math.Sqrt(Math.PI) * Erf(Math.Sqrt(x));
        }

        public static double GammaUpper05(double x)
        {
            return Math.Sqrt(Math.PI) * Erfc(Math.Sqrt(x));
        }

        public static double Kummer(double a, double b, double z)
        {
            double t = 1.0, sum = 0.0;

            for (int k = 0; k < 200; k++)
            {
                sum += t;
                t *= a + k;
                t /= b + k;
                t *= z;
                t /= (k + 1);
                if (Math.Abs(t) < 1e-40) break;
            }

            return sum;
        }

        // http://www.luschny.de/math/factorial/approx/SimpleCases.html
        public static double FactLn(double n)
        {
            double[] c = new[] { 1.0 / 24, 3.0 / 80, 18029.0 / 45360, 6272051.0 / 14869008 };
            double N = n + 0.5;
            double p = (N * N) / (N + c[0] / (N + c[1] / (N + c[2] / (N + c[3] / N))));
            double logFact = Math.Log(2 * Math.PI) / 2 + N * (Math.Log(p) - 1);
            return logFact;
        }

        public static double GammaLn(double z)
        {
            return z < 1 ? GammaLn(z + 1) - Math.Log(z) : FactLn(z - 1);
        }

        public static double Gamma(double z)
        {
            return Math.Exp(GammaLn(z));
        }

        public static double GammaLowerStar(double a, double x)
        {
            return Kummer(a, a + 1, -x) / Math.Exp(GammaLn(a + 1));
        }

        public static double GammaLower(double a, double x)
        {
            return Math.Exp(GammaLn(a)) * GammaLowerStar(a, x) / Math.Pow(x, -a);
        }


    }
}
