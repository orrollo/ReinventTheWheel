using System;

namespace Common.Special
{
    public static class Factorial
    {
        // luschnyCF4; http://www.luschny.de/math/factorial/approx/SimpleCases.html
        public static double LuschnyCf4FactorialLn(double n)
        {
            double[] c = new[] { 1.0 / 24, 3.0 / 80, 18029.0 / 45360, 6272051.0 / 14869008 };
            double N = n + 0.5;
            double p = (N * N) / (N + c[0] / (N + c[1] / (N + c[2] / (N + c[3] / N))));
            double logFact = Math.Log(2 * Math.PI) / 2 + N * (Math.Log(p) - 1);
            return logFact;
        }

        public static double StieltjesLnFactorial(double z)
        {
            double[] a = new[] { 1.0/12, 1.0/30, 53.0/210, 195.0/371, 22999.0/22737, 29944523.0/19733142, 109535241009.0/48264275462 };
            double zz = z + 1, cc = 0.0;
            for (int i = 6; i >= 0; i--) cc = zz + a[i] / cc;
            return Math.Log(2 * Math.PI) / 2.0 + (zz - 0.5) * Math.Log(zz) - zz + cc;
        }

        public static double StieltjesFactorial(double x)
        {
            double y = x, p = 1;
            while (y < 8) { p = p * y; y = y + 1; }
            double r = Math.Exp(StieltjesLnFactorial(y));
            if (x < 8) r = (r * x) / (p * y);
            return r;
        }
    }
}