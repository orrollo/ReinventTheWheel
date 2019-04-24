using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Common.Special;

namespace Common
{
    public class ChiSquare
    {
        public static double Pdf(double k, double x)
        {
            double k2 = k / 2, v = (k2 - 1.0) * Math.Log(x) - (x / 2.0) - k2 * Math.Log(2) - Functions.GammaLn(k2);
            return Math.Exp(v);
        }

        public static double Cdf(double k, double x)
        {
            return Functions.GammaLower(k / 2, x / 2) / Functions.Gamma(k / 2);
        }

        public static double Quantile(double k, double p)
        {
            double lower = 1e-30, upper = 1;
            while (Cdf(k, upper) < p) upper *= 2;
            while ((upper - lower) > 1e-5)
            {
                var middle = (lower + upper) / 2.0;
                double value = Cdf(k, middle);
                if (Math.Abs(value - p) < 1e-10) return middle;
                if (value < p) lower = middle;
                else upper = middle;
            }

            return (lower + upper) / 2.0;
        }
    }
}
