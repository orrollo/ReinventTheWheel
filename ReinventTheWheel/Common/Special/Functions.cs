using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Common.Special
{
    public static class Functions
    {
        public static double GammaLower05(double x)
        {
            return Math.Sqrt(Math.PI) * ErrorFunction.Erf(Math.Sqrt(x));
        }

        public static double GammaUpper05(double x)
        {
            return Math.Sqrt(Math.PI) * ErrorFunction.Erfc(Math.Sqrt(x));
        }

        // "A Computational Procedure for Incomplete Gamma Functions", Walter Gautschi
        // https://www.cs.purdue.edu/homes/wxg/selected_works/section_02/068.pdf
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
        public static double GammaLn(double z)
        {
            return z < 1 ? GammaLn(z + 1) - Math.Log(z) : Factorial.LuschnyCf4FactorialLn(z - 1);
        }

        public static double Gamma(double z)
        {
            return Math.Exp(GammaLn(z));
        }

        // "A Computational Procedure for Incomplete Gamma Functions", Walter Gautschi
        // https://www.cs.purdue.edu/homes/wxg/selected_works/section_02/068.pdf
        public static double GammaLowerStar(double a, double x)
        {
            return Kummer(a, a + 1, -x) / Math.Exp(GammaLn(a + 1));
        }

        // "A Computational Procedure for Incomplete Gamma Functions", Walter Gautschi
        // https://www.cs.purdue.edu/homes/wxg/selected_works/section_02/068.pdf
        public static double GammaLower(double a, double x)
        {
            return Math.Exp(GammaLn(a)) * GammaLowerStar(a, x) / Math.Pow(x, -a);
        }


    }
}
