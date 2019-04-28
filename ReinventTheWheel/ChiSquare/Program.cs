using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Common.Special;

namespace ChiSquare
{
    class Program
    {
        // "A Computational Procedure for Incomplete Gamma Functions", Walter Gautschi
        // https://www.cs.purdue.edu/homes/wxg/selected_works/section_02/068.pdf
        static double Kummer(double a, double b, double z)
        {
            double t = 1.0, sum = 0.0;

            for (int k = 0; k < 100; k++)
            {
                sum += t;
                t *= a + k;
                t /= b + k;
                t *= z;
                t /= (k + 1);
                if (Math.Abs(t) < 1e-30) break;
            }

            return sum;
        }

        // "A Computational Procedure for Incomplete Gamma Functions", Walter Gautschi
        // https://www.cs.purdue.edu/homes/wxg/selected_works/section_02/068.pdf
        static double GammaLowerStar(double a, double x)
        {
            return Kummer(a, a + 1, -x) / Math.Exp(GammaLn2(a + 1));
        }

        static double GammaLower4(double a, double x)
        {
            return Math.Exp(GammaLn2(a)) * GammaLowerStar(a, x) / Math.Pow(x, -a);
        }

        static double Quantile1(double alfa, int n)
        {
            double d, A, B, C, D, E;

            if (alfa < 0.5)
            {
                d = -2.0637 * Math.Pow((Math.Log(1 / alfa) - 0.16), 0.4274) + 1.5774;
            }
            else
            {
                d = 2.0637 * Math.Pow((Math.Log(1 / (1.0-alfa)) - 0.16), 0.4274) - 1.5774;
            }

            A = d * Math.Sqrt(2);
            B = 2.0 * (d * d - 1) / 3.0;
            C = d * (d * d - 7) / (9 * Math.Sqrt(2));
            D = -(6 * d * d * d * d + 14 * d * d - 32) / 405;
            E = d * (9 * d * d * d * d + 256 * d * d - 433) / (4860 * Math.Sqrt(2));

            double sn = Math.Sqrt(n);

            return n + A * sn + B + C / sn + D / n + E / (n * sn);
        }

        static double Quantile2(double alfa, int n)
        {
            double[] a = new[] { 1.0000886, 0.4713941, 0.0001348028, -0.008553069, 0.00312558, -0.0008426812, 0.00009780499 };
            double[] b = new[] { -0.2237368,0.02607083,0.01128186,-0.01153761,0.005169654,0.00253001,-0.001450117 };
            double[] c = new[] { -0.01513904,-0.008986007,0.02277679,-0.01323293,-0.006950356,0.001060438,0.001565326 };
            double d, result = 0.0;

            if (alfa < 0.5)
            {
                d = -2.0637 * Math.Pow((Math.Log(1 / alfa) - 0.16), 0.4274) + 1.5774;
            }
            else
            {
                d = 2.0637 * Math.Pow((Math.Log(1 / (1.0 - alfa)) - 0.16), 0.4274) - 1.5774;
            }


            for (int i = 0; i <= 6; i++) result += Math.Pow(n, -i / 2.0) * Math.Pow(d, i) * (a[i] + b[i] / n + c[i] / (n * n));
            return n * Math.Pow(result, 3);
        }

        static void check(Func<double, int, double> fn, double value, int n, double result, double border = 0.05)
        {
            var check = fn(value, n);
            var error = Math.Abs(result - check) / result;
            if (error > border) throw new Exception("unexact solution");
        }

        static double GammaLn(double z)
        {
            double emc = 0.5772156649, result = -emc*z - Math.Log(z);
            for (int n = 1; n < 10000; n++)
            {
                double relation = z / n;
                double value = relation - Math.Log(1.0 + relation);
                result += value;
            }
            return result;
        }

        static double Fact(double n)
        {
            return Math.Exp(FactLn(n));
        }

        // http://www.luschny.de/math/factorial/approx/SimpleCases.html
        private static double FactLn(double n)
        {
            double[] c = new[] {1.0 / 24, 3.0 / 80, 18029.0 / 45360, 6272051.0 / 14869008};
            double N = n + 0.5;
            double p = (N * N) / (N + c[0] / (N + c[1] / (N + c[2] / (N + c[3] / N))));
            double logFact = Math.Log(2 * Math.PI) / 2 + N * (Math.Log(p) - 1);
            return logFact;
        }

        static double GammaLn2(double z)
        {
            return z < 1 ?  GammaLn2(z+1) - Math.Log(z) : FactLn(z - 1);
        }

        static double PDF(double x, int k)
        {
            if (x < 1e-10) return 0;
            double k2 = k/2.0;
            double value = (k2-1)*Math.Log(x) - (x/2.0) -(k2 * Math.Log(2) + GammaLn2(k2));
            return Math.Exp(value);
        }

        static double subGv(double a, double x)
        {
            double sum = 0.0, t = 1.0;
            for (int k = 0; k < 100; k++)
            {
                sum += t;
                t *= -(a + k) * x;
                t /= (a + k + 1) * (k + 1);
                //if (Math.Abs(t) < 1e-30) break;
            }

            return Math.Pow(x, a + 1) * sum / (a + 1);
        }

        static double subGu(double a, double x)
        {
            return Math.Exp(GammaLn2(a)) - Math.Pow(x, a) / a;
        }

        static double GammaUpper(double a, double x)
        {
            return subGv(a, x) + subGu(a, x);
        }

        static double GammaUpper2(double a, double x)
        {
            double sum = 0.0, t = 1.0;
            for (int k = 0; k < 50; k++)
            {
                sum += t;
                t *= x;
                t /= (k + 1);
                //if (Math.Abs(t) < 1e-30) break;
            }

            return Math.Exp(-x) * Fact(a + 1) * sum;
        }

        static double GammaUpper3(double s, double z)
        {
            Func<int, double> chis = n => n == 0 ? s : n * (s - n);

            double t = 0.0;
            bool first = true;
            for (int n = 30; n >= 0; n--, first = false) t = (2 * n - 1) + z - s + (first ? 0 : chis(n) / t);
            return Math.Pow(z, s - 1) * Math.Exp(-z) * t;
        }

        static double GammaLower(double s, double z)
        {
            double sum = 0.0, t = 1.0 / s;
            for (int k = 0; k < 50; k++)
            {
                sum += t;
                t *= z;
                t /= (k + 1);
                //if (Math.Abs(t) < 1e-30) break;
            }
            return Math.Pow(z,s) * Math.Exp(-z) * sum;
        }

        static double IntegrPDF(int k, double x)
        {
            return Common.Quadrature.Integration.Adaptive(t => PDF(t, k), 0, x, 10, (e,v,a,b) => Math.Abs(e-v) < 0.05 * v);
        }

        static void Main(string[] args)
        {
            var igamma = new List<double>();
            for (int i = 1; i < 100; i++)
            {
                double x = i / 10.0;
                igamma.Add(GammaLower3(0.5, x));
            }


            //var r0 = new List<double>();
            //for (int i = 1; i <= 7; i++) r0.Add(Common.Quadrature.Integration.Trap(t => gammaSubIntegr(t, s), 0, x, 1 << i));

            //double m4 = 4;

            //while (r0.Count > 1)
            //{
            //    var r1 = new List<double>();
            //    for (int i = 1; i < r0.Count; i++) r1.Add((m4 * r0[i] - r0[i - 1]) / (m4 - 1));
            //    m4 *= 4;
            //    r0 = r1;
            //}

            //var q095 = 3.8415;
            //var q0975 = 5.0239;
            //var q0995 = 7.87944;

            ////check(Quantile1, 0.95, 1, q095, 0.05);
            ////check(Quantile2, 0.95, 1, q095, 0.05);

            ////check(Quantile1, 0.975, 1, q0975, 0.05);
            ////check(Quantile2, 0.975, 1, q0975, 0.05);

            ////check(Quantile1, 0.995, 1, q0995, 0.05);
            ////check(Quantile2, 0.995, 1, q0995, 0.05);

            //var g10 = Math.Exp(GammaLn(10));

            //var g10b = Math.Exp(GammaLn2(10));

            //var p0975 = IntegrPDF(1, q0975);

            // -0,887260933725335
            // 0,729059639748926



            //Console.WriteLine(GammaUpper(1.0, 1.0));
            //Console.WriteLine(GammaUpper(0.5, 0.5));

            //Console.WriteLine();

            //Console.WriteLine(GammaUpper2(1.0, 1.0));
            //Console.WriteLine(GammaUpper2(0.5, 0.5));

            //Console.WriteLine();

            //Console.WriteLine(Math.Exp(GammaLn2(0.5)));
            //Console.WriteLine(Math.Exp(GammaLn2(0.5)) - GammaLower(0.5, 0.5));

            //Console.WriteLine();

            //Console.WriteLine(GammaUpper3(1.0, 1.0));
            //Console.WriteLine(GammaUpper3(0.5, 0.5));

            //Console.WriteLine(Functions.Erf(0.02));
            //Console.WriteLine(Functions.Erf(0.2));
            //Console.WriteLine(Functions.Erf(0.3));
            //Console.WriteLine(Functions.Erf(0.5));
            //Console.WriteLine(Functions.Erf(0.9));
            //Console.WriteLine(Functions.Erf(1.2));
            //Console.WriteLine(Functions.Erf(1.5));
            //Console.WriteLine(Functions.Erf(1.7));
            //Console.WriteLine(Functions.Erf(2.0));
            //Console.WriteLine(Functions.Erf(3.0));
            //Console.WriteLine(Functions.Erf(3.5));

            //Console.WriteLine(Functions.GammaLower05(100));
            //Console.WriteLine(Functions.GammaLower05(10));
            //Console.WriteLine(Functions.GammaLower05(5));
            //Console.WriteLine(Functions.GammaLower05(3));
            //Console.WriteLine(Functions.GammaLower05(2));
            //Console.WriteLine(Functions.GammaLower05(1.5));
            //Console.WriteLine(Functions.GammaLower05(1));
            //Console.WriteLine(Functions.GammaLower05(0.7));
            //Console.WriteLine(Functions.GammaLower05(0.5));
            //Console.WriteLine(Functions.GammaLower05(0.25));
            //Console.WriteLine(Functions.GammaLower05(0.05));

            //Console.WriteLine(GammaLower4(0.5, 0.7));
            //Console.WriteLine(GammaLower4(0.5, 0.5));
            //Console.WriteLine(GammaLower4(0.5, 0.25));
            //Console.WriteLine(GammaLower4(0.5, 0.05));

            Console.WriteLine(Common.ChiSquare.Quantile(1, 0.95));
            Console.WriteLine(Common.ChiSquare.Quantile(1, 0.99));
            Console.WriteLine(Common.ChiSquare.Quantile(1, 0.995));
            Console.WriteLine(Common.ChiSquare.Quantile(1, 0.999));

            Console.WriteLine(Common.ChiSquare.Quantile(7, 0.95));

            Console.ReadLine();
        }

        private static double GammaLower3(double s, double x)
        {
            bool stop(double e, double v, double a, double b) => (b - a) < 1e-10 || Math.Abs(0.001 * v) >= Math.Abs(e);
            return Common.Quadrature.Integration.Adaptive(t => gammaSubIntegr(t, s), 0, x, 2, stop);
        }

        private static double gammaSubIntegr(double t, double s)
        {
            return t < 1e-30 ? 0.0 : Math.Pow(t, s - 1) * Math.Exp(-t);
        }
    }
}
