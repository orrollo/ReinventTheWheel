using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Common.Quadrature
{
    public static class Integration
    {
        public static double Trap(Func<double, double> fn, double a, double b, int n)
        {
            double h = (b - a) / n, result = fn(a) + fn(b), x = a + h;
            for (int i = 1; i < n; i++, x += h) result += 2 * fn(x);
            return h * result / 2.0;
        }

        public static double Simpson(Func<double, double> fn, double a, double b, int n)
        {
            double h = (b - a) / n, result = fn(a) + fn(b), x = a, h2 = h / 2;
            for (int i = 0; i < n; i++, x += h)
            {
                result += 4 * fn(x + h2);
                if (i > 0) result += 2 * fn(x);
            }
            return result * h / 6;
        }

        public static double Adaptive(Func<double, double> fn, double a, double b, int n, Func<double,double,double,double,bool> isAllowedFn)
        {
            double v1 = Simpson(fn, a, b, n), v2 = Simpson(fn, a, b, 2 * n), err = Math.Abs(v1 - v2) / 7;
            if (isAllowedFn(err,v2,a,b)) return v2;
            double c = (a + b) / 2;
            return Adaptive(fn, a, c, n, isAllowedFn) + Adaptive(fn, c, b, n, isAllowedFn);
        }

        public static double Adaptive(Func<double, double> fn, double a, double b, double error = 1e-3, int n = 10)
        {
            double v1 = Simpson(fn, a, b, n), v2 = Simpson(fn, a, b, 2 * n), err = Math.Abs(v1 - v2) / 7;
            if (err <= error || Math.Abs(b-a) < 1e-10) return v2;
            double c = (a + b) / 2;
            return Adaptive(fn, a, c, error / 2, n) + Adaptive(fn, c, b, error / 2, n);
        }
    }
}
