using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using cq = Common.Quadrature;

namespace Integration
{
    class Program
    {
        static double Gamma(double z)
        {
            if (z > 15) return Math.Pow(2, z-1) * Gamma(z / 2) * Gamma((z + 1) / 2) / Math.Sqrt(Math.PI);
            double error = Math.Max(z / Math.Pow(2, z + 2), 1e-5); //Math.Pow(10,-(z+1)/3);
            var integral = cq.Integration.Adaptive(y => y > 0 ? y * Math.Pow(-Math.Log(y), z) : 0, 0, 1, error);
            return Math.Pow(2, z + 1) * integral / z;
        }

        static void Main(string[] args)
        {
            var check1 = cq.Integration.Adaptive(Math.Sin, 0, Math.PI, 1e-7, 10);
            var error = Math.Abs(2.0 - check1);
            if (error > 1e-7) throw new Exception("unexact solution");

            var check2 = Gamma(10);
            error = Math.Abs(362880 - check2);
            if (error > 1) throw new Exception("unexact solution");

            var check3 = Gamma(15);
            error = Math.Abs(87178291200 - check3);
            if (error > 1) throw new Exception("unexact solution");

            var check4 = Gamma(20);
            var fact19 = 121645100408832000;
            error = Math.Abs(fact19 - check4) / fact19; // relative error
            if (error > 1e-4) throw new Exception("unexact solution");
        }
    }
}
