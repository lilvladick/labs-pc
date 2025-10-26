using System;
using System.Diagnostics;
using System.Threading.Tasks;

namespace sh_parallel;

class Program
{
    static double Sh(double x, double eps, out int terms)
    {
        double sum = 0;
        terms = 0;
        while (true)
        {
            int k = 2 * terms + 1;
            double term = Math.Pow(x, k) / Factorial(k);
            if (Math.Abs(term) < eps) break;
            sum += term;
            terms++;
        }
        return sum;
    }

    static double Factorial(int n)
    {
        double result = 1;
        for (int i = 2; i <= n; i++) result *= i;
        return result;
    }

    static double ShParallelNoLock(double x, double eps, int threads)
    {
        double sum = 0;
        int termIndex = 0;
        bool done = false;

        Parallel.For(0, threads, new ParallelOptions { MaxDegreeOfParallelism = threads }, () => 0.0,
            (_, state, localSum) =>
            {
                while (true)
                {
                    int k = Interlocked.Increment(ref termIndex) - 1;
                    k = 2 * k + 1;
                    double term = Math.Pow(x, k) / Factorial(k);

                    if (Math.Abs(term) < eps)
                    {
                        Interlocked.CompareExchange(ref done, true, false);
                        break;
                    }

                    if (done) break;

                    localSum += term;
                }
                return localSum;
            },
            localSum =>
            {
                double currentSum, newSum;
                do
                {
                    currentSum = sum;
                    newSum = currentSum + localSum;
                } while (Interlocked.CompareExchange(ref sum, newSum, currentSum) != currentSum);
            });

        return sum;
    }

    static double ShParallel(double x, double eps, int threads)
    {
        double sum = 0;
        int termIndex = 0;
        object locker = new();
        bool done = false;

        Parallel.For(0, threads, new ParallelOptions { MaxDegreeOfParallelism = threads }, _ =>
        {
            while (true)
            {
                int k;
                lock (locker)
                {
                    if (done) break;
                    k = 2 * termIndex + 1;
                    termIndex++;
                }

                double term = Math.Pow(x, k) / Factorial(k);

                if (Math.Abs(term) < eps)
                {
                    lock (locker)
                    {
                        done = true;
                    }

                    break;
                }

                lock (locker)
                {
                    sum += term;
                }
            }

        });

        return sum;
    }

    static void Main(string[] args)
    {
        const double x = 5;
        double[] epsilons = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10];
        int[] threadCounts = [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24];
        const int repeats = 100;

        foreach (double eps in epsilons)
        {
            Console.WriteLine("Точность: {0}", eps);
            double seqTimeTotal = 0.0;
            int terms = 0;
            for (int r = 0; r < repeats; r++)
            {
                var sw = Stopwatch.StartNew();
                double res = Sh(x, eps, out terms);
                sw.Stop();
                seqTimeTotal += sw.Elapsed.TotalMilliseconds;
            }
            double seqAvg = seqTimeTotal / repeats;
            Console.WriteLine($"Последовательный: {seqAvg:F4} мс");

            foreach (int threads in threadCounts)
            {
                double parTimeTotal = 0.0;
                for (int r = 0; r < repeats; r++)
                {
                    var sw = Stopwatch.StartNew();
                    double res = ShParallelNoLock(x, eps, threads);
                    sw.Stop();
                    parTimeTotal += sw.Elapsed.TotalMilliseconds;
                }

                double parAvg = parTimeTotal / repeats;
                Console.WriteLine($"Параллельный ({threads} потоков): {parAvg:F4} мс");
            }
        }
    }
}