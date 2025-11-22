using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Threading.Tasks;

class IntegralsFunctions
{
    public static double FunctionA(double x) // e^x / (e^x + 1) dx
    {   
        double expX = Math.Exp(x);
        return expX / (expX + 1);
    }

    public static double FunctionB(double x) // x^2 * cos^2(x) dx
    {
        double cosX = Math.Cos(x);
        return x * x * cosX * cosX;
    }

    public static double FunctionC(double x) // sin(x) / (cos^2(x) - 2*cos(x) + 5) dx
    {
        double cosX = Math.Cos(x);
        double denominator = cosX * cosX - 2 * cosX + 5;
        return Math.Sin(x) / denominator;
    }

    public static double FunctionD(double x) // (sqrt(x) - 1) * (6sqrt(x) + 1) / 3(sqrt(x²)) dx
    {
        if (x <= 0) 
            throw new ArgumentException("x должен быть положительным");

        double sqrtX = Math.Sqrt(x);          // x^(1/2)
        double sixthRootX = Math.Pow(x, 1.0 / 6.0); // x^(1/6)
        double cubeRootX2 = Math.Pow(x, 2.0 / 3.0); // x^(2/3)

        return (sqrtX - 1) * (sixthRootX + 1) / cubeRootX2;
    }

    public static double FunctionE(double x) // (x^3 - 2) / (x^2 - 5x + 6) dx
    {
        double denominator = x * x - 5 * x + 6;
        if (Math.Abs(denominator) < 1e-10) 
            throw new ArgumentException("знаменатель ноль =(");

        return (x * x * x - 2) / denominator;
    }
}

static class IntegralMethods
{
    public static double RightRect(Func<double, double> f, double a, double b, long n)
    {
        double h = (b - a) / n;
        double sum = 0;
        for (long i = 1; i <= n; i++)
            sum += f(a + i * h);
        return sum * h;
    }

    public static double Trapezoid(Func<double, double> f, double a, double b, long n)
    {
        double h = (b - a) / n;
        double sum = (f(a) + f(b)) / 2.0;
        for (long i = 1; i < n; i++)
            sum += f(a + i * h);
        return sum * h;
    }

    public static double Simpson(Func<double, double> f, double a, double b, long n)
    {
        if (n % 2 == 1) n++;
        double h = (b - a) / n;
        double sum = f(a) + f(b);
        for (long i = 1; i < n; i++)
            sum += f(a + i * h) * (i % 2 == 1 ? 4 : 2);
        return sum * h / 3.0;
    }
}

static class ParallelIntegralMethods
{
    public static double RightRectPar(Func<double, double> f, double a, double b, long n, int threads = 4)
    {
        double h = (b - a) / n;
        double sum = ParallelSum(i => f(a + i * h), 1, n + 1, threads);
        return sum * h;
    }

    public static double TrapezoidPar(Func<double, double> f, double a, double b, long n, int threads = 4)
    {
        double h = (b - a) / n;
        double sum = (f(a) + f(b)) / 2.0 + ParallelSum(i => f(a + i * h), 1, n, threads);
        return sum * h;
    }

    public static double SimpsonPar(Func<double, double> f, double a, double b, long n, int threads = 4)
    {
        if (n % 2 == 1) n++;
        double h = (b - a) / n;
        double sum = f(a) + f(b) + ParallelSum(i => f(a + i * h) * ((i & 1) == 1 ? 4.0 : 2.0), 1, n, threads);
        return sum * h / 3.0;
    }

   public static double ParallelSum(Func<long, double> term, long start, long end, int threads)
    {
        long n = end - start; // OK
        long chunk = (n + threads - 1) / threads; 


        double[] local = new double[threads];

        Parallel.For(0, threads, t =>
        {
            long s = start + t * chunk;
            long e = (t == threads - 1) ? end : s + chunk;

            double sum = 0;
            for (long i = s; i < e; i++)
                sum += term(i);

            local[t] = sum;
        });

        double total = 0;
        for (int i = 0; i < threads; i++)
            total += local[i];

        return total;
    }
}

static class IntegrationUtils
{
    public static long FindOptimalN(Func<double, double> f, double a, double b, double eps, Func<Func<double, double>, double, double, long, double> integrator)
    {
        long N = 100;
        int maxIterations = 20;
        
        for (int i = 0; i < maxIterations; i++)
        {
            double result1 = integrator(f, a, b, N);
            double result2 = integrator(f, a, b, N * 2);
            
            if (Math.Abs(result1 - result2) < eps)
            {
                Console.WriteLine($"Найдено N = {N * 2}, погрешность = {Math.Abs(result1 - result2):E10}");
                return N * 2;
            }
            
            N *= 2;
            if (N > 10000000) break;
        }
        
        Console.WriteLine($"Достигнуто максимальное N = {N}");
        return N;
    }
}

class Program
{
    static void Main(string[] args)
    {
        const int THREADS = 4;
        double[] accuracies = { 1e-4, 1e-6, 1e-8, 1e-10 };
        const int REPEATS = 100;

        var functions = new List<Func<double, double>> { 
            IntegralsFunctions.FunctionA, 
            IntegralsFunctions.FunctionB, 
            IntegralsFunctions.FunctionC, 
            IntegralsFunctions.FunctionD, 
            IntegralsFunctions.FunctionE 
        };

        double[][] intervals = new double[][]
        {
            new double[] {1.0, 5.0}, // funcA
            new double[] {0.0, 5.0}, // funcB
            new double[] {1.0, 5.0}, // funcC
            new double[] {1.0, 3.0}, // funcD
            new double[] {4.0, 6.0}  // funcE
        };

        string[] functionNames = { "A", "B", "C", "D", "E" };
        string[] methodNames = { "Прямоугольники", "Трапеции", "Симпсон" };

        var integratorsSeq = new Func<Func<double, double>, double, double, long, double>[] 
        {
            IntegralMethods.RightRect,
            IntegralMethods.Trapezoid,
            IntegralMethods.Simpson
        };
        
        var integratorsPar = new Func<Func<double, double>, double, double, long, int, double>[] 
        {
            ParallelIntegralMethods.RightRectPar,
            ParallelIntegralMethods.TrapezoidPar,
            ParallelIntegralMethods.SimpsonPar
        };

        // Для каждого метода
        for (int methodIdx = 0; methodIdx < methodNames.Length; methodIdx++)
        {
            Console.WriteLine($"\n=== МЕТОД: {methodNames[methodIdx]} ===");
            
            // Для каждой функции
            for (int funcIdx = 0; funcIdx < functions.Count; funcIdx++)
            {
                var f = functions[funcIdx];
                double a = intervals[funcIdx][0];
                double b = intervals[funcIdx][1];
                
                Console.WriteLine($"\n--- Функция {functionNames[funcIdx]} на [{a:F1}, {b:F1}] ---");
                
                // Для каждой точности
                foreach (double eps in accuracies)
                {
                    Console.WriteLine($"\nТочность: {eps:E0}");
                    
                    // оптимальное N
                    long optimalN = IntegrationUtils.FindOptimalN(f, a, b, eps, 
                        (func, start, end, n) => integratorsSeq[methodIdx](func, start, end, n));
                    
                    Console.WriteLine($"Оптимальное N: {optimalN:N0}");
            
                    var seqTimeNs = MeasureTime(() =>
                    {
                        for (int i = 0; i < REPEATS; i++)
                            integratorsSeq[methodIdx](f, a, b, optimalN);
                    });

                    double seqResult = integratorsSeq[methodIdx](f, a, b, optimalN);
                    double seqTimeMs = seqTimeNs / (double)REPEATS / 1_000_000.0;
                    Console.WriteLine($"Последовательный: {seqTimeMs:F3} мс, результат = {seqResult:F8}");

                    if (methodIdx == 2 && optimalN % 2 == 1)
                        optimalN++;

                    
                    var parTimeNs = MeasureTime(() => 
                    {
                        for (int i = 0; i < REPEATS; i++)
                            integratorsPar[methodIdx](f, a, b, optimalN, THREADS);
                    });
                    
                    double parResult = integratorsPar[methodIdx](f, a, b, optimalN, THREADS);
                    double parTimeMs = parTimeNs / (double)REPEATS / 1_000_000.0;
                    double acceleration = seqTimeNs / (double)parTimeNs;
                    double efficiency = acceleration / THREADS;
                    
                    Console.WriteLine($"Параллельный ({THREADS} потоков): {parTimeMs:F3} мс, " +
                                      $"ускорение = {acceleration:F2}, " +
                                      $"эффективность = {efficiency:F2}, " +
                                      $"результат = {parResult:F8}");
                }
            }
        }
        
        Console.WriteLine("\nЭксперимент завершён. Нажмите любую клавишу...");
        Console.ReadKey();
    }
    

    
    static long MeasureTime(Action action)
    {
        var sw = Stopwatch.StartNew();
        action();
        sw.Stop();
        return sw.ElapsedTicks * (1_000_000_000L / Stopwatch.Frequency);
    }

}