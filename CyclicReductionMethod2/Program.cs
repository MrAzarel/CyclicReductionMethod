class Program
{
    static void Main()
    {
        int N = 100;  // Количество узлов сетки
        double[] x, y;

        SolvePeriodicBVP(N, out x, out y);

        // Вывод результатов
        Console.WriteLine("x\ty(x)");
        for (int i = 0; i < N + 1; i++)
        {
            Console.WriteLine($"{x[i]}\t{y[i]}");
        }
    }

    static void SolvePeriodicBVP(int N, out double[] x, out double[] y)
    {
        double[] a, b, c, d;

        PeriodicBoundaryConditions(N, out a, out b, out c);

        double h = 1.0 / N;
        x = new double[N + 1];
        y = new double[N + 1];
        d = new double[N + 1];

        for (int i = 0; i <= N; i++)
        {
            x[i] = i * h;
            d[i] = h * h * F(x[i]);
        }

        // Решение системы методом циклической прогонки
        y = CyclicTridiagonalSolver(a, b, c, d);
    }

    static void PeriodicBoundaryConditions(int N, out double[] a, out double[] b, out double[] c)
    {
        a = new double[N];
        b = new double[N + 1];
        c = new double[N];

        for (int i = 0; i < N; i++)
        {
            a[i] = 1;
            b[i] = -2;
            c[i] = 1;
        }

        // Имитация периодических граничных условий
        a[0] = 1;
        c[N - 1] = 1;
    }

    static double F(double x)
    {
        // Заданная функция
        return Math.Sin(2 * Math.PI * x);
    }

    static double[] CyclicTridiagonalSolver(double[] a, double[] b, double[] c, double[] d)
    {
        int N = b.Length;
        double[] alpha = new double[N];
        double[] beta = new double[N];
        double[] y = new double[N];

        // Прямая прогонка
        alpha[0] = c[0] / b[0];
        beta[0] = d[0] / b[0];

        for (int i = 1; i < N - 1; i++)
        {
            alpha[i] = c[i] / (b[i] - a[i - 1] * alpha[i - 1]);
            beta[i] = (a[i - 1] * beta[i - 1] + d[i]) / (b[i] - a[i - 1] * alpha[i - 1]);
        }

        // Обратная прогонка
        y[N - 1] = beta[N - 1];

        for (int i = N - 2; i >= 0; i--)
        {
            y[i] = alpha[i] * y[(i + 1) % N] + beta[i];
        }

        return y;
    }
}
