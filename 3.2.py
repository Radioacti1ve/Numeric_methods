import matplotlib.pyplot as plt

def compute_spline_coeffs(x, f):
    n = len(x)
    h = [x[i+1] - x[i] for i in range(n-1)]

    # Правая часть системы
    alpha = [0.0] * n
    for i in range(2, n):
        alpha[i] = 3 * ((f[i] - f[i-1]) / h[i-1] - (f[i-1] - f[i-2]) / h[i-2])

    l = [0.0] * n ## изменённая главная диагональ
    mu = [0.0] * n ## вспомогательные коэффициенты для обратного хода
    z = [0.0] * n ## модифицированная правая часть
    l[1] = 1.0
    mu[1] = 0.0
    z[1] = 0.0

    # Прямой ход прогонки
    for i in range(2, n):
        l[i] = 2 * (h[i-2] + h[i-1]) - h[i-2] * mu[i-1]
        mu[i] = h[i-1] / l[i]
        z[i] = (alpha[i] - h[i-2] * z[i-1]) / l[i]

    l[n-1] = 1.0
    z[n-1] = 0.0
    a = f[:-1]
    c = [0.0] * n
    b = [0.0] * (n-1)
    d = [0.0] * (n-1)

    for j in range(n-2, 0, -1):
        c[j] = z[j] - mu[j] * c[j+1]

    c[0] = 0.0
    c[n-1] = 0.0

    for j in range(1, n):
        b[j-1] = (f[j] - f[j-1]) / h[j-1] - (h[j-1] * (c[j] + 2*c[j-1])) / 3
        d[j-1] = (c[j] - c[j-1]) / (3*h[j-1])

    return a, b, c, d

def evaluate_spline(x, a, b, c, d, X_):
    for i in range(len(x)-1):
        if x[i] <= X_ <= x[i+1]:
            return a[i] + b[i]*(X_ - x[i]) + c[i]*(X_ - x[i])**2 + d[i]*(X_ - x[i])**3
    print("Точка X* вне сплайна")    
    return None

def plot_spline(x, f, a, b, c, d, X_=None, f_star=None):
    xs = []
    ys = []
    step = 0.01
    cur = x[0]
    while cur <= x[-1]:
        xs.append(cur)
        ys.append(evaluate_spline(x, a, b, c, d, cur))
        cur += step

    plt.figure(figsize=(8,6))
    plt.plot(x, f, 'o', label='Исходные узлы')
    plt.plot(xs, ys, '-', label='Кубический сплайн')
    if X_ is not None and f_star is not None:
        plt.plot(X_, f_star, 's', markersize=8, label=f'f({X_:.2f}) = {f_star:.5f}')
    plt.grid(True)
    plt.title('Кубический сплайн интерполяции')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend()
    plt.show()

x = [1.0, 1.9, 2.8, 3.7, 4.6]
f = [2.8069, 1.8279, 1.6091, 1.5713, 1.5663]
X_ = 2.66666667  # X*

a, b, c, d = compute_spline_coeffs(x, f)
F_ = evaluate_spline(x, a, b, c, d, X_)
print(f"f(X*) = {F_:.6f}")

plot_spline(x, f, a, b, c, d, X_, F_)