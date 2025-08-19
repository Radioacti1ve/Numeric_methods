import matplotlib.pyplot as plt

def gauss_method(A, B):
    A = [row[:] for row in A]
    B = B[:]
    n = len(B)

    ## Нормализуем построчно
    for i in range(n):
        max = A[i][i]
        for j in range(i, n):
            A[i][j] /= max
        B[i] /= max

        ## приводим к верхнетреугольной матрице A, обнуляя столбцы
        for k in range(i + 1, n):
            factor = A[k][i]
            for j in range(i, n):
                A[k][j] -= factor * A[i][j]
            B[k] -= factor * B[i]

    X = [0] * n
    for i in range(n - 1, -1, -1):
        X[i] = B[i]
        for j in range(i + 1, n):
            X[i] -= A[i][j] * X[j]
    return X

def get_linear_coeffs(x, y):
    n = len(x)
    sx = sum(x)
    sx2 = sum(xi ** 2 for xi in x)
    sy = sum(y)
    sxy = sum(x[i] * y[i] for i in range(n))

    A = [[n, sx],
         [sx, sx2]]
    B = [sy, sxy]

    return gauss_method(A, B)

def get_square_coeffs(x, y):
    n = len(x)
    sx = sum(x)
    sx2 = sum(xi ** 2 for xi in x)
    sx3 = sum(xi ** 3 for xi in x)
    sx4 = sum(xi ** 4 for xi in x)
    sy = sum(y)
    sxy = sum(x[i] * y[i] for i in range(n))
    sx2y = sum((x[i] ** 2) * y[i] for i in range(n))

    A = [[n, sx, sx2],
         [sx, sx2, sx3],
         [sx2, sx3, sx4]]
    B = [sy, sxy, sx2y]

    return gauss_method(A, B)
    
x = [1.0, 1.9, 2.8, 3.7, 4.6, 5.5]
y = [3.4142, 2.9818, 3.3095, 3.8184, 4.3599, 4.8318]

a0, a1 = get_linear_coeffs(x, y)
b0, b1, b2 = get_square_coeffs(x, y)

def f_lin(x_): return a0 + a1 * x_
def f_sq(x_):  return b0 + b1 * x_ + b2 * x_ * x_

print(f"Многочлен 1-й степени: f(x) = {a0:.5f} + {a1:.5f} * x")
print(f"Многочлен 2-й степени: f(x) = {b0:.5f} + {b1:.5f} * x + {b2:.5f} * x^2")

error1 = sum((f_lin(x[i]) - y[i]) ** 2 for i in range(len(x)))
error2 = sum((f_sq(x[i]) - y[i]) ** 2 for i in range(len(x)))
print(f"Сумма квадратов ошибок для многочлена 1-й степени: {error1:.5f}")
print(f"Сумма квадратов ошибок для многочлена 2-й степени: {error2:.5f}")

# построим кривые вблизи диапазона узлов
lo, hi = min(x) - 0.2, max(x) + 0.2
x_val = [lo + i * (hi - lo) / 300 for i in range(301)]
y_val1 = [f_lin(xx) for xx in x_val]
y_val2 = [f_sq(xx) for xx in x_val]

import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.scatter(x, y, color='black', label='Табличные данные', zorder=5)
plt.plot(x_val, y_val1, label='Многочлен 1-й степени', linewidth=2)
plt.plot(x_val, y_val2, label='Многочлен 2-й степени', linestyle='--', linewidth=2)
plt.xlabel('x'); plt.ylabel('y'); plt.title('Аппроксимация МНК')
plt.grid(True); plt.legend(); plt.show()