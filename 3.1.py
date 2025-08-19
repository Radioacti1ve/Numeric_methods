import math
import matplotlib.pyplot as plt

# функция из условия
def f(x):
    return 1/math.tan(x) + x   # ctg(x) + x

def lagrange_interpolation(x, y, X_):
    polynom_str = 'L(x) ='
    polynom_X_ = 0
    for i in range(len(x)):
        cur_enum, cur_denom = 1, 1
        for j in range(len(x)):
            if i == j:
                continue
            cur_enum *= (X_[0] - x[j])
            cur_denom *= (x[i] - x[j])
        polynom_str += f' + {(y[i] / cur_denom):.5f}*' + ''.join([f'(x-{x[j]:.2f})' for j in range(len(x)) if j != i])
        polynom_X_ += y[i] * cur_enum / cur_denom
    return polynom_str, polynom_X_

def newton_interpolation(x, y, X_):
    n = len(x)
    difference = [y[i] for i in range(n)]
    for i in range(1, n):
        for j in range(n - 1, i - 1, -1):
            difference[j] = float(difference[j] - difference[j - 1]) / float(x[j] - x[j - i])

    polynom_str = 'P(x) = '
    polynom_X_ = 0
    cur_coeffs = 1
    for i in range(n):
        polynom_X_ += cur_coeffs * difference[i]
        if i == 0:
            polynom_str += f'{difference[i]:.5f}'
        else:
            polynom_str += ' + ' + ''.join([f'(x-{x[j]:.2f})' for j in range(i)]) + f'*{difference[i]:.5f}'
        cur_coeffs *= (X_[0] - x[i])
    return polynom_str, polynom_X_

def generate_linspace(start, end, num_points):
    step = (end - start) / (num_points - 1)
    return [start + i * step for i in range(num_points)]

def plot_interpolation(x, y, func, lagrange_func, newton_func, title, x_star):
    x_plot = generate_linspace(min(x)-0.5, max(x)+0.5, 500)
    y_exact = []
    y_lagrange = []
    y_newton = []
    for xi in x_plot:
        y_exact.append(func(xi))
        _, yi_l = lagrange_func(x, y, (xi, 0))
        _, yi_n = newton_func(x, y, (xi, 0))
        y_lagrange.append(yi_l)
        y_newton.append(yi_n)

    plt.figure()
    plt.plot(x_plot, y_exact, label='y = ctg(x)+x', linewidth=2)
    plt.plot(x_plot, y_lagrange, '--', label='Лагранж')
    plt.plot(x_plot, y_newton, ':', label='Ньютон')
    plt.scatter(x, y, color='red', label='Точки')
    plt.scatter(x_star, [func(x_star)], color='green', label='X*', zorder=5)
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    # Вариант а
    x_a = [math.pi/8, 2*math.pi/8, 3*math.pi/8, 4*math.pi/8]
    y_a = [f(xi) for xi in x_a]

    # Вариант б
    x_b = [math.pi/8, math.pi/3, 3*math.pi/8, math.pi/2]
    y_b = [f(xi) for xi in x_b]

    X_ = 3*math.pi/16
    Y_ = f(X_)

    print('-----------------Интерполяционный многочлен Лагранжа----------------')
    print('----Вариант а----')
    lagrange_polynom_a, lagrange_value_a = lagrange_interpolation(x_a, y_a, (X_, 0))
    print(lagrange_polynom_a)
    print(f'y(X*) = {Y_:.5f}')
    print(f'L(X*) = {lagrange_value_a:.5f}')
    print(f'Погрешность = {abs(Y_ - lagrange_value_a)}')

    print('----Вариант б----')
    lagrange_polynom_b, lagrange_value_b = lagrange_interpolation(x_b, y_b, (X_, 0))
    print(lagrange_polynom_b)
    print(f'y(X*) = {Y_:.5f}')
    print(f'L(X*) = {lagrange_value_b:.5f}')
    print(f'Погрешность = {abs(Y_ - lagrange_value_b)}')

    print('\n-----------------Интерполяционный многочлен Ньютона-----------------')
    print('----Вариант а----')
    newton_polynom_a, newton_value_a = newton_interpolation(x_a, y_a, (X_, 0))
    print(newton_polynom_a)
    print(f'y(X*) = {Y_:.5f}')
    print(f'P(X*) = {newton_value_a:.5f}')
    print(f'Погрешность = {abs(Y_ - newton_value_a)}')

    print('----Вариант б----')
    newton_polynom_b, newton_value_b = newton_interpolation(x_b, y_b, (X_, 0))
    print(newton_polynom_b)
    print(f'y(X*) = {Y_:.5f}')
    print(f'P(X*) = {newton_value_b:.5f}')
    print(f'Погрешность = {abs(Y_ - newton_value_b)}')

    plot_interpolation(x_a, y_a, f, lagrange_interpolation, newton_interpolation, 'Вариант а', X_)
    plot_interpolation(x_b, y_b, f, lagrange_interpolation, newton_interpolation, 'Вариант б', X_)
