# Функция умножения двух матриц A и B
def multiply_matrix(A, B):
    # Проверка, совместимы ли размеры для умножения
    if len(A[0]) != len(B):
        raise ValueError("Невозможно умножить матрицы такого размера")
    
    # Размер результата
    rows_A = len(A)
    cols_B = len(B[0])
    # Инициализация результирующей матрицы нулями
    C = [[0 for _ in range(cols_B)] for _ in range(rows_A)]
    
    # Умножение матриц по определению
    for i in range(rows_A):
        for j in range(cols_B):
            for k in range(len(B)):
                C[i][j] += A[i][k] * B[k][j]
    return C

# Функция транспонирования матрицы A
def transpose_matrix(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]

# Функция для создания единичной матрицы размера n x n
def identity_matrix(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]

# Функция для выполнения QR-разложения методом отражений Хаусхолдера
def qr_decomposition(A):
    n = len(A)
    Q = identity_matrix(n)  # Начинаем с единичной матрицы Q
    A_k = [row.copy() for row in A]  # Работаем с копией матрицы A
    
    # Итерация по каждому столбцу, кроме последнего
    for i in range(n - 1):
        v = [0] * n 
        sum_sq = 0
        # Считаем норму вектора столбца начиная с элемента i
        for j in range(i, n):
            sum_sq += A_k[j][i] ** 2
        norm = sum_sq ** 0.5  # Евклидова норма
        
        # Выбор знака для повышения устойчивости
        sign = 1 if A_k[i][i] >= 0 else -1
        v[i] = A_k[i][i] + sign * norm  # Формируем вектор отражения
        
        # Заполняем оставшиеся элементы вектора
        for j in range(i + 1, n):
            v[j] = A_k[j][i]
        
        # Строим матрицу v * v^T (матрицу ранга 1)
        v_matrix = [[v[i] * v[j] for j in range(n)] for i in range(n)]
        v_norm = sum(v[i]**2 for i in range(n))  # Квадрат нормы вектора v
        
        # Строим матрицу Хаусхолдера H
        H = identity_matrix(n)
        for x in range(n):
            for y in range(n):
                H[x][y] -= 2 * v_matrix[x][y] / v_norm
        
        Q = multiply_matrix(Q, H)
        A_k = multiply_matrix(H, A_k)
    
    # Возвращаем ортогональную матрицу Q и верхнетреугольную A_k (R)
    return Q, A_k

# Функция для вычисления собственных значений через QR-алгоритм
def calculate_eigenvalues(A, eps=1e-10, max_iter=1000):
    A_k = [row.copy() for row in A]  # Копируем исходную матрицу для работы
    
    # Основной цикл QR-алгоритма
    for _ in range(max_iter):
        converged = True  # Флаг сходимости
        # Проверка, что все элементы ниже главной диагонали малы
        for i in range(1, len(A)):
            for j in range(i):
                if abs(A_k[i][j]) > eps:
                    converged = False
                    break
            if not converged:
                break
        
        if converged:
            break  # Если сходимость достигнута, выходим
        
        # QR-разложение текущей матрицы
        Q, R = qr_decomposition(A_k)
        # Пересобираем матрицу как R * Q (а не Q * R!)
        A_k = multiply_matrix(R, Q)
    
    # Собственные значения — диагональные элементы A_k после сходимости
    eigenvalues = [A_k[i][i] for i in range(len(A_k))]
    return eigenvalues

# === Пример использования ===

# Исходная матрица
matrix = [
    [1, 7, -1],
    [-2, 2, -2],
    [9, -7, 3]
]

# Вызываем алгоритм для вычисления собственных значений
eigenvalues = calculate_eigenvalues(matrix)

# Выводим результат
print(f'Полученные собственные значения: {eigenvalues}')