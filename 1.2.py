def is_tridiagonal(matrix):
    n = len(matrix)
    for i in range(n):
        for j in range(n):
            if abs(i - j) > 1 and matrix[i][j] != 0:
                return False
    return True

def get_diagonal(matrix, offset):
    n = len(matrix)
    if offset == 0:  # главная диагональ
        return [matrix[i][i] for i in range(n)]
    elif offset == -1:  # поддиагональ
        return [0] + [matrix[i][i-1] for i in range(1, n)]
    elif offset == 1:  # наддиагональ
        return [matrix[i][i+1] for i in range(n-1)] + [0]
    else:
        raise ValueError("Неподдерживаемый offset")

def run_through_algorithm(matrix, d):
    if not is_tridiagonal(matrix):
        raise ValueError("матрица не является трехдиагональной")
    
    n = len(matrix)
    a = get_diagonal(matrix, -1) # поддиагональ
    b = get_diagonal(matrix, 0)# главная диагональ
    c = get_diagonal(matrix, 1) # наддиагональ
    d = d.copy()
    
    if any(b_i == 0 for b_i in b):
        raise ValueError("главная диагональ содержит нули")
    
    # Прямой ход
    p = [0] * n
    q = [0] * n
    p[0] = -c[0] / b[0]
    q[0] = d[0] / b[0]
    
    for i in range(1, n):
        denom = b[i] + a[i] * p[i-1]
        if denom == 0:
            raise ValueError("деление на ноль")
        p[i] = -c[i] / denom
        q[i] = (d[i] - a[i] * q[i-1]) / denom
    
    x = [0] * n
    x[-1] = q[-1]
    for i in range(n-2, -1, -1):
        x[i] = p[i] * x[i+1] + q[i]
    
    return x

def matrix_multiply(A, B):
    n = len(A)
    m = len(B[0])
    result = [[0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            result[i][j] = sum(A[i][k] * B[k][j] for k in range(len(B)))
    return result

A = [
        [16, -8, 0, 0, 0],
        [-7, -16, 5, 0, 0],
        [0, 4, 12, 3, 0],
        [0, 0, -4, 12, -7],
        [0, 0, 0, -1, 7]
    ]

B = [0, -123, -68, 104, 20]

solution = run_through_algorithm(A, B)
print("Рещение СЛАУ методом прогонки:", solution)
print("\nпроверка:", matrix_multiply(A, [[x] for x in solution]))

