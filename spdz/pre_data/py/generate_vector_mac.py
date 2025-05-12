import random
import math


def read_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    matrix = [line.strip().split() for line in lines]
    return matrix


def mac_matrix(matrix1, matrix2):
    transposed = []
    sum_scalar_key = []
    for i in range(len(matrix1[0])):
        col = [row[i] for row in matrix1]
        transposed.append(col)
    for i in range(len(matrix2[0])):
        col = [row[i] for row in matrix2]
        sum_scalar_key.append(col)
    transposed = sum_columns(transposed)
    transposed = transposed[1:]
    print('sum of matrix a:', transposed)
    sum_scalar_key = sum_columns(sum_scalar_key)
    sum_scalar_key = sum_scalar_key[1:]
    print('sum of vector key:', sum_scalar_key)
    result = matrix_vector_mult(transposed, sum_scalar_key)
    print('sum of vector mac:', result)
    return result


def generate_mac_share_matrix(mat, m, n):
    matrix = []
    k = 20000 * m * m
    for i in range(m - 1):
        row = [str(i + 1)] + [str(random.randint(1, k)) for _ in range(n)]
        matrix.append(row)
    row = [str(m)]
    for i in range(n):
        sum = 0
        for j in range(m - 1):
            sum += int(matrix[j][i + 1])
            if mat[i] - sum <= 0:
                sub = math.ceil((sum - mat[i] + 0.05 * k) / (m - 1))
                for l in range(m - 1):
                    matrix[l][i + 1] = str(int(matrix[l][i]) - sub)
                sum1 = sub * (m - 1)
                sum -= sum1
        row += [str(mat[i] - sum)]
    matrix.append(row)
    return matrix


def matrix_vector_mult(matrix, vector):
    m = len(vector)
    result = [0] * m
    for i in range(m):
        for j in range(m):
            result[i] += matrix[i * m + j] * vector[j]
    return result


def save_matrix_to_file(matrix, filename):
    with open(filename, 'w') as file:
        for row in matrix:
            line = ' '.join(row) + '\n'
            file.write(line)


def sum_columns(matrix):
    return [sum(map(int, col)) for col in matrix]


filename1 = 'a.txt'
filename2 = 'vector_key.txt'
matrix1 = read_file(filename1)
matrix2 = read_file(filename2)
m = len(matrix1)  # nP
n = len(matrix2[0]) - 1  # number of key_rows(i.e. m)
sum_vector_mac = mac_matrix(matrix1, matrix2)
mac_share_matrix = generate_mac_share_matrix(sum_vector_mac, m, n)
filename3 = 'vector_mac_a.txt'
save_matrix_to_file(mac_share_matrix, filename3)
print("Output successfully")
