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
    # print(transposed)
    sum_scalar_key = sum_columns(sum_scalar_key)
    sum_scalar_key1 = sum_scalar_key[1]
    # print(sum_scalar_key1)
    transposed1 = [x * sum_scalar_key1 for x in transposed]
    return transposed1


def generate_mac_share_matrix(mat, m, n):
    matrix = []
    k = 20000 * m
    for i in range(m - 1):
        row = [str(i + 1)] + [str(random.randint(1, k)) for _ in range(n)]
        matrix.append(row)
    row = [str(m)]
    for i in range(n):
        sum = 0
        for j in range(m - 1):
            sum += int(matrix[j][i + 1])
            if mat[i + 1] - sum <= 0:
                sub = math.ceil((sum - mat[i + 1] + 0.05 * k) / (m - 1))
                for l in range(m - 1):
                    matrix[l][i + 1] = str(int(matrix[l][i + 1]) - sub)
                sum1 = sub * (m - 1)
                sum -= sum1
        row += [str(mat[i + 1] - sum)]
    matrix.append(row)
    return matrix


def save_matrix_to_file(matrix, filename):
    with open(filename, 'w') as file:
        for row in matrix:
            line = ' '.join(row) + '\n'
            file.write(line)


def sum_columns(matrix):
    return [sum(map(int, col)) for col in matrix]


filename1 = 'a.txt'
filename2 = 'scalar_key.txt'
matrix1 = read_file(filename1)
matrix2 = read_file(filename2)
m = len(matrix1)  # number of rows
n = len(matrix1[0]) - 1  # number of columns
transposed = mac_matrix(matrix1, matrix2)
# print(transposed)
mac_share_matrix = generate_mac_share_matrix(transposed, m, n)
filename3 = 'scalar_mac_a.txt'
save_matrix_to_file(mac_share_matrix, filename3)
print("Output successfully")
