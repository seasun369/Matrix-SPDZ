import pandas as pd
import random
import numpy as np

from openpyxl.utils import rows_from_range

# 读取 CSV 文件
file_path = 'Student_performance_data _.csv'
df = pd.read_csv(file_path)

# 选择需要的列
columns_to_process = [
    'GPA'
]

k = 2305843009213693951

# 处理样本
y = df[columns_to_process].values
y = y / np.max(y)
y = y * pow(2,15)
y = y.astype(int)
y = y % k

train_ratio = 0.8
split_idx = int(len(y) * train_ratio)
y = y[:split_idx]

matrix = []
for i in range(0, split_idx):
    matrix = matrix + [int(y[i][j]) for j,col in enumerate(columns_to_process)]


def generate_share_matrix(mat, m, n):
    matrix = []
    for i in range(m - 1):
        row = [str(i + 1)] + [str(random.randint(1, k)) for _ in range(n)]
        matrix.append(row)
    row = [str(m)]
    for i in range(n):
        sum = 0
        for j in range(m - 1):
            sum = (sum + int(matrix[j][i + 1])) % k
        row += [str((mat[i] - sum) % k)]
    matrix.append(row)
    return matrix


def save_matrix_to_file(matrix, filename):
    with open(filename, 'w') as file:
        for row in matrix:
            line = ' '.join(row) + '\n'
            file.write(line)

print(matrix)
matrix= generate_share_matrix(matrix, 2, split_idx * len(columns_to_process))
filename = 'y.txt'
save_matrix_to_file(matrix, filename)
