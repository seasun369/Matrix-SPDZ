import pandas as pd
import random
import numpy as np

from openpyxl.utils import rows_from_range

# 读取 CSV 文件
file_path = 'Student_performance_data _.csv'
df = pd.read_csv(file_path)

# 选择需要的列
columns_to_process = [
    'Age', 'Gender', 'Ethnicity', 'ParentalEducation', 'StudyTimeWeekly',
    'Absences', 'Tutoring', 'ParentalSupport', 'Extracurricular',
    'Sports', 'Music', 'Volunteering', 'GradeClass'
]

k = 2305843009213693951

# 处理样本
X = df[columns_to_process].values
X_mean = np.mean(X, axis=0)
X_std = np.std(X, axis=0)
X = (X - X_mean) / X_std
X = X * pow(2,15)
X = np.c_[np.ones(X.shape[0]), X]
X[:, 0] *= pow(2, 15)
X = X.astype(int)
X = X % k

# 训练集、测试集划分(80% 训练集，20% 测试集)
train_ratio = 0.8
split_idx = int(len(X) * train_ratio)
X = X[:split_idx]

matrix = []
for i in range(0, split_idx):
    matrix = matrix + [int(X[i][0])]
    matrix = matrix + [int(X[i][j + 1]) for j,col in enumerate(columns_to_process)]

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
matrix= generate_share_matrix(matrix, 2, split_idx * (len(columns_to_process) + 1))
filename = 'x.txt'
save_matrix_to_file(matrix, filename)
