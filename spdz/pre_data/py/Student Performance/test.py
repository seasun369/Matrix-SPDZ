import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

inverse_accuracy=70368744177664
ppp=2305843009213693951

# 读取数据
file_path = 'Student_performance_data _.csv'
df = pd.read_csv(file_path)

feature_cols = ['Age', 'Gender', 'Ethnicity', 'ParentalEducation', 'StudyTimeWeekly', 'Absences', 'Tutoring', 'ParentalSupport', 'Extracurricular', 'Sports', 'Music', 'Volunteering', 'GradeClass']
target_col = 'GPA'

X = df[feature_cols].values
y = df[target_col].values

# 数据预处理，标准化 X
X_mean = np.mean(X, axis=0)
X_std = np.std(X, axis=0)
X = (X - X_mean) / X_std
X = X * pow(2, 15)
X = X.astype(int)

# 归一化 y
y = y / np.max(y)
y = y * pow(2, 15)
y = y.astype(int)
# y_mean = np.mean(y, axis=0)
# y_std = np.std(y, axis=0)
# y = (y - y_mean) / y_std

# X 前加一列 1
X = np.c_[np.ones(X.shape[0]), X]
X[:, 0] *= pow(2, 15)

# 训练集、测试集划分(80% 训练集，20% 测试集)
train_ratio = 0.8
split_idx = int(len(X) * train_ratio)
X_train, X_test = X[:split_idx], X[split_idx:]
y_train, y_test = y[:split_idx], y[split_idx:]

# 梯度下降
def compute_loss(X, y, w):
    m = len(y)
    predictions = X @ w
    predictions = predictions * pow(2, -15)
    loss = (1 / (2 * m)) * np.sum((predictions - y) ** 2) * pow(2,-15)  # 1/(2m)*(Xw-y)^2
    return loss

def decode(encoding):
    encoding = encoding * pow(2,-15)
    return encoding

def gradient_descent(X, y, num_epochs):
    m, n = X.shape
    w = np.zeros(n)
    losses = []

    for epoch in range(num_epochs):
        predictions = X @ w
        predictions = predictions * pow(2,-15)
        predictions = [int(x) for x in predictions]
        predictions = predictions - y
        #print("predictions:", predictions)
        predictions = X.T @ predictions
        #print("predictions:", predictions)
        predictions = predictions * pow(2, -15)
        predictions = predictions.astype(int)
        #print("predictions:", predictions)
        gradients = predictions * pow(2,-13)
        #print("gradients:", gradients)
        gradients = gradients.astype(int)
        w -= gradients  # 更新参数

        loss = compute_loss(X, y, w)
        losses.append(loss)

        if epoch % 50 == 0:
            print(f"Epoch {epoch}, Loss: {loss:.4f}")

    return w, losses

# 训练模型
w, losses = gradient_descent(X_train, y_train, num_epochs = 50)
print("w:", w)
print("loss:", losses)
y_pred = X_test @ w
y_pred = y_pred * pow(2,-15)

# 计算 R²
def r2_score(y_true, y_pred):
    ss_total = np.sum((y_true - np.mean(y_true)) ** 2)
    ss_residual = np.sum((y_true - y_pred) ** 2)
    return 1 - (ss_residual / ss_total)

# 衡量模型的预测准确性
r2 = r2_score(y_test, y_pred)
print(f"测试集 R² 准确率: {r2:.4f}")

# 计算MSE
mse = np.mean((y_pred - y_test) ** 2) * pow(2, -15)
print(f"测试集均方误差MSE: {mse:.4f}")

'''
# 绘制散点图（预测结果）
plt.figure(figsize=(8, 6))
y_test = y_test
y_pred = y_pred
plt.scatter(y_test * pow(2, -15), y_pred * pow(2, -15), alpha=0.5, color='b')  # alpha 设置透明度
plt.plot([0, 1], [0, 1], color='r', linestyle='--')  # 添加一条红色虚线表示完美预测

plt.title('Predicted vs Actual values')
plt.xlabel('Actual Values (y_test)')
plt.ylabel('Predicted Values (y_pred)')
plt.grid(True)
plt.show()
'''

# 创建图表Losses
plt.plot(losses, label='Loss')  # 绘制损失值曲线
plt.xlabel('Epoch')  # x轴标签
plt.ylabel('Loss')  # y轴标签
plt.title('Loss Curve')  # 图表标题
plt.legend()  # 显示图例
plt.grid(True)  # 显示网格线
plt.show()  # 展示图表