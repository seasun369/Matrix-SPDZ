#include <iostream>

class FiniteFieldElement {
private:
    int value;  // 有限域Fq中的元素值
    int modulo; // 有限域的模数

public:
    // 构造函数，初始化有限域元素
    FiniteFieldElement(int val, int mod) : value(val % mod), modulo(mod) {}

    // 获取有限域元素值
    int getValue() const {
        return value;
    }

    // 有限域加法
    FiniteFieldElement operator+(const FiniteFieldElement& other) const {
        return FiniteFieldElement((value + other.value) % modulo, modulo);
    }

    // 有限域减法
    FiniteFieldElement operator-(const FiniteFieldElement& other) const {
        return FiniteFieldElement((value - other.value + modulo) % modulo, modulo);
    }

    // 有限域乘法
    FiniteFieldElement operator*(const FiniteFieldElement& other) const {
        return FiniteFieldElement((value * other.value) % modulo, modulo);
    }

    // 有限域除法
    FiniteFieldElement operator/(const FiniteFieldElement& other) const {
        // 计算乘法逆元
        int inverse = calculateMultiplicativeInverse(other.value, modulo);
        return FiniteFieldElement((value * inverse) % modulo, modulo);
    }

    // 打印有限域元素
    void print() const {
        std::cout << value;
    }

private:
    // 计算乘法逆元
    int calculateMultiplicativeInverse(int a, int m) const {
        // 使用扩展欧几里得算法计算乘法逆元
        int m0 = m;
        int y = 0, x = 1;

        if (m == 1) return 0;

        while (a > 1) {
            int q = a / m;
            int t = m;

            m = a % m;
            a = t;
            t = y;

            y = x - q * y;
            x = t;
        }

        if (x < 0) x += m0;

        return x;
    }
};

