#include <stdio.h>
#include <mpir.h>

int main() {
    // 初始化MPIR
    //mpir_init();

    // 定义一个多精度整数数组
    mp_limb_t a[3] = {0xFFFFFFFF, 0x12345678, 0xABCDEF01}; // 原始数据
    mp_limb_t result[3]; // 用于存储结果

    // 使用mpn_neg进行取反操作
    mpn_neg(result, a, 3);

    // 打印结果
    printf("Result: ");
    for (int i = 2; i >= 0; i--) { // 反向输出，以符合大端顺序
        printf("%08lx ", result[i]); // 输出每个 limb
    }
    printf("\n");

    // 清理MPIR
    //mpir_clear();

    return 0;
}

