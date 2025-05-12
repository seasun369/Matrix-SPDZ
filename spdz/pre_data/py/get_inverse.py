def extended_gcd(a, b):
    if b == 0:
        return a, 1, 0
    gcd, x1, y1 = extended_gcd(b, a % b)
    x = y1
    y = x1 - (a // b) * y1
    return gcd, x, y

def mod_inverse(a, m):
    gcd, x, _ = extended_gcd(a, m)
    if gcd != 1:
        raise Exception(f"{a} 与 {m} 不互素, 无逆元")
    return x % m

def main():
    mod = 2305843009213693951
    try:
        a = int(input("请输入元素 a: "))
        inv_a = mod_inverse(a, mod)
        print(f"元素 {a} 在模 {mod} 的群中的逆元为: {inv_a}")
    except Exception as e:
        print(f"出错了: {e}")

if __name__ == "__main__":
    main()
    print((822005582270216360 + 1491160415918477591)%2305843009213693951)
    print((1767726482876860042 + 545723813411833909) % 2305843009213693951)