int bsgs(int a, int b, int p) { //BSGS算法
    unordered_map<int, int> f;
    int m = ceil(sqrt(p));
    b %= p;
    for(int i = 1; i <= m; i++) {
        b = b * a % p;
        f[b] = i;
    }
    int tmp = ksm(a, m, p);
    b = 1;
    for(int i = 1; i <= m; i++) {
        b = b * tmp % p;
        if(f[b]) {
            return (i * m - f[b] + p) % p;
        }
    }
    return -1;
}
int exbsgs(int a, int b, int p) {
    b %= p;
    a %= p;
    if(b == 1 || p == 1) {
        return 0;    //特殊情况，x=0时最小解
    }
    int g = __gcd(a, p), k = 0, na = 1;
    while(g > 1) {
        if(b % g != 0) {
            return -1;    //无法整除则无解
        }
        k++;
        b /= g;
        p /= g;
        na = na * (a / g) % p;
        if(na == b) {
            return k;    //na=b说明前面的a的次数为0，只需要返回k
        }
        g = __gcd(a, p);
    }
    int f = bsgs(a, b * inv(na, p) % p, p);
    if(f == -1) {
        return -1;
    }
    return f + k;
}