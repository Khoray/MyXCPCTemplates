bool is_prime(ll n) {
    if(n < 2 || n % 6 % 4 != 1) return (n | 1) == 3;
    int A[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022},
              s = __builtin_ctzll(n - 1), d = n >> s;
    for(int a : A) {  // ^ count t ra i l in g zeroes
        int p = ksm(a % n, d, n), i = s;
        while(p != 1 && p != n - 1 && a % n && i--)
            p = (__int128) p * p % n;
        if(p != n - 1 && i != s) return 0;
    }
    return 1;
}