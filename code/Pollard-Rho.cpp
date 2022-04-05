inline int pollard_rho(int x) {
    auto f = [&](int x, int c, int n) {
        return ((__int128) x * x + c) % n;
    };
    int s = 0, t = 0, c = 1ll * rand() % (x - 1) + 1;
    int stp = 0, goal = 1;
    int val = 1;
    for(goal = 1;; goal <<= 1, s = t, val = 1) {
        for(stp = 1; stp <= goal; ++stp) {
            t = f(t, c, x);
            val = (__int128)val * abs(t - s) % x;
            if((stp % 127) == 0) {
                int d = __gcd(val, x);
                if(d > 1)
                    return d;
            }
        }
        int d = __gcd(val, x);
        if(d > 1)
            return d;
    }
}

inline void get_factor_a(int x, vector<int> &fac) {
    if(x < 2) return;
    if(is_prime(x)) {
        fac.push_back(x);
        return;
    }
    int p = x;
    while(p >= x) p = pollard_rho(x);
    while((x % p) == 0) x /= p;
    get_factor_a(x, fac), get_factor_a(p, fac);
}