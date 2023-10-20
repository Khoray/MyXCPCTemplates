inline ll pollard_rho(ll x) {
    auto f = [&](ll x, ll c, ll n) {
        return ((__int128) x * x + c) % n;
    };
    ll s = 0, t = 0, c = 1ll * rand() % (x - 1) + 1;
    ll stp = 0, goal = 1;
    ll val = 1;
    for(goal = 1;; goal <<= 1, s = t, val = 1) {
        for(stp = 1; stp <= goal; ++stp) {
            t = f(t, c, x);
            val = (__int128)val * abs(t - s) % x;
            if((stp % 127) == 0) {
                ll d = __gcd(val, x);
                if(d > 1)
                    return d;
            }
        }
        ll d = __gcd(val, x);
        if(d > 1)
            return d;
    }
}

inline void get_factor_a(ll x, vector<ll> &fac) {
    if(x < 2) return;
    if(is_prime(x)) {
        fac.push_back(x);
        return;
    }
    ll p = x;
    while(p >= x) p = pollard_rho(x);
    while((x % p) == 0) x /= p;
    get_factor_a(x, fac), get_factor_a(p, fac);
}