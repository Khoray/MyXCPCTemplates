#include<bits/stdc++.h>
using namespace std;
#define int long long
using ll = int64_t;
using lll = __int128;
#ifdef _KHORAY
#include<debughr.h>
#else
#define out(args...) 42
#endif

#define dec

int ksm(int a, int b, int MOD_KSM) {
    int ret = 1;
    while(b) {
        if(b & 1) {
            ret = (__int128) ret * a % MOD_KSM;
        }
        a = (__int128) a * a % MOD_KSM;
        b >>= 1;
    }
    return ret;
}


bool is_prime(int n) {
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
int pollard(int n) {
    auto f = [n](int x) {
        return (__int128) x * x % (n + 1);
    };
    int x = 0, y = 0, t = 30, prd = 2, i = 1, q;
    while(t++ % 40 || __gcd(prd, n) == 1) {
        if(x == y) x = ++i, y = f(x);
        if((q = (__int128) prd * (max(x, y) - min(x, y)) % n)) prd = q;
        x = f(x), y = f(f(y));
    }
    return __gcd(prd, n);
}
vector<int> factor(int n) {
    if(n == 1) return {};
    if(is_prime(n)) return {n};
    int x = pollard(n);
    auto l = factor(x), r = factor(n / x);
    l.insert(l.end(), r.begin(), r.end());
    return l;
}
mt19937 mt(time(0));

int pollard_rho(int n, int c) {
    int x = uniform_int_distribution<int>(1, n - 1)(mt), y = x;
    auto f = [&](int v) {
        int t = (__int128) v * v % n + c;
        return t < n ? t : t - n;
    };
    while(1) {
        x = f(x);
        y = f(f(y));
        if(x == y) return n;
        int d = __gcd(abs(x - y), n);
        if(d != 1) return d;
    }
}
int max_factor = 0;

void get_fac(int n, vector<int> &fac, int cc = 15653413) {
    if(n == 4) {
        max_factor = max(max_factor, 2ll);
        return;
    }
    if(is_prime(n)) {
        max_factor = max(max_factor, n);
        return;
    }
    int p = n;
    while(p == n) p = pollard_rho(n, --cc);
    get_fac(p, fac);
    get_fac(n / p, fac);
}
vector<int> go_fac(int n) {
    vector<int> ret;
    if(n > 1) get_fac(n, ret);
    return ret;
}
inline ll PR(ll x) {
    auto f = [&](int x, int c, int n) {
        return ((lll) x * x + c) % n;
    };
    ll s = 0, t = 0, c = 1ll * rand() % (x - 1) + 1;
    int stp = 0, goal = 1;
    ll val = 1;
    for(goal = 1;; goal <<= 1, s = t, val = 1) {
        for(stp = 1; stp <= goal; ++stp) {
            t = f(t, c, x);
            val = (lll)val * abs(t - s) % x;
            if((stp % 127) == 0) {
                ll d = gcd(val, x);
                if(d > 1)
                    return d;
            }
        }
        ll d = gcd(val, x);
        if(d > 1)
            return d;
    }
}

inline void fac(ll x) {
    if(x <= max_factor || x < 2)
        return;
    if(is_prime(x)) {
        max_factor = max_factor > x ? max_factor : x;
        return;
    }
    ll p = x;
    while(p >= x)
        p = PR(x);
    while((x % p) == 0)
        x /= p;
    fac(x), fac(p);
}

void solve() {
    int n;
    cin >> n;
    int mx = 0;
    max_factor = 0;
    auto f = go_fac(n);
    // fac(n);
    // for(int i = 0; i < f.size(); i++) {
    //     mx = max(mx, f[i]);
    // }
    out(f);
    if(max_factor == n) {
        cout << "Prime\n";
    } else {
        cout << max_factor << '\n';
    }
}
signed main() {
    // cin.tie(0)->sync_with_stdio(false);
    int t;
    cin >> t;
    while(t--) solve();
}
