// sieve first!!

namespace min25 {

ll n, sqr, w[N];
int c;
ll g[N];

inline int id(ll x) {
    return x ? x <= sqr ? c - x + 1 : n / x : 0;
}

void cal_g(ll n_) {  // 这里的 g 只能是一个没有系数的单项式！！不是极性函数f
    n = n_; sqr = sqrt(n_); c = 0;
    for (ll l = 1, r; l <= n; l = r + 1) {
        ll v = w[++c] = n / l; r = n / v;
        g[c] = (v - 1) % mod;
    }
    for (int i = 1; i <= prinum; i++) {
        int p = pri[i];
        if (1ll * p * p > n) break;
        for (int j = 1; 1ll * p * p <= w[j]; ++j)
            g[j] -= g[id(w[j] / p)] - g[id(p - 1)];
    }
}
int cal_s(int n, ll x, int y) {   // 用 g 的单项式乘以系数加起来。binom那个地方是填f(p^e)
    if(x <= pri[y]) return 0;
    int ans = (g[id(x)] - g[id(pri[y])] + mod) * n % mod;
    
    for(int i = y + 1; i <= prinum && pri[i] * pri[i] <= x; i++) {
        ll P = pri[i];
        for(int e = 1; P <= x; e++, P *= pri[i]) {
            ans = (ans + (ll) binom(n + e - 1, n - 1) * (cal_s(n, x / P, i) + (e != 1)) % mod) % mod;
        }
    }
    return ans;
}

}