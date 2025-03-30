namespace Comb {
    using ll = long long;
    const int N = 1e7 + 5;
    const int mod = 1e9 + 7;
    int ksm(int a, int b = mod - 2, int MOD_KSM = mod) {
        int ret = 1;
        while(b) {
            if(b & 1) {
                ret = (ll) ret * a % MOD_KSM;
            }
            a = (ll) a * a % MOD_KSM;
            b >>= 1;
        }
        return ret;
    }
    int facinv[N], fac[N];
    void init_fac() {
        fac[0] = fac[1] = 1;
        for(int i = 2; i < N; i++) {
            fac[i] = (ll) fac[i - 1] * i % mod;
        }
        facinv[N - 1] = ksm(fac[N - 1], mod - 2);
        for(int i = N - 2; i >= 0; i--) {
            facinv[i] = (ll) facinv[i + 1] * (i + 1) % mod;
        }
    }
    int binom(int n, int k) {
        if(n < 0 || k < 0 || k > n) { return 0; }
        return (ll) fac[n] * facinv[n - k] % mod * facinv[k] % mod;
    }

    int inv[N];
    void init_inv() {
        inv[0] = inv[1] = 1;
        for(int i = 2; i < N; i++) {
            inv[i] = (ll) (mod - mod / i) * inv[mod % i] % mod;
        }
    }

    int _ = [] () { init_fac(); init_inv(); return 0; } ();
}

using namespace Comb;
