#define fp(i, a, b) for (int i = (a); i <= (b); i++)
#define fd(i, a, b) for (int i = (a); i >= (b); i--)
const int N = 3e5 + 5, mod = 998244353; // (N = 4 * n)

using ll = int64_t;
using Poly = vector<int>;
/*---------------------------------------------------------------------------*/
// 二次剩余
class Cipolla {
    int mod, I2{};
    using pll = pair<ll, ll>;
#define X first
#define Y second
    ll mul(ll a, ll b) const { return a * b % mod; }
    pll mul(pll a, pll b) const { return {(a.X * b.X + I2 * a.Y % mod * b.Y) % mod, (a.X * b.Y + a.Y * b.X) % mod}; }
    template<class T> T ksm(T a, int b, T x) { for (; b; b >>= 1, a = mul(a, a)) if (b & 1) x = mul(x, a); return x; }
public:
    Cipolla(int p = 0) : mod(p) {}
    pair<int, int> sqrt(int n) {
        int a = rand(), x;
        if (!(n %= mod)) return {0, 0};
        if (ksm(n, (mod - 1) >> 1, 1ll) == mod - 1) return {-1, -1};
        while (ksm(I2 = ((ll) a * a - n + mod) % mod, (mod - 1) >> 1, 1ll) == 1) a = rand();
        x = (int) ksm(pll{a, 1}, (mod + 1) >> 1, {1, 0}).X;
        if (2 * x > mod) x = mod - x;
        return {x, mod - x};
    }
#undef X
#undef Y
};
/*--------------------------------------Modular----------------------------------*/
#define mul(a, b) (ll(a) * (b) % mod)
#define add(a, b) (((a) += (b)) >= mod ? (a) -= mod : 0) // (a += b) %= P
#define dec(a, b) (((a) -= (b)) < 0 ? (a) += mod: 0)  // ((a -= b) += P) %= P
Poly getInv(int L) { Poly inv(L); inv[1] = 1; fp(i, 2, L - 1) inv[i] = mul((mod - mod / i), inv[mod % i]); return inv; }
int ksm(ll a, int b = mod - 2, ll x = 1) { for (; b; b >>= 1, a = a * a % mod) if (b & 1) x = x * a % mod; return x; }
auto inv = getInv(N); // NOLINT
/*----------------------------------NTT--------------------------------------*/
namespace NTT {
    const int g = 3;
    Poly Omega(int L) {
        int wn = ksm(g, mod / L);
        Poly w(L); w[L >> 1] = 1;
        fp(i, L / 2 + 1, L - 1) w[i] = mul(w[i - 1], wn);
        fd(i, L / 2 - 1, 1) w[i] = w[i << 1];
        return w;
    }
    auto W = Omega(1 << 20); // NOLINT
    void DIF(int *a, int n) {
        for (int k = n >> 1; k; k >>= 1)
            for (int i = 0, y; i < n; i += k << 1)
                for (int j = 0; j < k; ++j)
                    y = a[i + j + k], a[i + j + k] = mul(a[i + j] - y + mod, W[k + j]), add(a[i + j], y);
    }
    void IDIT(int *a, int n) {
        for (int k = 1; k < n; k <<= 1)
            for (int i = 0, x, y; i < n; i += k << 1)
                for (int j = 0; j < k; ++j)
                    x = a[i + j], y = mul(a[i + j + k], W[k + j]),
                    a[i + j + k] = x - y < 0 ? x - y + mod : x - y, add(a[i + j], y);
        int Inv = mod - (mod - 1) / n;
        fp(i, 0, n - 1) a[i] = mul(a[i], Inv);
        reverse(a + 1, a + n);
    }
}
/*-------------------------------------Polynomial 全家桶-----------------------------------*/
namespace Polynomial {
    // basic operator
    int norm(int n) { return 1 << (__lg(n - 1) + 1); }
    void norm(Poly &a) { if (!a.empty()) a.resize(norm(a.size()), 0); else a = {0}; }
    void DFT(Poly &a) { NTT::DIF(a.data(), a.size()); }
    void IDFT(Poly &a) { NTT::IDIT(a.data(), a.size()); }
    Poly &dot(Poly &a, Poly &b) { fp(i, 0, a.size() - 1) a[i] = mul(a[i], b[i]); return a; }
    
    // mul / div int
    Poly &operator*=(Poly &a, int b) { for (auto &x : a) x = mul(x, b); return a; }
    Poly operator*(Poly a, int b) { return a *= b; }
    Poly operator*(int a, Poly b) { return b * a; }
    Poly &operator/=(Poly &a, int b) { return a *= ksm(b); }
    Poly operator/(Poly a, int b) { return a /= b; }

    // Poly add / sub
    Poly &operator+=(Poly &a, Poly b) {
        a.resize(max(a.size(), b.size()));
        fp(i, 0, b.size() - 1) add(a[i], b[i]);
        return a;
    }
    Poly operator+(Poly a, Poly b) { return a += b; }
    Poly &operator-=(Poly &a, Poly b) {
        a.resize(max(a.size(), b.size()));
        fp(i, 0, b.size() - 1) dec(a[i], b[i]);
        return a;
    }
    Poly operator-(Poly a, Poly b) { return a -= b; }

    // Poly mul
    Poly operator*(Poly a, Poly b) {
        int n = a.size() + b.size() - 1, L = norm(n);
        if (a.size() <= 8 || b.size() <= 8) {
            Poly c(n);
            fp(i, 0, a.size() - 1) fp(j, 0, b.size() - 1)
                c[i + j] = (c[i + j] + (ll) a[i] * b[j]) % mod;
            return c;
        }
        a.resize(L), b.resize(L);
        DFT(a), DFT(b), dot(a, b), IDFT(a);
        return a.resize(n), a;
    }
    
    // Poly inv
    Poly Inv2k(Poly a) { // |a| = 2 ^ k
        int n = a.size(), m = n >> 1;
        if (n == 1) return {ksm(a[0])};
        Poly b = Inv2k(Poly(a.begin(), a.begin() + m)), c = b;
        b.resize(n), DFT(a), DFT(b), dot(a, b), IDFT(a);
        fp(i, 0, n - 1) a[i] = i < m ? 0 : mod - a[i];
        DFT(a), dot(a, b), IDFT(a);
        return move(c.begin(), c.end(), a.begin()), a;
    }
    Poly Inv(Poly a) {
        int n = a.size();
        norm(a), a = Inv2k(a);
        return a.resize(n), a;
    }

    // Poly div / mod
    Poly operator/(Poly a,Poly b){
        int k = a.size() - b.size() + 1;
        if (k < 0) return {0};
        reverse(a.begin(), a.end());
        reverse(b.begin(), b.end());
        b.resize(k), a = a * Inv(b);
        a.resize(k), reverse(a.begin(), a.end());
        return a;
    }
    pair<Poly, Poly> operator%(Poly a, const Poly& b) {
        Poly c = a / b;
        a -= b * c, a.resize(b.size() - 1);
        return {c, a};
    }
    
    // Poly calculus
    Poly deriv(Poly a) {
        fp(i, 1, a.size() - 1) a[i - 1] = mul(i, a[i]);
        return a.pop_back(), a;
    }
    Poly integ(Poly a) {
        a.push_back(0);
        fd(i, a.size() - 1, 1) a[i] = mul(inv[i], a[i - 1]);
        return a[0] = 0, a;
    }

    // Poly ln
    Poly Ln(Poly a) {
        int n = a.size();
        a = deriv(a) * Inv(a);
        return a.resize(n - 1), integ(a);
    }
    
    // Poly exp
    Poly Exp(Poly a) {
        int n = a.size(), k = norm(n);
        Poly b = {1}, c, d; a.resize(k);
        for (int L = 2; L <= k; L <<= 1) {
            d = b, b.resize(L), c = Ln(b), c.resize(L);
            fp(i, 0, L - 1) c[i] = a[i] - c[i] + (a[i] < c[i] ? mod : 0);
            add(c[0], 1), DFT(b), DFT(c), dot(b, c), IDFT(b);
            move(d.begin(), d.end(), b.begin());
        }
        return b.resize(n), b;
    }
    
    // Poly sqrt
    Poly Sqrt(Poly a) {
        int n = a.size(), k = norm(n); a.resize(k);
        Poly b = {(new Cipolla(mod))->sqrt(a[0]).first, 0}, c;
        for (int L = 2; L <= k; L <<= 1) {
            b.resize(L), c = Poly(a.begin(), a.begin() + L) * Inv2k(b);
            fp(i, L / 2, L - 1) b[i] = mul(c[i], (mod + 1) / 2);
        }
        return b.resize(n), b;
    }

    // Poly pow
    Poly Pow(Poly &a, int b) { return Exp(Ln(a) * b); } // a[0] = 1
    Poly Pow(Poly a, int b1, int b2) { // b1 = b % mod, b2 = b % phi(mod) and b >= n iff a[0] > 0
        int n = a.size(), d = 0, k;
        while (d < n && !a[d]) ++d;
        if ((ll) d * b1 >= n) return Poly(n);
        a.erase(a.begin(), a.begin() + d);
        k = ksm(a[0]), norm(a *= k);
        a = Pow(a, b1) * ksm(k, mod - 1 - b2);
        a.resize(n), d *= b1;
        fd(i, n - 1, 0) a[i] = i >= d ? a[i - d] : 0;
        return a;
    }
    
}
using namespace Polynomial;