#define fp(i, a, b) for (int i = (a); i <= (b); i++)
#define fd(i, a, b) for (int i = (a); i >= (b); i--)
const int N = 3e5 + 5, mod = 998244353; // (N = 4 * n)

using ll = int64_t;
using Poly = vector<int>;
/*--------------------------------------Modular----------------------------------*/
#define MUL(a, b) (ll(a) * (b) % mod)
#define ADD(a, b) (((a) += (b)) >= mod ? (a) -= mod : 0) // (a += b) %= P
#define DEC(a, b) (((a) -= (b)) < 0 ? (a) += mod: 0)  // ((a -= b) += P) %= P
Poly getInv(int L) {
    Poly inv(L);
    inv[1] = 1;
    fp(i, 2, L - 1) inv[i] = MUL((mod - mod / i), inv[mod % i]);
    return inv;
}
int ksm(ll a, int b = mod - 2, ll x = 1) {
    for(; b; b >>= 1, a = a * a % mod) if(b & 1) x = x * a % mod;
    return x;
}
auto inv = getInv(N); // NOLINT
/*----------------------------------NTT--------------------------------------*/
namespace NTT {
const int g = 3;
Poly Omega(int L) {
    int wn = ksm(g, mod / L);
    Poly w(L);
    w[L >> 1] = 1;
    fp(i, L / 2 + 1, L - 1) w[i] = MUL(w[i - 1], wn);
    fd(i, L / 2 - 1, 1) w[i] = w[i << 1];
    return w;
}
auto W = Omega(1 << 20); // NOLINT
void DIF(int *a, int n) {
    for(int k = n >> 1; k; k >>= 1)
        for(int i = 0, y; i < n; i += k << 1)
            for(int j = 0; j < k; ++j)
                y = a[i + j + k], a[i + j + k] = MUL(a[i + j] - y + mod, W[k + j]), ADD(a[i + j], y);
}
void IDIT(int *a, int n) {
    for(int k = 1; k < n; k <<= 1)
        for(int i = 0, x, y; i < n; i += k << 1)
            for(int j = 0; j < k; ++j)
                x = a[i + j], y = MUL(a[i + j + k], W[k + j]),
                a[i + j + k] = x - y < 0 ? x - y + mod : x - y, ADD(a[i + j], y);
    int Inv = mod - (mod - 1) / n;
    fp(i, 0, n - 1) a[i] = MUL(a[i], Inv);
    reverse(a + 1, a + n);
}
}
/*-------------------------------------Polynomial 全家桶-----------------------------------*/
namespace Polynomial {
    // basic operator
    int norm(int n) {
        return 1 << (__lg(n - 1) + 1);
    }
    void norm(Poly &a) {
        if(!a.empty()) a.resize(norm(a.size()), 0);
        else a = {0};
    }
    void DFT(Poly &a) {
        NTT::DIF(a.data(), a.size());
    }
    void IDFT(Poly &a) {
        NTT::IDIT(a.data(), a.size());
    }
    Poly &dot(Poly &a, Poly &b) {
        fp(i, 0, a.size() - 1) a[i] = MUL(a[i], b[i]);
        return a;
    }

    // MUL / div int
    Poly &operator*=(Poly &a, int b) {
        for(auto &x : a) x = MUL(x, b);
        return a;
    }
    Poly operator*(Poly a, int b) {
        return a *= b;
    }
    Poly operator*(int a, Poly b) {
        return b * a;
    }
    Poly &operator/=(Poly &a, int b) {
        return a *= ksm(b);
    }
    Poly operator/(Poly a, int b) {
        return a /= b;
    }

    // Poly ADD / sub
    Poly &operator+=(Poly &a, Poly b) {
        a.resize(max(a.size(), b.size()));
        fp(i, 0, b.size() - 1) ADD(a[i], b[i]);
        return a;
    }
    Poly operator+(Poly a, Poly b) {
        return a += b;
    }
    Poly &operator-=(Poly &a, Poly b) {
        a.resize(max(a.size(), b.size()));
        fp(i, 0, b.size() - 1) DEC(a[i], b[i]);
        return a;
    }
    Poly operator-(Poly a, Poly b) {
        return a -= b;
    }

    // Poly MUL
    Poly operator*(Poly a, Poly b) {
        int n = a.size() + b.size() - 1, L = norm(n);
        if(a.size() <= 8 || b.size() <= 8) {
            Poly c(n);
            fp(i, 0, a.size() - 1) fp(j, 0, b.size() - 1)
            c[i + j] = (c[i + j] + (ll) a[i] * b[j]) % mod;
            return c;
        }
        a.resize(L), b.resize(L);
        DFT(a), DFT(b), dot(a, b), IDFT(a);
        return a.resize(n), a;
    }
}
using namespace Polynomial;