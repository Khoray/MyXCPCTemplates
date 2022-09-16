#define fp(i, a, b) for (int i = a, i##_ = (b) + 1; i < i##_; ++i)
#define fd(i, a, b) for (int i = a, i##_ = (b) - 1; i > i##_; --i)
using namespace std;
#define int long long
using ll = int64_t;
using db = double;
using Poly = vector<ll>;
struct cp {
    db x, y;
    cp(db real = 0, db imag = 0) : x(real), y(imag){};
    cp operator+(cp b) const { return {x + b.x, y + b.y}; }
    cp operator-(cp b) const { return {x - b.x, y - b.y}; }
    cp operator*(cp b) const { return {x * b.x - y * b.y, x * b.y + y * b.x}; }
};
using vcp = vector<cp>;
const int mod = 1e9 + 7;
namespace FFT {
    const db pi = acos(-1);
    vcp Omega(int L) { // In order to reduce the accuracy error
        vcp w(L); w[1] = 1;
        for (int i = 2; i < L; i <<= 1) {
            auto w0 = w.begin() + i / 2, w1 = w.begin() + i;
            cp wn(cos(pi / i), sin(pi / i));
            for (int j = 0; j < i; j += 2)
                w1[j] = w0[j >> 1], w1[j + 1] = w1[j] * wn;
        }
        return w;
    }
    auto W = Omega(1 << 23); // NOLINT
    void DIF(cp *a, int n) {
        cp x, y;
        for (int k = n >> 1; k; k >>= 1)
            for (int i = 0; i < n; i += k << 1)
                for (int j = 0; j < k; ++j)
                    x = a[i + j], y = a[i + j + k],
                    a[i + j + k] = (x - y) *  W[k + j], a[i + j] = x + y;
    }
    void IDIT(cp *a, int n) {
        cp x, y;
        for (int k = 1; k < n; k <<= 1)
            for (int i = 0; i < n; i += k << 1)
                for (int j = 0; j < k; ++j)
                    x = a[i + j], y = a[i + j + k] * W[k + j],
                    a[i + j + k] = x - y, a[i + j] = x + y;
        const db Inv = 1. / n;
        fp(i, 0, n - 1) a[i].x *= Inv, a[i].y *= Inv;
        reverse(a + 1, a + n);
    }
}
namespace Polynomial {
    void DFT(vcp &a) { FFT::DIF(a.data(), a.size()); }
    void IDFT(vcp &a) { FFT::IDIT(a.data(), a.size()); }
    int norm(int n) { return 1 << (__lg(n - 1) + 1); }
     
    // Poly mul
    vcp &dot(vcp &a, vcp &b) { fp(i, 0, a.size() - 1) a[i] = a[i] * b[i]; return a; }
    Poly &operator+=(Poly &a, Poly b) {
        a.resize(max(a.size(), b.size()));
        fp(i, 0, a.size() - 1) a[i] += b[i];
        return a;
    };
    Poly &operator-=(Poly &a, Poly b) {
        a.resize(max(a.size(), b.size()));
        fp(i, 0, a.size() - 1) DEC(a[i], b[i]);
        return a;
    };
    Poly operator * (Poly A, Poly B) {
        int n = A.size() + B.size() - 1;
        int L = norm(n);
        Poly res(n); vcp a(L), b(L), c(L), d(L);
        fp(i, 0, A.size() - 1) a[i] = cp(A[i] & 0x7fff, A[i] >> 15);
        fp(i, 0, B.size() - 1) b[i] = cp(B[i] & 0x7fff, B[i] >> 15);
        FFT::DIF(a.data(), L), FFT::DIF(b.data(), L);
        for (int k = 1, i = 0, j = 0; k < L; j ^= k, k <<= 1)
            for (; i < k * 2; i++, j = i ^ (k - 1)) {
                c[i] = cp(a[i].x + a[j].x, a[i].y - a[j].y) * b[i] * 0.5,
                d[i] = cp(a[i].y + a[j].y, -a[i].x + a[j].x) * b[i] * 0.5;
            }
        FFT::IDIT(c.data(), L), FFT::IDIT(d.data(), L);
        for (int i = 0; i < n; i++) {
            ll x = c[i].x + 0.5, y = c[i].y + 0.5, z = d[i].x + 0.5, w = d[i].y + 0.5;
            x %= mod;
            y %= mod;
            z %= mod;
            w %= mod;
            res[i] = (x + ((y + z) << 15) % mod + (w << 30) % mod) % mod;
        }
        return res;
    }
    vcp operator * (vcp a, vcp b) {
        int n = a.size() + b.size() - 1;
        int L = norm(n);
        a.resize(L), b.resize(L);
        DFT(a), DFT(b);
        dot(a, b);
        IDFT(a);
        return a;
    }
}
using namespace Polynomial;