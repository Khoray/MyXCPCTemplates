// #pragma GCC optimize("-Ofast","-funroll-all-loops","-ffast-math")
// #pragma GCC optimize("-fno-math-errno")
// #pragma GCC optimize("-funsafe-math-optimizations")
// #pragma GCC optimize("-freciprocal-math")
// #pragma GCC optimize("-fno-trapping-math")
// #pragma GCC optimize("-ffinite-math-only")
// #pragma GCC optimize("-fno-stack-protector")
// #pragma GCC target ("avx2","sse4.2","fma")
// #include <immintrin.h>
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using lll = __int128_t;
using ull = unsigned long long;

namespace __POLY__ {
	#define inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) inline
	const int mod = 998244353;
	const int proot = 3, pool_siz = 1 << 23;
	typedef long long i64;
	typedef unsigned int u32;
	typedef unsigned long long u64;
	typedef vector<u32> vu32;
	namespace math {
		inline int qp(long long x, int y, int ans = 1) {
			for (y < 0 ? y += mod - 1 : 0; y; y >>= 1, x = x * x % mod)
				if (y&  1) ans = ans * x % mod;
			return ans;
		}
		inline constexpr int lg(u32 x) { return x == 0 ? 0 : ((int)sizeof(int) * __CHAR_BIT__ - __builtin_clz(x)); }
		inline u32 fst_mul(u32 x, u64 p, u64 q) { return x * p - (q * x >> 32) * mod; }
		inline u32 Norm(u32 m) { return m >= mod ? m - mod : m; }
		const u32 modm2 = mod + mod;
		vu32 __fac({1, 1}), __ifc({1, 1}), __inv({0, 1});
		inline void __prep(int n) {
			static int i = 2;
			if (i < n) for (__fac.resize(n), __ifc.resize(n), __inv.resize(n); i < n; i++)
				__fac[i] = 1ll * i * __fac[i - 1] % mod, __inv[i] = 1ll * (mod - mod / i) * __inv[mod % i] % mod, __ifc[i] = 1ll * __inv[i] * __ifc[i - 1] % mod;
		}
		inline u32 gfac(u32 x) { return __prep(x + 1), __fac[x]; }
		inline u32 gifc(u32 x) { return __prep(x + 1), __ifc[x]; }
		inline u32 ginv(u32 x) { return __prep(x + 1), __inv[x]; }
		inline u32 gC(u32 n, u32 m) {
			if (n < m) return 0;
			return 1ll * gfac(n) * gifc(m) % mod * gifc(n - m) % mod;
		}
		u32 I = 0;
		struct cpl {
			u32 x, y;
			cpl(u32 _x = 0, u32 _y = 0) : x(_x), y(_y) {}
			inline cpl operator*(const cpl& a) const { return cpl((1ull * x * a.x + 1ull * I * y % mod * a.y) % mod, (1ull * x * a.y + 1ull * y * a.x) % mod); }
		};
		inline cpl cplpow(cpl a, int y, cpl b = cpl(1, 0)) {
			for (; y; y >>= 1, a = a * a) if (y&  1) b = b * a;
			return b;
		}
		inline u32 isqrt(u32 x) {
			static mt19937 rnd(mod);
			if (mod == 2 || !x || x == 1) return x;
			u32 a = 0;
			do {
				a = rnd() % mod;
			} while (qp((1ull * a * a + mod - x) % mod, mod >> 1) != mod - 1);
			I = (1ll * a * a + mod - x) % mod;
			a = cplpow(cpl(a, 1), (mod + 1) >> 1).x;
			return min(a, mod - a);
		}
	} using namespace math;
	namespace polynomial {
		const int maxbit = 23;
		namespace fast_number_theory_transform {
			u32 pool[(pool_siz) * 4] __attribute__((aligned(64))), *ptr = pool;
			u32 *p0[(pool_siz)], *p1[(pool_siz)], *q0[(pool_siz)], *q1[(pool_siz)];
			inline void bit_flip(u32 *p, int t) {
				for (int i = 0, j = 0; i < t; ++i) {
					if (i > j) swap(p[i], p[j]);
					for (int l = t >> 1; (j ^= l) < l; l >>= 1) ;
				}
			}
			inline void prep(int n) {
				static int t = 1;
				for (; t < n; t <<= 1) {
					int g = qp(proot, (mod - 1) / (t * 2));
					u32 *p, *q;
					p = p0[t] = ptr;
					ptr += max(t, 16);
					p[0] = 1;
					for (int m = 1; m < t; ++m)
						p[m] = p[m - 1] * (u64)g % u32(mod);
					bit_flip(p, t);
					q = q0[t] = ptr;
					ptr += max(t, 16);
					for (int i = 0; i < t; ++i)
						q[i] = (u64(p[i]) << 32) / mod;
					g = qp(g, mod - 2);
					p = p1[t] = ptr;
					ptr += max(t, 16);
					p[0] = 1;
					for (int m = 1; m < t; ++m)
						p[m] = p[m - 1] * (u64)g % u32(mod);
					bit_flip(p, t);
					q = q1[t] = ptr;
					ptr += max(t, 16);
					for (int i = 0; i < t; ++i)
						q[i] = (u64(p[i]) << 32) / mod;
				}
			}
			inline u32 my_mul(u32 a, u32 b, u32 c) { return b * (u64)a - ((u64(a) * c) >> 32) * u64(mod); }
			inline __m128i my_mullo_epu32(const __m128i& a, const __m128i& b) { return (__m128i)((__v4su)a * (__v4su)b); }
			inline __m128i my_mulhi_epu32(const __m128i& a, const __m128i& b) {
				__m128i a13 = _mm_shuffle_epi32(a, 0xF5);
				__m128i b13 = _mm_shuffle_epi32(b, 0xF5);
				__m128i prod02 = _mm_mul_epu32(a, b);
				__m128i prod13 = _mm_mul_epu32(a13, b13);
				__m128i prod01 = _mm_unpacklo_epi32(prod02, prod13);
				__m128i prod23 = _mm_unpackhi_epi32(prod02, prod13);
				__m128i prod = _mm_unpackhi_epi64(prod01, prod23);
				return prod;
			}
			void ntt(u32 *__restrict__ x, int bit) {
				int n = 1 << bit, t = n;
				prep(n);
				for (int m = 1; m < n; m <<= 1) {
					t >>= 1;
					u32 *__restrict__ p = p0[m];
					u32 *__restrict__ q = q0[m];
					if (t == 1 or t == 2) {
						u32 *xa = x, *xb = x + t;
						for (int i = 0; i < m; ++i, xa += t + t, xb += t + t)
							for (int j = 0; j < t; ++j) {
								u32 u = xa[j] - (xa[j] >= modm2) * modm2;
								u32 v = my_mul(xb[j], p[i], q[i]);
								xa[j] = u + v;
								xb[j] = u - v + modm2;
							}
					}
					else if (t == 4) {
						u32 *xa = x, *xb = x + t;
						for (int i = 0; i < m; ++i, xa += t + t, xb += t + t) {
							const __m128i p4 = _mm_set1_epi32(p[i]), q4 = _mm_set1_epi32(q[i]), mm = _mm_set1_epi32(mod + mod), m0 = _mm_set1_epi32(0), m1 = _mm_set1_epi32(mod);
							for (int j = 0; j < t; j += 4) {
								__m128i u = _mm_loadu_si128((__m128i *)(xa + j));
								u = _mm_sub_epi32(u, _mm_and_si128(_mm_or_si128(_mm_cmpgt_epi32(u, mm), _mm_cmpgt_epi32(m0, u)), mm));
								__m128i v = _mm_loadu_si128((__m128i *)(xb + j));
								v = _mm_sub_epi32(my_mullo_epu32(v, p4), my_mullo_epu32(my_mulhi_epu32(v, q4), m1));
								_mm_storeu_si128((__m128i *)(xa + j), _mm_add_epi32(u, v));
								_mm_storeu_si128((__m128i *)(xb + j), _mm_add_epi32(_mm_sub_epi32(u, v), mm));
							}
						}
					}
					else {
						u32 *xa = x, *xb = x + t;
						for (int i = 0; i < m; ++i, xa += t + t, xb += t + t) {
							const __m128i p4 = _mm_set1_epi32(p[i]), q4 = _mm_set1_epi32(q[i]), mm = _mm_set1_epi32(mod + mod), m0 = _mm_set1_epi32(0), m1 = _mm_set1_epi32(mod);
							for (int j = 0; j < t; j += 8) {
								__m128i u0 = _mm_loadu_si128((__m128i *)(xa + j));
								__m128i u1 = _mm_loadu_si128((__m128i *)(xa + j + 4));
								__m128i v0 = _mm_loadu_si128((__m128i *)(xb + j));
								__m128i v1 = _mm_loadu_si128((__m128i *)(xb + j + 4));
								u0 = _mm_sub_epi32(u0, _mm_and_si128(_mm_or_si128(_mm_cmpgt_epi32(u0, mm), _mm_cmpgt_epi32(m0, u0)), mm));
								u1 = _mm_sub_epi32(u1, _mm_and_si128(_mm_or_si128(_mm_cmpgt_epi32(u1, mm), _mm_cmpgt_epi32(m0, u1)), mm));
								v0 = _mm_sub_epi32(my_mullo_epu32(v0, p4), my_mullo_epu32(my_mulhi_epu32(v0, q4), m1));
								v1 = _mm_sub_epi32(my_mullo_epu32(v1, p4), my_mullo_epu32(my_mulhi_epu32(v1, q4), m1));
								_mm_storeu_si128((__m128i *)(xa + j), _mm_add_epi32(u0, v0));
								_mm_storeu_si128((__m128i *)(xa + j + 4), _mm_add_epi32(u1, v1));
								_mm_storeu_si128((__m128i *)(xb + j), _mm_add_epi32(_mm_sub_epi32(u0, v0), mm));
								_mm_storeu_si128((__m128i *)(xb + j + 4), _mm_add_epi32(_mm_sub_epi32(u1, v1), mm));
							}
						}
					}
				}
				for (int i = 0; i < n; ++i) x[i] -= (x[i] >= modm2) * modm2, x[i] -= (x[i] >= u32(mod)) * u32(mod);
			}
			void intt(u32 *__restrict__ x, int bit) {
				int n = 1 << bit, t = 1;
				prep(n);
				for (int m = (n >> 1); m; m >>= 1) {
					u32 *__restrict__ p = p1[m];
					u32 *__restrict__ q = q1[m];
					if (t == 1 or t == 2) {
						u32 *xa = x, *xb = x + t;
						for (int i = 0; i < m; ++i, xa += t + t, xb += t + t)
							for (int j = 0; j < t; ++j) {
								u32 u = xa[j], v = xb[j];
								xa[j] = u + v - (u + v >= modm2) * modm2;
								xb[j] = my_mul(u - v + modm2, p[i], q[i]);
							}
					} else if (t == 4) {
						u32 *xa = x, *xb = x + t;
						for (int i = 0; i < m; ++i, xa += t + t, xb += t + t) {
							const __m128i p4 = _mm_set1_epi32(p[i]), q4 = _mm_set1_epi32(q[i]), mm = _mm_set1_epi32(mod + mod), m0 = _mm_set1_epi32(0), m1 = _mm_set1_epi32(mod);
							for (int j = 0; j < t; j += 4) {
								__m128i u = _mm_loadu_si128((__m128i *)(xa + j));
								__m128i v = _mm_loadu_si128((__m128i *)(xb + j));
								__m128i uv = _mm_add_epi32(u, v);
								_mm_storeu_si128((__m128i *)(xa + j), _mm_sub_epi32(uv, _mm_and_si128(_mm_or_si128(_mm_cmpgt_epi32(uv, mm), _mm_cmpgt_epi32(m0, uv)), mm)));
								uv = _mm_add_epi32(_mm_sub_epi32(u, v), mm);
								_mm_storeu_si128((__m128i *)(xb + j), _mm_sub_epi32(my_mullo_epu32(uv, p4), my_mullo_epu32(my_mulhi_epu32(uv, q4), m1)));
							}
						}
					} else {
						u32 *xa = x, *xb = x + t;
						for (int i = 0; i < m; ++i, xa += t + t, xb += t + t) {
							const __m128i p4 = _mm_set1_epi32(p[i]), q4 = _mm_set1_epi32(q[i]), mm = _mm_set1_epi32(mod + mod), m0 = _mm_set1_epi32(0), m1 = _mm_set1_epi32(mod);
							for (int j = 0; j < t; j += 8) {
								__m128i u0 = _mm_loadu_si128((__m128i *)(xa + j));
								__m128i u1 = _mm_loadu_si128((__m128i *)(xa + j + 4));
								__m128i v0 = _mm_loadu_si128((__m128i *)(xb + j));
								__m128i v1 = _mm_loadu_si128((__m128i *)(xb + j + 4));
								__m128i uv0 = _mm_add_epi32(u0, v0);
								__m128i uv1 = _mm_add_epi32(u1, v1);
								_mm_storeu_si128((__m128i *)(xa + j), _mm_sub_epi32(uv0, _mm_and_si128(_mm_or_si128(_mm_cmpgt_epi32(uv0, mm), _mm_cmpgt_epi32(m0, uv0)), mm)));
								_mm_storeu_si128((__m128i *)(xa + j + 4), _mm_sub_epi32(uv1, _mm_and_si128(_mm_or_si128(_mm_cmpgt_epi32(uv1, mm), _mm_cmpgt_epi32(m0, uv1)), mm)));
								uv0 = _mm_add_epi32(_mm_sub_epi32(u0, v0), mm);
								uv1 = _mm_add_epi32(_mm_sub_epi32(u1, v1), mm);
								_mm_storeu_si128((__m128i *)(xb + j), _mm_sub_epi32(my_mullo_epu32(uv0, p4), my_mullo_epu32(my_mulhi_epu32(uv0, q4), m1)));
								_mm_storeu_si128((__m128i *)(xb + j + 4), _mm_sub_epi32(my_mullo_epu32(uv1, p4), my_mullo_epu32(my_mulhi_epu32(uv1, q4), m1)));
							}
						}
					} t <<= 1;
				}
				u32 rn = qp(n, mod - 2);
				for (int i = 0; i < n; ++i) x[i] = x[i] * (u64)rn % mod;
			}
		}
		using fast_number_theory_transform::intt;
		using fast_number_theory_transform::ntt;
		struct poly {
			vu32 f;
			template <typename _Tp = size_t, typename _Tv = u32>
			poly(_Tp len = 1, _Tv same_val = 0) : f(len, same_val) {}
			poly(const vu32& _f) : f(_f) {}
			poly(const vector<int>& _f) {
				f.resize(((int)_f.size()));
				for (int i = 0; i < ((int)_f.size()); i++) 
					f[i] = _f[i] + ((_f[i] >> 31)&  mod);
			}
			template <typename T> poly(initializer_list<T> _f) : poly(vector<T>(_f)) {}
			template <typename T> poly(T *__first, T *__last) : poly(vector<typename iterator_traits<T>::value_type>(__first, __last)) {}
			inline operator vu32() const { return f; }
			inline vu32::iterator begin() { return f.begin(); }
			inline vu32::iterator end() { return f.end(); }
			inline const vu32::const_iterator begin() const { return f.begin(); }
			inline const vu32::const_iterator end() const { return f.end(); }
			inline void swap(poly& _f) { f.swap(_f.f); }
			inline int degree() const { return (int)f.size() - 1; }
			inline int size() const { return (int)f.size(); }
			inline poly& resize(int x) { return f.resize(x), *this; }
			inline poly& redegree(int x) { return f.resize(x + 1), *this; }
			inline void clear() { f.resize(1), f[0] = 0; }
			inline void shrink() { int ndeg = f.size() - 1; while (ndeg > 0 and f[ndeg] == 0) ndeg--; f.resize(ndeg + 1); }
			inline void rev() { reverse(f.begin(), f.end()); }
			inline poly split(int n) const { return n <= 0 ? poly(1, 1) : (n < (int)f.size() ? poly(f.begin(), f.begin() + n + 1) : poly(*this).redegree(n)); }
			inline u32& operator[](u32 x) { return f[x]; }
			inline u32 operator[](u32 x) const { return f[x]; }
			inline u32 get(u32 x) const { return x < f.size() ? f[x] : 0; }
			inline int eval(int x) { int ret = f[degree()]; for(int i = degree() - 1; i >= 0; -- i) ret = (1ll * ret * x + f[i]) % mod; return ret; }
			inline friend istream& operator>>(istream& in, poly& x) {
				for (int i = 0, _buf; i < x.size(); i++) in >> _buf, _buf %= mod, _buf += (_buf < 0) * mod, x[i] = _buf; 
				return in;
			}
			inline friend ostream& operator<<(ostream& out, const poly& x) {
				out << x[0];
				for (int i = 1; i < x.size(); i++) out << ' ' << x[i];
				return out;
			}
			inline u32 *data() { return f.data(); }
			inline const u32 *data() const { return f.data(); }
			inline poly& operator+=(const poly& a) {
				f.resize(max(f.size(), a.f.size()));
				for (int i = 0; i < a.f.size(); i++) f[i] = f[i] + a.f[i] - (f[i] + a.f[i] >= mod) * mod;
				return *this;
			}
			inline poly& operator-=(const poly& a) {
				f.resize(max(f.size(), a.f.size()));
				for (int i = 0; i < a.f.size(); i++) f[i] = f[i] - a.f[i] + (f[i] < a.f[i]) * mod;
				return *this;
			}
			inline poly& operator+=(const u32& b) { f[0] = f[0] + b - mod * (f[0] + b >= mod); return *this; }
			inline poly& operator-=(const u32& b) { f[0] = f[0] - b + mod * (f[0] < b); return *this; }
			inline poly operator+(const poly& a) const { return (poly(*this) += a); }
			inline poly operator-(const poly& a) const { return (poly(*this) -= a); }
			friend inline poly operator+(u32 a, const poly& b) { return (poly(1, a) += b); }
			friend inline poly operator-(u32 a, const poly& b) { return (poly(1, a) -= b); }
			friend inline poly operator+(const poly& a, u32 b) { return (poly(a) += poly(1, b)); }
			friend inline poly operator-(const poly& a, u32 b) { return (poly(a) -= poly(1, b)); }
			inline poly operator-() const {
				poly _f;
				_f.f.resize(f.size());
				for (int i = 0; i < _f.f.size(); i++) _f.f[i] = (f[i] != 0) * mod - f[i];
				return _f;
			}
			inline poly shiftvar(int k) const {
				poly ret(size());
				for (int i = 0; i * k <= degree(); ++i) ret[i * k] = f[i];
				return ret;
			}
			inline poly amp(int k) const {
				poly ret(size());
				for (int i = 0; i * k <= degree(); ++i) ret[i * k] = f[i];
				return ret;
			}
			inline poly& operator*=(const poly& a) {
				int n = degree(), m = a.degree();
				if (n <= 32 || m <= 32) {
					f.resize(n + m + 1);
					for (int i = n + m; i >= 0; i--) {
						f[i] = 1ll * f[i] * a.f[0] % mod;
						for (int j = max(1, i - n), j_up = min(m, i); j <= j_up; j++) f[i] = (f[i] + 1ll * f[i - j] * a.f[j]) % mod;
					} return *this;
				}
				vu32 _f(a.f);
				int bit = lg(n + m);
				f.resize(1 << bit), _f.resize(1 << bit);
				ntt(f.data(), bit), ntt(_f.data(), bit);
				for (int i = 0; i < (1 << bit); i++) f[i] = 1ll * f[i] * _f[i] % mod;
				intt(f.data(), bit), f.resize(n + m + 1);
				return *this;
			}
			inline poly operator*(const poly& a) const { return (poly(*this) *= a); }
			template <typename T> inline friend poly operator*(const poly& a, const T& b) {
				poly ret(a);
				for (int i = 0; i < ret.f.size(); ++i) ret[i] = 1ll * ret[i] * b % mod;
				return ret;
			}
			template <typename T> inline friend poly operator*(const T& b, const poly& a) {
				poly ret(a);
				for (int i = 0; i < ret.f.size(); ++i) ret[i] = 1ll * ret[i] * b % mod;
				return ret;
			}
			template <typename T> inline poly& operator*=(const T& b) { for (int i = 0; i < f.size(); ++i) f[i] = 1ll * f[i] * b % mod; return *this; }
			inline poly& operator>>=(int x) { return f.resize(f.size() + x), memmove(f.data() + x, f.data(), 4 * (f.size() - x)), memset(f.data(), 0, 4 * x), *this; }
			inline poly operator>>(int x) const { return (poly(*this) >>= x); }
			inline poly& operator<<=(int x) { return x >= f.size() ? (clear(), *this) : (memmove(f.data(), f.data() + x, 4 * (f.size() - x)), f.resize(f.size() - x), *this); }
			inline poly operator<<(int x) const { return (poly(*this) <<= x); }
			inline poly& shiftindexwith(int x) { return x >= f.size() ? (memset(f.data(), 0, 4 * f.size()), *this) : (memmove(f.data(), f.data() + x, 4 * (f.size() - x)), memset(f.data(), 0, 4 * x), *this); }
			inline poly shiftindex(int x) const { return (poly(*this).shiftindexwith(x)); }
			inline poly inv() const;
			inline poly quo(const poly& g) const;
			inline poly operator/(const poly& g) { return f.size() == 1 ? poly(1, qp(g[0], -1, f[0])) : quo(g); }
			inline poly& quowith(const poly& g) { return f.size() == 1 ? (f[0] = qp(g[0], -1, f[0]), *this) : (*this = quo(g)); }
			friend pair<poly,poly> operator % (const poly&  x, const poly&  y) {
				poly A = x, B = y;
				A.rev(), B.rev(); B.resize(x.degree() - y.degree() + 1);
				B = B.inv(); B *= A;
				B.resize(x.degree() - y.degree() + 1); B.rev();
				poly R = x - B * y; R.resize(y.degree() - 1);
				return {B, R};
			}
			inline poly deri() const {
				int n = degree(); poly res(n);
				for (int i = 1; i <= n; i++) res[i - 1] = 1ll * f[i] * i % mod;
				return res;
			}
			inline poly intg(u32 C = 0) const {
				int n = degree(); poly res(n + 2); res[0] = C;
				for (int i = 0; i <= n; i++) res[i + 1] = 1ll * ginv(i + 1) * f[i] % mod;
				return res;
			}
			inline poly pow(u32 x, u32 modphix = -1) {
				if (modphix == -1) modphix = x;
				int n = size() - 1;
				long long empt = 0;
				while (empt <= n and !f[empt]) ++empt;
				if (1ll * empt * x > n) return poly(size());
				poly res(size());
				for (int i = 0; i <= n - empt; ++i) res[i] = f[i + empt];
				int val_0 = res[0], inv_0 = qp(val_0, mod - 2), pow_0 = qp(val_0, modphix);
				for (int i = 0; i <= n - empt; ++i) res[i] = 1ll * res[i] * inv_0 % mod;
				res = (res.ln() * x).exp();
				empt *= x;
				for (int i = n; i >= empt; --i) res[i] = 1ll * res[i - empt] * pow_0 % mod;
				for (int i = empt - 1; i >= 0; --i) res[i] = 0;
				return res;
			}
			inline poly ivsqrt() const {
				int nsize = f.size(), mxb = lg(f.size() - 1);
				vu32 a(1 << mxb), _f(f);
				_f.resize(1 << mxb);
				a[0] = qp(isqrt(f[0]), mod - 2);
				for (int nb = 0; nb < mxb; nb++) {
					vu32 _a(a.begin(), a.begin() + (1 << nb)), _b(_f.begin(), _f.begin() + (2 << nb));
					_a.resize(4 << nb), _b.resize(4 << nb);
					ntt(_a.data(), nb + 2), ntt(_b.data(), nb + 2);
					for (int i = 0; i < (4 << nb); i++)
						_a[i] = 1ull * (mod - _a[i]) * _a[i] % mod * _a[i] % mod * _b[i] % mod, _a[i] = (_a[i] + (_a[i]&  1) * mod) >> 1;
					intt(_a.data(), nb + 2), memcpy(a.data() + (1 << nb), _a.data() + (1 << nb), 4 << nb);
				}
				return a.resize(nsize), a;
			}
			inline poly sqrt() const {
				if (f.size() == 1) return poly(1, isqrt(f[0]));
				if (f.size() == 2 and f[0] == 1)
					return poly(vector<int>{1, (int)(1ll * f[1] * (mod + 1) / 2 % mod)});
				int nsize = f.size(), mxb = lg(nsize - 1);
				vu32 a(1 << mxb), _f(f), _b;
				_f.resize(1 << mxb);
				a[0] = qp(isqrt(f[0]), mod - 2);
				for (int nb = 0; nb < mxb - 1; nb++) {
					vu32 _a(a.begin(), a.begin() + (1 << nb));
					_b = vu32(_f.begin(), _f.begin() + (2 << nb));
					_a.resize(4 << nb), _b.resize(4 << nb);
					ntt(_a.data(), nb + 2), ntt(_b.data(), nb + 2);
					for (int i = 0; i < (4 << nb); i++)
						_a[i] = 1ull * (mod - _a[i]) * _a[i] % mod * _a[i] % mod * _b[i] % mod, _a[i] = (_a[i] + (_a[i]&  1) * mod) >> 1;
					intt(_a.data(), nb + 2);
					memcpy(a.data() + (1 << nb), _a.data() + (1 << nb), 4 << nb);
				}
				ntt(a.data(), mxb);
				vu32 _a(a);
				for (int i = 0; i < (1 << mxb); i++) a[i] = 1ll * a[i] * _b[i] % mod;
				intt(a.data(), mxb), memset(a.data() + (1 << (mxb - 1)), 0, 2 << mxb);
				vu32 g0(a);
				ntt(a.data(), mxb), ntt(_f.data(), mxb);
				for (int i = 0; i < (1 << mxb); i++)
					a[i] = (1ll * a[i] * a[i] + mod - _f[i]) % mod * (mod - _a[i]) % mod, a[i] = (a[i] + (a[i]&  1) * mod) >> 1;
				intt(a.data(), mxb);
				memcpy(g0.data() + (1 << (mxb - 1)), a.data() + (1 << (mxb - 1)), 2 << mxb);
				return g0;
			}
			inline poly czt(int c, int m) const {
				poly ret(f);
				int inv = qp(c, mod - 2), n = ret.size();
				ret.resize(m);
				poly F(n), G(n + m);
				for (int i = 0, p1 = 1, p2 = 1; i < n; ++i) {
					F[n - i - 1] = 1ll * ret[i] * p1 % mod;
					if (i > 0) p2 = 1ll * p2 * inv % mod, p1 = 1ll * p1 * p2 % mod;
				}
				for (int i = 0, p1 = 1, p2 = 1; i < n + m; ++i) {
					G[i] = p1;
					if (i > 0) p2 = 1ll * p2 * c % mod, p1 = 1ll * p1 * p2 % mod;
				}
				F = F * G;
				for (int i = 0, p1 = 1, p2 = 1; i < m; ++i) {
					ret[i] = 1ll * F[i + n - 1] * p1 % mod;
					if (i > 0) p2 = 1ll * p2 * inv % mod, p1 = 1ll * p1 * p2 % mod;
				}
				return ret;
			}
			inline poly ChirpZ(int c, int m) const { return czt(c, m); }
			inline poly shift(int c) const {
				c %= mod;
				c = c + (c < 0) * mod;
				if (c == 0) return *this;
				poly A(size()), B(size()), ret(size());
				for (int i = 0; i < size(); ++i) A[size() - i - 1] = 1ll * f[i] * gfac(i) % mod;
				for (int i = 0, pc = 1; i < size(); ++i, pc = 1ll * pc * c % mod)
					B[i] = 1ll * pc * gifc(i) % mod;
				A *= B, A.resize(size());
				for (int i = 0; i < size(); ++i) ret[i] = 1ll * A[size() - i - 1] * gifc(i) % mod;
				return ret;
			}
			inline poly fdt() const {
				poly F(*this), E(size());
				for (int i = 0; i < size(); ++i) E[i] = gifc(i);
				F *= E, F.resize(size());
				for (int i = 0; i < size(); ++i) F[i] = 1ll * F[i] * gfac(i) % mod;
				return F;
			}
			inline poly ifdt() const {
				poly F(*this), E(size());
				for (int i = 0; i < size(); ++i) F[i] = 1ll * F[i] * gifc(i) % mod;
				for (int i = 0; i < size(); ++i)
					if (i&  1) E[i] = mod - gifc(i);
					else E[i] = gifc(i);
				return (F * E).split(degree());
			}
			inline poly ln() const;
			inline poly exp() const;
			inline poly eval(int n, poly a) const;
			inline poly intp(int n, const poly& x, const poly& y);
			inline poly sin() const {
				int omega_4 = qp(proot, (mod - 1) >> 2);
				poly F = ((*this) * omega_4).exp();
				return qp(omega_4 * 2, mod - 2) * (F - F.inv());
			}
			inline poly cos() const {
				int omega_4 = qp(proot, (mod - 1) >> 2);
				poly F = ((*this) * omega_4).exp();
				return qp(2, mod - 2) * (F + F.inv());
			}
			inline poly tan() const { return sin() / cos(); }
			inline poly asin() const {
				poly A = deri(), B = (*this) * (*this);
				B.resize(size());
				B = (1 - B).ivsqrt();
				return (A * B).intg().split(degree());
			}
			inline poly acos() const {
				poly A = (mod - 1) * deri(), B = (*this) * (*this);
				B.resize(size());
				B = (1 - B).ivsqrt();
				return (A * B).intg().split(degree());
			}
			inline poly atan() const {
				poly A = deri(), B = 1 + (*this) * (*this);
				B.resize(size());
				B = B.inv();
				return (A * B).intg().split(degree());
			}
			inline poly composite_inv(int n = -1) const {
				if (n == -1) n = size();
				auto enum_kth = [&](const poly& f, const poly& g, int k, int n) {
					/*return f(y) = [x^k](g(x) / (1 - y* f(x))) = \sum_{i = 0}^{n - 1} [x^k] g(x) f^i(x) y^i*/
					if (k < 0 or n <= 0) return poly();
					poly P(k + 1), Q((k + 1) << 1);
					copy_n(g.f.cbegin(), min(P.size(), g.size()), P.f.begin());
					Q.f.front() = 1;
					if (f.size()) for (int i = k + 1, j = 0; i < Q.size() and j < f.size();) Q[i ++] = (f[j] == 0 ? 0 : mod - f[j]), ++ j;
					auto quad_nonres = [&](){ for(int i = 2; ; ++ i) if (qp(i, mod >> 1) == mod - 1) return i; };
					auto sylow2_subgroup_gen = [&](){ return qp(quad_nonres(), mod >> __builtin_ctz(mod - 1)); };
					auto get_root = [&](int n) {
						vu32 root = {ginv(2)};
						vector<int> irt(__builtin_ctz(mod - 1) - 1);
						irt.back() = qp(sylow2_subgroup_gen(), mod - 2);
						for(int i = __builtin_ctz(mod - 1) - 3; i >= 0; -- i) irt[i] = 1ll * irt[i + 1] * irt[i + 1] % mod;
						int s = (int)root.size();
						if (s < n) {
							root.resize(n);
							for (int i = __builtin_ctz(s), j; (1 << i) < n; ++ i) {
								root[j = (1 << i)] = irt[i];
								for (int k = j + 1; k < (j << 1); ++ k)
									root[k] = 1ll * root[k - j] * root[j] % mod;
								root[j] = 1ll * root[j] * root.front() % mod;
							}
						} return root;
					};
					for (int d = 1; k != 0; d <<= 1, k >>= 1) {
						const int lg_len = lg((2 * d + 1) * (2 * k + 2) - 1), len = 1 << lg_len; 
						poly P_(len), Q_(len), U(len / 2), V(len / 2);
						for (int i = 0; i <  d; ++ i) copy_n(P.f.cbegin() + i * (k + 1), k + 1, P_.f.begin() + i * (2 * k + 2));
						for (int i = 0; i <= d; ++ i) copy_n(Q.f.cbegin() + i * (k + 1), k + 1, Q_.f.begin() + i * (2 * k + 2)); 
						ntt(P_.data(), lg_len); ntt(Q_.data(), lg_len);
						if (k&  1) {
							auto root = get_root(len >> 1);
							for (int i = 0; i < len; i += 2) {
								U[i / 2] = 1ll * (1ll * P_[i] * Q_[i + 1] % mod - 1ll * P_[i + 1] * Q_[i] % mod + mod) * root[i / 2] % mod;
								V[i / 2] = 1ll * Q_[i] * Q_[i + 1] % mod;
							}
						} else {
							auto root = get_root(1);
							for (int i = 0; i < len; i += 2) {
								U[i / 2] = 1ll * (1ll * P_[i] * Q_[i + 1] + 1ll * P_[i + 1] * Q_[i]) % mod * root[0] % mod;
								V[i / 2] = 1ll * Q_[i] * Q_[i + 1] % mod;
							}
						} 
						intt(U.data(), lg_len - 1), intt(V.data(), lg_len - 1);
						P.f.assign((2 * d) * (k / 2 + 1), 0);
						Q.f.assign((2 * d + 1) * (k / 2 + 1), 0);
						for (int i = 0; i <  (d << 1); ++ i) copy_n(U.f.cbegin() + i * (k + 1), k / 2 + 1, P.f.begin() + i * (k / 2 + 1));
						for (int i = 0; i <= (d << 1); ++ i) copy_n(V.f.cbegin() + i * (k + 1), k / 2 + 1, Q.f.begin() + i * (k / 2 + 1)); 
					} P.resize(n), Q.resize(n);
					return (P / Q).resize(n);
				};
				if (n <= 0 or f.size() < 2) return poly(0);
				if (n == 1) return poly(1);
				poly F = *this; F.resize(n);
				int f1_inv = qp(F[1], mod - 2), _c = f1_inv;
				for (int i = 1; i < n; ++ i) F[i] = 1ll * F[i] * _c % mod, _c = 1ll * _c * f1_inv % mod;
			   
				auto a = enum_kth(F, (poly){1}, n - 1, n);
				for (int i = 1; i < n; ++ i) a[i] = 1ll * a[i] * (n - 1) % mod * ginv(i) % mod;
				poly a_(a.size());
				for (int i = 0; i < a.size(); ++ i)  a_[i] = a[a.degree() - i];
				a_ = a_.pow(mod - qp(n - 1, mod - 2));
				poly B(2); B[0] = 0, B[1] = f1_inv;
				return (a_ * B).resize(n);
			} 
			inline poly composite(const poly& g, int n = -1) const {
				if (n == -1) n = size();
				if (n <= 0) return poly();
				if (g.size() == 0) return poly(n + 1);
				poly Q(n * 2);
				int g0_ = g[0];
				Q[0] = 1; 
				for (int i = n, j = 0; j < g.size() and i < 2 * n;) Q[i ++] = (g[j] == 0 ? 0 : mod - g[j]), ++ j;
				function<poly(const poly&, int, int)> rec = [&](const poly& Q, int d, int n) {
					if (n == 0) {
						poly P(d), Qinv(d);
						for(int i = d - 1, j = 0; j < f.size() and i >= 0; ) P[i --] = f[j ++];
						for(int i = 0, e = 1; i < d; ++ i) Qinv[i] = 1ll * gC(d + i - 1, i) * e % mod, e = 1ll * e * g0_ % mod;
						return (P * Qinv).resize(d);
					}
					const int lg_len = lg((2 * d + 1) * (2 * n + 2) - 1), len = 1 << lg_len;
					poly Q_(len), VV(1 << (lg_len - 1));
					for (int i = 0; i <= d; ++ i) copy_n(Q.f.begin() + i * (n + 1), n + 1, Q_.f.begin() + i * (n * 2 + 2));
					ntt(Q_.data(), lg_len);
					for (int i = 0; i < len; i += 2) VV[i / 2] = 1ll * Q_[i] * Q_[i + 1] % mod;
					intt(VV.data(), lg_len - 1);
					poly V((d * 2 + 1) * (n / 2 + 1));
					for (int i = 0; i <= 2 * d; ++ i) copy_n(VV.f.begin() + i * (n + 1), n / 2 + 1, V.f.begin() + i * (n / 2 + 1));
					const poly T = rec(V, 2 * d, n / 2);
					poly T_(len / 2), UU(len);
					for (int i = 0; i < 2 * d; ++ i) copy_n(T.f.begin() + i * (n / 2 + 1), n / 2 + 1, T_.f.begin() + i * (n + 1));
					ntt(T_.data(), lg_len - 1);
					for (int i = 0; i < len; i += 2) UU[i] = 1ll * T_[i / 2] * Q_[i + 1] % mod, UU[i + 1] = 1ll * T_[i / 2] * Q_[i] % mod;
					intt(UU.data(), lg_len);
					poly U(d * (n + 1));
					for (int i = 0; i < d; ++ i) copy_n(UU.f.begin() + (i + d) * (n * 2 + 2), n + 1, U.f.begin() + i * (n + 1));
					return U;
				};
				return rec(Q, 1, max(n - 1, size() - 1)).resize(n);
			}
		};
		inline poly operator""_p(const char *str, size_t len) {
			poly ans(2);
			int sgn = 1, phase = 0, coeff = 0, touch = 0, cnum = 0;
			auto clean = [&]() {if(sgn==-1)coeff=(coeff==0?coeff:mod-coeff);if(phase==-1)ans[1]+=coeff;else if(phase==0){if(sgn==-1)cnum=(cnum==0?0:mod-cnum);ans[0]+=(int)cnum;}else if(phase==1)ans.resize(max(cnum+1,ans.size())),ans[cnum]+=coeff;else assert(0);phase=cnum=touch=0; };
			for (int i = 0; i < (int)len; ++i) {
				if (str[i] == '+') clean(), sgn = 1;
				else if (str[i] == '-') clean(), sgn = -1;
				else if ('0' <= str[i] and str[i] <= '9') {
					assert(phase == 0 || phase == 1);
					if (phase == 0) touch = 1, cnum = (10ll * cnum + str[i] - 48) % mod;
					else cnum = 10ll * cnum + str[i] - 48, assert(cnum < 1e8);
				} else if (str[i] == 'x') {
					while (str[i + 1] == ' ') ++i;
					assert(str[i + 1] == '^' || str[i + 1] == '+' || str[i + 1] == '-' || str[i + 1] == 0);
					phase = -1;
					coeff = touch ? cnum : 1;
					cnum = 0;
				} else if (str[i] == '^') {
					assert(phase == -1);
					phase = 1;
				}
			} clean();
			return ans;
		}
		namespace __semiconvol__ {
			const int logbr = 4, br = 1 << logbr, maxdep = (maxbit - 1) / logbr + 1, __bf = 7, bf = max(__bf, logbr - 1), pbf = 1 << bf;
			inline void src(poly& f, const poly& g, const function<void(const int& , poly& , const poly& )>& relax) {
				int nsize = g.size(), mxb = lg(nsize - 1);
				f.resize(1 << mxb);
				vu32 __prentt[maxdep][br];
				for (int i = 0, k = mxb; k > bf; k -= logbr, i++) {
					for (int j = 0; j < br - 1; j++) {
						if ((j << (k - logbr)) >= nsize)
							break;
						__prentt[i][j].resize(2 << (k - logbr));
						int nl = (j << (k - logbr)), nr = min(((j + 2) << (k - logbr)), nsize) - nl;
						memcpy(__prentt[i][j].data(), g.data() + nl, nr * 4);
						ntt(__prentt[i][j].data(), k - logbr + 1);
					}
				}
				function<void(int, int, int)> __div = [&](int x, int l, int r) {if(r-l<=pbf){for(int i=l;i<r;i++){relax(i,f,g);if(i+1<r)for(int j=i+1;j<r;j++)f[j]=(f[j]+1ll*f[i]*g[j-i])%mod;}return;}int nbit=mxb-logbr*(x+1),nbr=0;vu32 __tmp[br];while(l+(nbr<<nbit)<r){__tmp[nbr].resize(2<<nbit);nbr++;}for(int i=0;i<nbr;i++){if(i!=0){intt(__tmp[i].data(),nbit+1);for(int j=0;j<(1<<nbit);j++){u32&x=f[l+(i<<nbit)+j],&y=__tmp[i][j+(1<<nbit)];x=x+y-(x+y>=mod)*mod,y=0;}}__div(x+1,l+(i<<nbit),min(l+((i+1)<<nbit),r));if(i!=nbr-1){memcpy(__tmp[i].data(),f.data()+l+(i<<nbit),4<<nbit);ntt(__tmp[i].data(),nbit+1);for(int j=i+1;j<nbr;j++)for(int k=0;k<(2<<nbit);k++)__tmp[j][k]=(__tmp[j][k]+1ll*__tmp[i][k]*__prentt[x][j-i-1][k])%mod;}} };
				__div(0, 0, nsize);
				f.resize(nsize);
			}
		}
		using __semiconvol__::src;
		inline poly poly::ln() const {
			poly ret;
			src(ret, *this, [&](const int& i, poly& f, const poly& g) {if(i==0)f[i]=0;else f[i]=(1ll*g[i]*i+mod-f[i])%mod; });
			for (int i = degree(); i >= 1; -- i) ret[i] = 1ll * ginv(i) * ret[i] % mod;
			return ret;
		}
		inline poly poly::exp() const {
			poly ret, tmp(*this);
			for (int i = 0; i < size(); ++ i) tmp[i] = 1ll * tmp[i] * i % mod;
			src(ret, tmp, [&](const int& i, poly& f, const poly& g) {if(i==0)f[i]=1;else f[i]=1ll*f[i]*ginv(i)%mod; });
			return ret;
		}
		inline poly poly::inv() const {
			poly ret, tmp(*this);
			int ivf0 = qp(f[0], mod - 2);
			tmp[0] = 0;
			src(ret, tmp, [&](const int& i, poly& f, const poly& g) {if(i==0)f[i]=ivf0;else f[i]=1ll*ivf0*(mod-f[i])%mod; });
			return ret;
		}
		inline poly poly::quo(const poly& g) const {
			using namespace __semiconvol__;
			int nsize = f.size(), mxb = lg(nsize - 1);
			vu32 res(1 << mxb), __prentt[maxdep][br], _f(g.f);
			u32 ivf0 = qp(_f[0], -1);
			_f[0] = 0, _f.resize(nsize);
			for (int i = 0, k = mxb; k > bf; k -= logbr, i++) {
				for (int j = 0; j < br - 1; j++) {
					if ((j << (k - logbr)) >= nsize) break;
					__prentt[i][j].resize(2 << (k - logbr));
					int nl = (j << (k - logbr)), nr = min(((j + 2) << (k - logbr)), nsize) - nl;
					memcpy(__prentt[i][j].data(), _f.data() + nl, nr * 4);
					ntt(__prentt[i][j].data(), k - logbr + 1);
				}
			}
			function<void(int, int, int)> __div = [=,& res,& __prentt,& _f,& mxb,& __div,& ivf0](int x, int l, int r) {if(r-l<=pbf){for(int i=l;i<r;i++){res[i]=1ll*ivf0*(i==0?f[0]:f[i]+mod-res[i])%mod;if(i+1<r){u64 __tmp=res[i];for(int j=i+1;j<r;j++)res[j]=(res[j]+__tmp*_f[j-i])%mod;}}return;}int nbit=mxb-logbr*(x+1),nbr=0;vu32 __tmp[br];while(l+(nbr<<nbit)<r){__tmp[nbr].resize(2<<nbit);nbr++;}for(int i=0;i<nbr;i++){if(i!=0){intt(__tmp[i].data(),nbit+1);for(int j=0;j<(1<<nbit);j++){u32&x=res[l+(i<<nbit)+j],&y=__tmp[i][j+(1<<nbit)];x=x+y-(x+y>=mod)*mod,y=0;}}__div(x+1,l+(i<<nbit),min(l+((i+1)<<nbit),r));if(i!=nbr-1){memcpy(__tmp[i].data(),res.data()+l+(i<<nbit),4<<nbit);ntt(__tmp[i].data(),nbit+1);for(int j=i+1;j<nbr;j++)for(int k=0;k<(2<<nbit);k++)__tmp[j][k]=(__tmp[j][k]+1ll*__tmp[i][k]*__prentt[x][j-i-1][k])%mod;}} };
			__div(0, 0, nsize);
			return res.resize(nsize), res;
		}
		namespace __multipoint_operation__ {
			vector<poly> __Q;
			poly _E_Mul(poly A, poly B) {
				int n = A.size(), m = B.size();
				B.rev(), B = A * B;
				for (int i = 0; i < n; ++i) A[i] = B[i + m - 1];
				return A;
			}
			void _E_Init(int p, int l, int r, poly& a) {
				if (l == r) {
					__Q[p].resize(2);
					__Q[p][0] = 1, __Q[p][1] = (a[l] ? mod - a[l] : a[l]);
					return;
				} int mid = l + r >> 1;
				_E_Init(p << 1, l, mid, a), _E_Init(p << 1 | 1, mid + 1, r, a);
				__Q[p] = __Q[p << 1] * __Q[p << 1 | 1];
			}
			void _E_Calc(int p, int l, int r, const poly& F, poly& g) {
				if (l == r) return void(g[l] = F[0]);
				poly __F(r - l + 1);
				for (int i = 0, ed = r - l + 1; i < ed; ++i) __F[i] = F[i];
				int mid = l + r >> 1;
				_E_Calc(p << 1, l, mid, _E_Mul(__F, __Q[p << 1 | 1]), g);
				_E_Calc(p << 1 | 1, mid + 1, r, _E_Mul(__F, __Q[p << 1]), g);
			}
			vector<poly> __P;
			void _I_Init(int p, int l, int r, const poly& x) {
				if (l == r) {
					__P[p].resize(2), __P[p][0] = (x[l] ? mod - x[l] : 0), __P[p][1] = 1;
					return;
				} int mid = l + r >> 1;
				_I_Init(p << 1, l, mid, x), _I_Init(p << 1 | 1, mid + 1, r, x);
				__P[p] = __P[p << 1] * __P[p << 1 | 1];
			}
			poly _I_Calc(int p, int l, int r, const poly& t) {
				if (l == r) return poly(1, t[l]);
				int mid = l + r >> 1;
				poly L(_I_Calc(p << 1, l, mid, t)), R(_I_Calc(p << 1 | 1, mid + 1, r, t));
				L = L * __P[p << 1 | 1], R = R * __P[p << 1];
				for (int i = 0; i < (int)R.size(); ++i) {
					L[i] = L[i] + R[i];
					if (L[i] >= mod) L[i] -= mod;
				} return L;
			}
		}
		inline poly poly::eval(int n, poly a) const {
			using namespace __multipoint_operation__;
			n = max(n, size());
			poly v(n), F(f);
			__Q.resize(n << 2);
			F.resize(n + 1), a.resize(n);
			_E_Init(1, 0, n - 1, a);
			__Q[1].resize(n + 1);
			_E_Calc(1, 0, n - 1, _E_Mul(F, __Q[1].inv()), v);
			return v;
		}
		inline poly poly::intp(int n, const poly& x, const poly& y) {
			using namespace __multipoint_operation__;
			__P.resize(n << 2);
			_I_Init(1, 0, n - 1, x);
			__P[1] = __P[1].deri();
			poly t = __P[1].eval(n, x);
			for (int i = 0; i < n; ++i) t[i] = 1ll * y[i] * qp(t[i], mod - 2) % mod;
			f = _I_Calc(1, 0, n - 1, t);
			return *this;
		}
	} using namespace polynomial;
    	namespace plugins {
		poly _S1R_solve(const int& n) {
			if (n == 0) return poly(1, 1);
			int mid = n >> 1;
			poly f = _S1R_solve(mid);
			f.resize(mid + 1);
			poly A = f.shift(mid), B(mid + 1);
			for (int i = 0; i <= mid; ++i) B[i] = f[i];
			A *= B, f.resize(n + 1);
			A.redegree(n);
			if (n&  1) for (int i = 0; i <= n; ++i)
				f[i] = ((i ? A[i - 1] : 0) + 1ll * (n - 1) * A[i]) % mod;
			else for (int i = 0; i <= n; ++i)
				f[i] = A[i];
			return f;
		}
		inline poly stirling2_row(const int& n) {
			poly A(n + 1), B(n + 1);
			for (int i = 0; i <= n; i++)
				A[i] = (i&  1 ? mod - gifc(i) : gifc(i)), B[i] = 1ll * qp(i, n) * gifc(i) % mod;
			return (A * B).split(n);
		}
		inline poly stirling2_col(const int& n, const int& k) {
			poly f(n + 1);
			for(int i = 1; i <= n; ++ i) f[i] = gifc(i);
			f = f.pow(k);
			for(int i = 1; i <= n; ++ i) f[i] = 1ll * f[i] * gfac(i) % mod * gifc(k) % mod;
			return f;
		}
		inline poly stirling1_row(const int& n) { return _S1R_solve(n); }
		inline poly stirling1_col(const int& n, const int& k) {
			poly f(n + 1);
			f[0] = 1, f[1] = mod - 1;
			f = f.inv().ln().pow(k);
			for(int i = 1; i <= n; ++ i) f[i] = 1ll * f[i] * gfac(i) % mod * gifc(k) % mod;
			return f;
		}
		inline poly eulerian_col(const int& n, const int& k) {
			poly f1(max(k, n) + 1), f2(max(k, n) + 1), xex(n + 1), ekx(n + 1), ek1x(n + 1);
			for(int i = 0; i < n; ++ i) {
				if (i&  1) xex[i + 1] = mod - gifc(i);
				else xex[i + 1] = gifc(i);
			}
			int _c1 = 1, _c2 = 1;
			for(int i = 0; i <= n; ++ i) {
				ekx[i] = 1ll * _c1 * gifc(i) % mod;
				ek1x[i] = 1ll * _c2 * gifc(i) % mod;
				_c1 = 1ll * _c1 * k % mod;
				_c2 = 1ll * _c2 * (k + 1) % mod;
			} 
			for(int i = 0; i <= k; ++ i) f1[i] = 1ll * qp(mod + i - k - 1, i) * gifc(i) % mod;
			for(int i = 0; i < k; ++ i) f2[i] = 1ll * qp(mod - k + i, i) * gifc(i) % mod;
			f1 = ek1x * f1.composite(xex);
			f2 = ekx * f2.composite(xex);
			f1 -= f2; f1.resize(n); 
			for(int i = 0; i < n; ++ i) f1[i] = 1ll * f1[i] * gfac(i) % mod;
			return f1;
		}
		inline poly SEQ(const poly& A) { return (1 - A).inv(); }
		inline poly MSET(const poly& A) {
			poly ret(A);
			for (int i = 2; i < A.size(); ++i)
				for (int j = 1; i * j < A.size(); ++j)
					ret[i * j] = (ret[i * j] + 1ll * ginv(i) * A[j]) % mod;
			return ret.exp();
		}
		inline poly Exp(const poly& A) { return MSET(A); }
		inline poly Ln(poly f) {
			f[0] = 1, f = f.ln();
			for (int i = 1; i < f.size(); ++i) f[i] = 1ll * f[i] * i % mod;
			poly prime(f.size() + 1), mu(f.size() + 1), ret(f.size()); mu[1] = 1;
			vector<bool> vis(f.size() + 1); int cnt = 0; 
			for (int i = 2; i <= f.size(); ++i) {
				if (!vis[i]) prime[++cnt] = i, mu[i] = mod - 1;
				for (int j = 1; j <= cnt and i * prime[j] <= f.size(); ++j) {
					vis[i * prime[j]] = 1;
					if (i % prime[j] == 0) break;
					mu[i * prime[j]] = (mu[i] == 0 ? 0 : mod - mu[i]);
				}
			}
			for (int i = 1; i < f.size(); ++i)
				for (int j = 1; i * j < f.size(); ++j)
					ret[i * j] = (ret[i * j] + 1ll * f[i] * mu[j]) % mod;
			for (int i = 1; i < ret.size(); ++i) if (ret[i]) ret[i] = 1ll * ret[i] * ginv(i) % mod;
			ret[0] = 0;
			return ret;
		}
		inline int _R_Div(poly F, poly G, u64 n) {
			int len = max(F.size(), G.size());
			int m = 1, o = 0;
			while (m < len) m <<= 1, ++o;
			F.resize(1 << o + 1), G.resize(1 << o + 1);
			while (n > m) {
				ntt(F.data(), o + 1), ntt(G.data(), o + 1);
				for (int i = 0; i < m * 2; ++i) F[i] = 1ll * F[i] * G[i ^ 1] % mod;
				for (int i = 0; i < m; ++i) G[i] = 1ll * G[i << 1] * G[i << 1 | 1] % mod;
				intt(F.data(), o + 1), intt(G.data(), o);
				for (int i = 0, j = n&  1; i < len; i++, j += 2) F[i] = F[j];
				for (int i = len, ed = 1 << o + 1; i < ed; ++i) F[i] = G[i] = 0;
				n >>= 1;
			} G.resize(m), G = G.inv();
			int s = n; n = F.size() - 1, m = G.size() - 1;
			int a = max(0, s - m), b = min(s, (int) n), ans = 0;
			for (int i = a; i <= b; ++i)
				ans = (ans + 1ll * F[i] * G[s - i]) % mod;
			return ans;
		}
		inline int linear_recur(u64 n, int k, poly f, poly a) {
			poly F(k + 1); F[k] = 1; assert(a.size() >= k); a.resize(k);
			for (int i = 1; i <= k; ++i) F[k - i] = (f[i] == 0 ? 0 : mod - f[i]);
			F.rev(); f = (a * F).split(a.degree());
			return _R_Div(f, F, n);
		}
		inline poly BM(int n, poly a) {
			a.f.emplace(a.f.begin(), 0);
			poly ans, lst; int w = 0; u64 delt = 0;
			for (int i = 1; i <= n; ++i) {
				u64 tmp = 0;
				for (int j = 0; j < ans.size(); ++j) tmp = (tmp + 1ll * a[i - j - 1] * ans[j]) % mod;
				if ((a[i] - tmp + mod) % mod == 0) continue;
				if (!w) {
					w = i, delt = Norm(a[i] - tmp + mod);
					for (int j = i; j; --j) ans.f.emplace_back(0);
					continue;
				} poly now = ans;
				u64 mult = 1ll * (a[i] - tmp + mod) * qp(delt, mod - 2) % mod;
				if (ans.size() < lst.size() + i - w) ans.resize(lst.size() + i - w);
				ans[i - w - 1] += mult;
				if (ans[i - w - 1] >= mod) ans[i - w - 1] -= mod;
				for (int j = 0; j < lst.size(); ++j) ans[i - w + j] = (ans[i - w + j] - 1ll * mult * lst[j] % mod + mod) % mod;
				if (now.size() - i < lst.size() - w) lst = now, w = i, delt = (a[i] - tmp + mod) % mod;
			} return ans.f.emplace(ans.f.begin(), 0), ans;
		}
		inline int recur_by_bm(u64 n, int k, poly a) {
			poly f = BM(k, a);
			if (f.size() == 2)
				return 1ll * a[0] * qp(f[1], n - 1) % mod;
			return linear_recur(n, f.degree(), f, a);
		}

		namespace mod_poly {
			struct FastMod {
				int m; ll b;
				inline void init(int _m = 1) { m = _m; if (m == 0) m = 1; b = ((lll)1 << 64) / m; } 
				FastMod(int _m = 1) { init(_m); }
				inline int operator()(ll a) { ll q = ((lll)a * b) >> 64; a -= q * m; if (a >= m) a -= m; return a; }
			} Mod(mod);
			struct Z {
				u32 v;
				Z(u32 v = 0) : v(v) {}
				Z(int v) : v(Norm(Mod(v) + mod)) { }
				Z(long long v) : v(Norm(Mod(v) + mod)) { }
				inline friend Z operator+(const Z& lhs, const Z &rhs) { return Norm(lhs.v + rhs.v); }
				inline friend Z operator+(const Z& lhs, const int &rhs) { return Norm(lhs.v + rhs); }
				inline friend Z operator+(const int& lhs, const Z& rhs) { return Norm(lhs + rhs.v); }
				inline friend Z operator-(const Z& lhs, const Z &rhs) { return Norm(lhs.v + mod - rhs.v); }
				inline friend Z operator-(const Z& lhs, const int &rhs) { return Norm(lhs.v - rhs + mod); }
				inline friend Z operator-(const int& lhs, const Z& rhs) { return Norm(lhs - rhs.v + mod); }
				Z operator-() const { return Norm(mod - v); }
				inline Z inv() const { return qp(*this, mod - 2); };
				inline friend Z operator*(const Z& lhs, const Z& rhs) { return Mod(1ll * lhs.v * rhs.v); }
				inline friend Z operator*(const Z& lhs, const int &rhs) { return Mod(1ll * lhs.v * rhs); }
				inline friend Z operator*(const int& lhs, const Z& rhs) { return Mod(1ll * lhs * rhs.v); }
				inline friend Z operator/(const Z& lhs, const Z &rhs) { return Mod(1ll * lhs.v * rhs.inv().v); }
				inline friend Z operator/(const Z& lhs, const int &rhs) { return Mod(1ll * lhs.v * qp(rhs, mod - 2)); }
				inline friend Z operator/(const int& lhs, const Z& rhs) { return Mod(1ll * lhs * rhs.inv().v); }
				operator u32() const { return v; }
				inline friend ostream &operator<< (ostream &out, const Z &x) { return out << x.v; }
			};
			inline Z &operator+=(Z &lhs, const Z &rhs) { return lhs = lhs + rhs; }
			inline Z &operator-=(Z &lhs, const Z &rhs) { return lhs = lhs - rhs; }
			inline Z &operator*=(Z &lhs, const Z &rhs) { return lhs = lhs * rhs; }
			inline Z &operator/=(Z &lhs, const Z &rhs) { return lhs = lhs / rhs; }

			struct opoly : vector<Z> {
				opoly(const int& n = 1, const Z& val = 0) : vector(n, val) {}
				opoly(const initializer_list<value_type> &il) : vector(il) {}
				opoly(const vector<Z> &il) : vector(il) {}
				inline int degree() const { return (int) size() - 1; }
				inline bool shrink() { int k = size(); while (k && !at(k - 1)) --k; resize(k); return k; }
				inline void redegree(const size_t& n) { resize(n + 1); }
				inline void rev() { reverse(begin(), end()); }
				inline opoly operator-() const { opoly ret(size()); for (int i = 0; i < size(); ++i) ret[i] = -at(i); return ret; }
				inline operator poly() const { poly ret(size()); for (int i = 0; i < (int)size(); ++ i) ret[i] = u32(at(i)); return ret; }
				inline friend opoly operator*(const opoly &a, const opoly &b) {
					int n = a.degree(), m = b.degree();
					if (n == -1 || m == -1) return opoly();
					opoly c(n + m + 1);
					for (int i = 0; i <= n + m; ++i)
						for (int j = max(0, i - m); j <= min(i, n); ++j)
							c[i] += a[j] * b[i - j];
					return c;
				}
				inline friend opoly operator*(const opoly &a, const Z &z) { opoly c(a); for (Z &x : c) x *= z; return c; }
				inline friend opoly operator*(const Z &z, const opoly &a) { opoly c(a); for (Z &x : c) x *= z; return c; }
				inline friend opoly operator+(const opoly &a, const opoly &b) {
					opoly c(max(a.size(), b.size()));
					for (int i = 0; i < a.size(); ++i) c[i] += a[i];
					for (int i = 0; i < b.size(); ++i) c[i] += b[i];
					return c;
				}
				inline friend opoly operator-(const opoly &a, const opoly &b) { return a + -b; }
				inline friend opoly operator+(const opoly &a, const Z &z) { opoly c(a); c[0] += z; return c; }
				inline friend opoly operator+(const Z &z, const opoly &a) { opoly c(a); c[0] += z; return c; }
				inline friend opoly operator-(const opoly &a, const Z &z) { opoly c(a); c[0] -= z; return c; }
				inline friend opoly operator-(const Z &z, const opoly &a) { opoly c(a); c[0] -= z; return c; }
				inline friend bool operator==(const opoly &a, const opoly &b) { 
					if (a.size() != b.size()) return false; 
					for (int i = 0; i < a.size(); ++ i)  if (a[i] != b[i]) return false;
					return true;
				}
				opoly deri() const {
					if (empty()) return *this;
					opoly a(*this);
					for (int i = 1; i < a.size(); ++i) a[i - 1] = a[i] * i;
					a.pop_back();
					return a;
				}
				Z eval(const Z &z) const {
					Z v = 0;
					for (int i = degree(); i >= 0; --i) v = v * z + at(i);
					return v;
				}
				poly trans_poly() {
					poly ret(size());
					for (int i = 0; i < size(); ++ i) ret[i] = u32(at(i));
					return ret;
				}
			};
			template <typename _vi>
			inline opoly topoly(const _vi& a) {
				opoly ret(a.size());
				for (int i = 0; i < a.size(); ++ i) 
					ret[i].v = a[i];
				return ret;
			}
			inline opoly gcd(opoly a, opoly b) {
				if (!a.shrink()) return b;
				if (!b.shrink()) return a;
				if (a.size() < b.size()) swap(a, b);
				while (b.shrink()) {
					Z in = b.back().inv();
					for (Z &x : b) x *= in;
					int n = a.degree(), m = b.degree();
					for (int i = n; i >= m; --i) {
						for (int j = 1; j <= m; ++j)
							a[i - j] -= a[i] * b[m - j];
						a[i] = 0;
					} swap(a, b);
				} return a.shrink(), a;
			}
			inline opoly div(const opoly& _a, const opoly& _b) {
				opoly a(_a), b(_b);
				Z in = b.back().inv();
				for (Z &x : b) x *= in;
				int n = a.degree(), m = b.degree();
				opoly ret(n - m + 1);
				for (int i = n; i >= m; --i) {
					ret[i - m] = a[i] * b[m];
					for (int j = 1; j <= m; ++j)
						a[i - j] -= a[i] * b[m - j];
				} for (Z &x : ret) x = x * in;
				return ret;
			}
			inline opoly operator/ (const opoly& _a, const opoly& _b) { return div(_a, _b); }
			inline opoly& operator/= (opoly& _a, const opoly& _b) { _a = div(_a, _b); return _a; }
			ostream &operator<<(ostream &out, const opoly &p) {
				if (p.empty()) return out << 0;
				out << p[0]; for (int i = 1; i < p.size(); ++i)
					out << " + " << p[i] << "x^" << i;
				return out;
			}
			template<class T> istream &operator>>(istream &is, vector<T> &v) { for (T &x : v) is >> x; return is; }
		}
		
		namespace _PR_base {
			using namespace mod_poly;
			int order, mxdeg;

			vector<opoly> P;

			opoly lagrange(opoly a0, const Z& v) {
				int n = a0.size() - 1;
				assert(!(v.v <= n)); assert(!((v + n).v <= n));
				for (int i = 0; i <= n; i++) {
					a0[i] = a0[i] * Z(gifc(i)) * Z(gifc(n - i));
					if ((n - i) & 1) a0[i] = - a0[i];
				} opoly prf(2 * n + 2), ivp(2 * n + 2);
				prf[0] = 1;
				for (int i = 0; i <= 2 * n; i++) prf[i + 1] = prf[i] * (v - n + i);
				ivp[2 * n + 1] = qp(prf[2 * n + 1], mod - 2);
				for (int i = 2 * n; i >= 0; i--) ivp[i] = ivp[i + 1] * (v - n + i);
				opoly f(2 * n + 1), ret(n + 1);
				for (int i = 0; i <= 2 * n; i++) f[i] = ivp[i + 1] * prf[i];
				f = topoly(f.trans_poly() * a0.trans_poly());
				for (int i = 0; i <= n; i++) ret[i] = f[i + n];
				for (int i = 0; i <= n; i++) ret[i] = ret[i] * ivp[i] * prf[i + n + 1];
				return ret;
			}

			vector<vector<opoly>> B, tmpM;
			opoly B2, tmpZ;
			int calc(const u64& n, opoly a) {
				int s = (int)ceil(sqrt(1. * (n - order) / mxdeg)); s = 1 << lg(s);
				// while (order + s * s * mxdeg <= n) s <<= 1;
				a.rev();
				for (int i = 0; i < order; ++ i) for (int j = 0; j < order; ++ j) B[i][j].redegree(mxdeg);
				poly tmpv(mxdeg + 1, 0);
				for (int i = 0; i <= mxdeg; ++ i) tmpv[i] = - P[0].eval(order + i);
				for (int i = 0; i < order; ++ i) for (int j = 0; j <= mxdeg; ++ j) {
					B[i][0][j] = P[i + 1].eval(order + j);
					if (i > 0) B[i - 1][i][j] = tmpv[j];
				}
				tmpM.resize(order); for (int i = 0; i < order; ++ i) tmpM[i].resize(order);
				for (int t = 2; t <= s; t <<= 1) {
					for (int i = 0; i < order; ++ i) for (int j = 0; j < order; ++ j) {
						auto tmp = lagrange(B[i][j], Z(mxdeg) * Z(t >> 1) + 1);
						B[i][j].insert(B[i][j].end(), tmp.begin(), tmp.end());
						tmp = lagrange(B[i][j], Z(mxdeg) * Z(t) + 2);
						B[i][j].insert(B[i][j].end(), tmp.begin(), tmp.end());
						tmpM[i][j].clear(); tmpM[i][j].resize(mxdeg * t + 1);
					}
					for (int i = 0; i < order; ++ i) for (int j = 0; j < order; ++ j)
						for (int k = 0; k < order; ++ k) for (int v = 0; v <= mxdeg * t; ++ v) 
							tmpM[i][k][v] += B[i][j][v << 1] * B[j][k][v << 1 | 1];
					for (int i = 0; i < order; ++ i) for (int j = 0; j < order; ++ j) {
						B[i][j].resize(tmpM[i][j].size());
						for (int k = 0; k < B[i][j].size(); ++ k) B[i][j][k] = tmpM[i][j][k]; 
					}
				} 
				for (int i = 0; i < (n - order) / s; ++ i) {
					opoly tmp(order, 0);
					for (int j = 0; j < order; ++ j) for (int k = 0; k < order; ++ k)
						tmp[k] += B[j][k][i] * a[j];
					for (int j = 0; j < order; ++ j) a[j] = tmp[j];
				} 
				for (int i = (n - order) / s * s + order; i <= n; ++ i) {
					Z tmp = 0;
					for (int j = 0; j < order; ++ j) tmp += P[j + 1].eval(i) * a[j];
					Z tmp2 = - P[0].eval(i);
					for (int j = order - 1; j; -- j) a[j] = a[j - 1] * tmp2;
					a[0] = tmp;
				}

				B2.redegree(mxdeg);
				for (int i = 0; i <= mxdeg; ++ i) B2[i] = - P[0].eval(order + i);
				for (int t = 2; t <= s; t <<= 1) {
					auto tmp = lagrange(B2, Z(mxdeg) * Z(t >> 1) + 1);
					B2.insert(B2.end(), tmp.begin(), tmp.end());
					tmp = lagrange(B2, Z(mxdeg) * Z(t) + 2);
					B2.insert(B2.end(), tmp.begin(), tmp.end());
					tmpZ.resize(0); tmpZ.resize(mxdeg * t + 1);
					for (int v = 0; v <= mxdeg * t; ++ v)
						tmpZ[v] = B2[v << 1] * B2[v << 1 | 1];
					B2 = tmpZ;
				} 
				Z ans = 1;
				for (int i = 0; i < (n - order) / s; ++ i) ans = ans * B2[i];
				for (int i = (n - order) / s * s + order; i <= n; ++ i) ans = ans * (- P[0].eval(i));
				return (a[0] * qp(ans, mod - 2)).v;
			}

			inline int p_recur(const u64& n, const int& m, const int& d, const poly& a, const poly* P) {
				/* m 阶整式递推，max deg P_i = d，满足 \sum_{k = 0}^m a_{n - k} P_k(n) = 0，初值给 m 个：a_0, ..., a_{m - 1}，要求算 a_n */
				if (n <= a.degree()) return a[n];
				if (d == 0) {
					poly f(m + 1); 
					for (int i = 0; i <= m; ++ i) f[i] = P[i][0];
					int ivv = qp(mod - f[0], mod - 2);
					for (int i = 1; i <= m; ++ i) f[i] = 1ll * f[i] * ivv % mod;
					return linear_recur(n, m, f, a);
				}
				_PR_base :: order = m, _PR_base :: mxdeg = d, _PR_base :: P.resize(m + 1), _PR_base :: B.resize(m + 1);
				for (int i = 0; i <= m; ++ i) _PR_base :: B[i].resize(m + 1);
				for (int i = 0; i <= m; ++ i) _PR_base :: P[i] = topoly(P[i]), _PR_base :: P[i].shrink();
				return _PR_base :: calc(n, topoly(a.split(m - 1)));
			}
			inline int p_recur(const u64& n, const int& m, const poly& a, const poly* P) {
				int d = 0;
				for (int i = 0; i <= m; ++ i) d = max(d, P[i].size() - 1);
				return p_recur(n, m, d, a, P);
			}
			inline int p_recur(const u64& n, const poly& a, const vector<poly>& P) {
				int m = (int)P.size() - 1, d = 0;
				for (int i = 0; i <= m; ++ i) d = max(d, P[i].size() - 1);
				return p_recur(n, m, d, a, P.data());
			}
		} using _PR_base :: p_recur;
		
		namespace ODE_Automaton {
		using namespace mod_poly;

		struct pfrac {
			opoly x, y;
			pfrac(const opoly &x = opoly(), const opoly &y = {Z(1)}) : x(x), y(y) {}
			void shrink() {
				y.shrink();
				if (!x.shrink()) { y = opoly{Z(1)}; return; }
				opoly g = gcd(x, y);
				x /= g, y /= g;
			}
			inline pfrac operator+(const pfrac &rhs) const { pfrac ret = pfrac(x * rhs.y + y * rhs.x, y * rhs.y); ret.shrink(); return ret; }
			inline pfrac operator-() const { pfrac ret = pfrac(-x, y); ret.shrink(); return ret; }
			inline pfrac operator-(const pfrac &rhs) const { pfrac ret = *this + -rhs; ret.shrink(); return ret; }
			inline pfrac operator*(const pfrac &rhs) const { pfrac ret = pfrac(x * rhs.x, y * rhs.y); ret.shrink(); return ret; }
			// inline pfrac operator*(const Z &rhs) const { pfrac ret = pfrac(x * rhs, y * rhs); ret.shrink(); return ret; } // ?
			inline pfrac inv() const { pfrac ret = pfrac(y, x); ret.shrink(); return ret; }
			inline pfrac operator/(const pfrac &rhs) const { pfrac ret = *this * rhs.inv(); ret.shrink(); return ret; }
			bool operator==(const pfrac &rhs) const { return x * rhs.y == y * rhs.x; }
			bool operator!=(const pfrac &rhs) const { return !operator==(rhs); }
			pfrac deri() const { pfrac ret = pfrac(x.deri() * y - y.deri() * x, y * y); ret.shrink(); return ret; }
		};

		struct Q_Basis {
			int dim, id;
			vector<vector<pfrac> > basis, augment;
			Q_Basis(int dim) : dim(dim), id(), basis(dim), augment(dim) {}
			vector<pfrac> insert(vector<pfrac> vec) {
				vector<pfrac> tmp(dim + 1);
				tmp[id++] = pfrac({Z(1)});
				for (int i = 0; i < dim; ++i) {
					if (vec[i] != pfrac()) {
						if (basis[i].empty()) {
							for (int j = i + 1; j < dim; ++j) {
								vec[j] = vec[j] / vec[i];
								vec[j].shrink();
							}
							for (int j = 0; j < id; ++j) {
								tmp[j] = tmp[j] / vec[i];
								tmp[j].shrink();
							}
							vec[i] = pfrac({Z(1)});
							basis[i] = vec;
							augment[i] = tmp;
							return {};
						} else {
							for (int j = i + 1; j < dim; ++j)
								vec[j] = vec[j] - vec[i] * basis[i][j];
							for (int j = 0; j < id; ++j)
								tmp[j] = tmp[j] - vec[i] * augment[i][j];
							vec[i] = pfrac();
						}
					}
				} return tmp;
			}
		};

		using PRec = vector<opoly>;
		using ODE_base = vector<opoly>;

		struct ODE {
			ODE_base ode;
			ODE(const int& n = 1, const opoly& val = {0}) : ode(n, val) {}
			ODE(const initializer_list<opoly> &il) : ode(il) {}
			inline int size() const { return (int)ode.size(); }
			inline int degree() const { return size() - 1; }
			inline opoly& operator[](int x) { return ode[x]; }
			inline opoly operator[](int x) const { return ode[x]; }
			inline ODE& resize(int x) { return ode.resize(x), *this; }
			inline ODE& redegree(int x) { return ode.resize(x + 1), *this; }
			inline ODE_base::iterator begin() { return ode.begin(); }
			inline ODE_base::iterator end() { return ode.end(); }
			void shrink() {
				opoly g = ode[0];
				for (opoly &x : ode) x.shrink(), g = gcd(g, x);
				for (opoly &x : ode) if (!x.empty()) x = div(x, g);
			}
			ODE deri() const {
				ODE _ode(*this);
				_ode.resize(_ode.size() + 1);
				for (int i = (int) _ode.size() - 2; i >= 0; --i) {
					_ode[i + 1] = _ode[i + 1] + _ode[i];
					_ode[i] = _ode[i].deri();
				} return _ode;
			}
			ODE theta() const {
				ODE _ode(deri());
				for (int i = 0; i < _ode.size(); ++i) _ode[i] = _ode[i] * opoly{Z(), Z(1)};
				return _ode;
			}

			PRec _rec;
			vector<Z> coef;
			inline PRec getPRec() {
				shrink();
				int tmp = numeric_limits<int>::max();
				int n = degree(), m = numeric_limits<int>::min();
				for (int i = 0; i <= n; ++i) {
					if (ode[i].empty()) continue;
					m = max(m, (int) ode[i].size() - 1 - i);
					int j = 0;
					while (ode[i][j] == 0) ++j;
					tmp = min(tmp, j - i);
				} m -= tmp;
				PRec rec(m + 1, opoly(n + 1)); opoly fall{Z(1)};
				for (int i = 0; i <= n; ++i) {
					opoly coef = fall;
					for (int j = 0; j < (int) ode[i].size() - i - tmp; ++j) {
						if (j + i + tmp >= 0) rec[j] = rec[j] + coef * ode[i][j + i + tmp];
						coef = div(coef * opoly{-Z(i + j), Z(1)}, opoly{-Z(j), Z(1)});
					} fall = fall * opoly{-Z(i), Z(1)};
				} 
				for (int i = 0; i <= m; ++i) rec[i].shrink();
				return _rec = rec;
			}
			int prf(int n) {
				coef.resize(n + 1);
				for (int i = 0; i <= n; ++i)
					coef[i] = _rec[0].eval(i);
				int r = 0;
				for (int i = n; i; --i)
					if (coef[i] == 0) {
						r = i;
						break;
					}
				return r;
			}
			opoly post(opoly init) {
				int m = init.size();
				auto invs = [&](const opoly &vec) {
					opoly prf(vec.size()), ret(vec.size());
					prf[0] = 1;
					for (int i = 1; i < vec.size(); ++i) prf[i] = prf[i - 1] * (vec[i - 1] == 0 ? Z(1) : vec[i - 1]);
					// Z tot = accumulate(vec.begin(), vec.end(), Z(1), multiplies<Z>()).inv();
					// for (int i = (int) vec.size() - 1; i >= 0; --i) ret[i] = tot * prf[i], tot *= vec[i];
                    Z tot = Z(1); 
                    for (auto v : vec) if (v.v != 0) tot *= v;
                    tot = tot.inv();
                    for (int i = (int) vec.size() - 1; i >= 0; --i) 
                        if (vec[i].v != 0) ret[i] = tot * prf[i], tot *= vec[i];
					return ret;
				};
				auto nvs = invs(vector<Z>(coef.begin() + m, coef.end()));
				init.resize(coef.size());
				for (int i = m; i < coef.size(); ++i) {
					for (int j = 1; j < min(i + 1, (int) _rec.size()); ++j)
						init[i] += init[i - j] * _rec[j].eval(i);
					init[i] = init[i] * -nvs[i - m];
				} return init;
			}
			opoly recur(const int& n, const opoly& a0) {
				getPRec();
				prf(n);
				opoly ret(a0.size());
				for (int i = 0; i < a0.size(); ++ i) ret[i].v = a0[i].v;
				if (n <= ret.degree()) return ret.redegree(n), ret;
				return post(ret);
			}
			poly recur(const int& n, const poly& a0) {
				getPRec();
				prf(n);
				opoly ret(a0.size());
				for (int i = 0; i < a0.size(); ++ i) ret[i].v = a0[i];
				if (n <= ret.degree()) return ret.redegree(n), ret;
				return post(ret);
			}
			vector<poly> trans_poly() {
				vector<poly> ret(_rec.size(), poly());
				for (int i = 0; i < _rec.size(); ++ i) 
					_rec[i].shrink(), ret[i] = _rec[i];
				return ret;
			}

			inline friend ODE operator+(const ODE &op, const ODE &oq) {
				int n = op.degree(), m = oq.degree();
				Q_Basis basis(n + m);
				vector<pfrac> pd(n + 1), qd(m + 1);
				pd[0] = qd[0] = pfrac({Z(1)});
				for (int dim = 0; dim <= n + m; ++dim) {
					vector<pfrac> vec(n + m);
					copy(pd.begin(), pd.begin() + n, vec.begin());
					copy(qd.begin(), qd.begin() + m, vec.begin() + n);
					auto ret = basis.insert(vec);
					if (!ret.empty()) {
						ODE ode(dim + 1); opoly prod = {Z(1)};
						for (int i = 0; i < dim; ++i) prod = prod * ret[i].y;
						ode[dim] = prod;
						for (int i = 0; i < dim; ++i) ode[i] = ret[i].x * div(prod, ret[i].y);
						ode.shrink();
						return ode;
					} pd[n] = qd[m] = pfrac();
					for (int j = n - 1; j >= 0; --j) pd[j + 1] = pd[j + 1] + pd[j], pd[j] = pd[j].deri();
					for (int j = m - 1; j >= 0; --j) qd[j + 1] = qd[j + 1] + qd[j], qd[j] = qd[j].deri();
					for (int j = 0; j < n; ++j) pd[j] = pd[j] - pd[n] * op[j] / op[n], pd[j].shrink();
					for (int j = 0; j < m; ++j) qd[j] = qd[j] - qd[m] * oq[j] / oq[m], qd[j].shrink();
				} return assert(0), ODE();
			}

			inline friend ODE operator*(const ODE &op, const ODE &oq) {
				int n = op.size() - 1, m = oq.size() - 1;
				Q_Basis basis(n * m);
				vector<vector<pfrac> > p(n + 1, vector<pfrac>(m + 1));
				p[0][0] = pfrac({Z(1)});
				for (int dim = 0; dim <= n * m; ++dim) {
					vector<pfrac> vec(n * m);
					for (int i = 0; i < n; ++i)
						for (int j = 0; j < m; ++j)
							vec[i * m + j] = p[i][j];
					auto ret = basis.insert(vec);
					if (!ret.empty()) {
						ODE ode(dim + 1); opoly prod = {Z(1)};
						for (int i = 0; i < dim; ++i) prod = prod * ret[i].y;
						ode[dim] = prod;
						for (int i = 0; i < dim; ++i) ode[i] = ret[i].x * div(prod, ret[i].y);
						ode.shrink();
						return ode;
					}
					for (int i = 0; i < n; ++i) p[i][m] = pfrac();
					for (int j = 0; j < m; ++j) p[n][j] = pfrac();
					for (int i = n - 1; i >= 0; --i)
						for (int j = m - 1; j >= 0; --j) {
							p[i + 1][j] = p[i + 1][j] + p[i][j];
							p[i][j + 1] = p[i][j + 1] + p[i][j];
							p[i][j] = p[i][j].deri();
						}
					for (int i = 0; i < n; ++i)
						for (int j = 0; j < m; ++j) 
							p[i][j] = p[i][j] - p[n][j] * op[i] / op[n] - p[i][m] * oq[j] / oq[m], p[i][j].shrink();
				} return assert(0), ODE();
			}

			inline ODE composite(pfrac q) const {
				q.shrink();
				int n = ode.size() - 1;
				vector<vector<pfrac>> tri(n + 1, vector<pfrac>(n + 1));
				tri[0][0] = pfrac({Z(1)});
				pfrac d = q.deri();
				d.shrink();
				for (int i = 0; i < n; ++i) {
					for (int j = 0; j <= i; ++j) {
						tri[i + 1][j] = tri[i + 1][j] + tri[i][j].deri();
						tri[i + 1][j + 1] = tri[i][j] * d;
					} for (int j = 0; j <= i + 1; ++j) tri[i + 1][j].shrink();
				} vector<pfrac> vec(n + 1);
				for (int i = 0; i <= n; ++i) 
					for (int j = ode[i].degree(); j >= 0; --j) vec[i] = vec[i] * q + pfrac({ode[i][j]}), vec[i].shrink();
				for (int i = n; i >= 0; --i) {
					vec[i] = vec[i] / tri[i][i], vec[i].shrink();
					for (int j = 0; j < i; ++j)  vec[j] = vec[j] - vec[i] * tri[i][j];
				} opoly prod = opoly{Z(1)};
				for (int i = 0; i <= n; ++i) prod = prod * vec[i].y;
				ODE ret(n + 1);
				for (int i = 0; i <= n; ++i) ret[i] = vec[i].x * div(prod, vec[i].y);
				return ret.shrink(), ret;
			}
		};
		ODE linear_add(ODE a, ODE b) {
			if (a.size() < b.size()) swap(a, b);
			for (int i = 0; i < b.size(); ++i) a[i] = a[i] + b[i];
			return a;
		}
		ODE scalar_mul(ODE ode, const Z &z) {
			for (int i = 0; i < ode.size(); ++i)
				ode[i] = ode[i] * z;
			return ode;
		}

		ostream &operator<<(ostream &out, const pfrac &q) { return out << '(' << q.x << ") / (" << q.y << ')'; }
		ostream &operator<<(ostream &out, const ODE &ode) {
			out << "(" << ode[0] << ")F";
			for (int i = 1; i < ode.size(); ++i) {
				out << " + (" << ode[i] << ")F^{(" << i << ")}";
			} return out << " = 0";
		}
		template<class T> ostream &operator<<(ostream &out, const vector<T> &v) {
			if (!v.empty()) {
				out << v.front();
				for (int i = 1; i < v.size(); ++i) out << ' ' << v[i];
			} return out;
		}

		const ODE ODE_EXP = {{-Z(1)}, {Z(1)} };
		const ODE ODE_LN = {{Z()}, {-Z(1)}, {Z(1), -Z(1)} };
		ODE ode_power(const Z &k) { return ODE{{-Z(k)}, {Z(), Z(1)}}; }
		ODE ode_pfrac(const pfrac &a) { return ode_power(1).composite(a); }
		ODE ode_pFq(const vector<Z> &a, const vector<Z> &b) {
			ODE l = {{Z(1)}}, r = l;
			for (int x : a) l = linear_add(l.theta(), scalar_mul(l, x));
			for (int x : b) r = linear_add(r.theta(), scalar_mul(r, x - 1));
			ODE ret = linear_add(r.deri(), scalar_mul(l, - Z(1)));
			return ret.shrink(), ret;
		} 
		} using namespace ODE_Automaton;
	} using namespace plugins;


} using namespace __POLY__;


int main() {
    ODE myode1 = {
        {Z(), Z(), -Z(1)}, 
        {Z(), Z(), -Z(1)}, 
        {Z(), Z(), Z(4), -Z(3)},
        {Z(), Z(), Z(), Z(32), -Z(1)},
        {Z(), Z(), Z(), Z(), Z(38)},
        {Z(), Z(), Z(), Z(), Z(), Z(12)},
        {Z(), Z(), Z(), Z(), Z(), Z(), Z(1)}
    };
    ODE myode2 = {
        {Z(), Z(-1)},
        {Z(), Z(1)},
        {Z(), Z(), Z(3)},
        {Z(), Z(), Z(), Z(1)}
    };
    cout << (myode1 * myode2).getPRec() << '\n';

}