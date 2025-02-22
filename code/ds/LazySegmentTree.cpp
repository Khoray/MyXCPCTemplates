/*
https://www.luogu.com.cn/problem/P3373
区间修改、区间查询
*/
#include<bits/stdc++.h>
using namespace std;

int mod;

struct S {
    int sum;
    int sz;
};

struct F {
    int add, mul;
};

inline S op(S a, S b) {
    return {(a.sum + b.sum) % mod, (a.sz + b.sz) % mod};
}

inline S mapping(F f, S a) {
    return {(1ll * a.sum * f.mul % mod + 1ll * f.add * a.sz % mod) % mod, a.sz};
}

// 先操作 g 再操作 f
inline F comp(F f, F g) {
    return {(1ll * g.add * f.mul % mod + f.add) % mod, 1ll * f.mul * g.mul % mod};
}

inline S e() {
    return {0, 0};
}

inline F id() {
    return {0, 1};
}

template<class S,
         auto op,
         auto e,
         class F,
         auto mapping,
         auto composition,
         auto id>
struct LazySegmentTree {
    static_assert(std::is_convertible_v<decltype(op), std::function<S(S, S)>>,
                  "op must work as S(S, S)");
    static_assert(std::is_convertible_v<decltype(e), std::function<S()>>,
                  "e must work as S()");
    static_assert(
        std::is_convertible_v<decltype(mapping), std::function<S(F, S)>>,
        "mapping must work as S(F, S)");
    static_assert(
        std::is_convertible_v<decltype(composition), std::function<F(F, F)>>,
        "composition must work as F(F, F)");
    static_assert(std::is_convertible_v<decltype(id), std::function<F()>>,
                  "id must work as F()");

    int _n;
    vector<S> a;
    vector<F> tag;

    // 下标从 1 开始，v[1..n]
    LazySegmentTree(int n): a(n << 2, e()), tag(n << 2, id()), _n(n) {}

    // 传入的 v 是 [0..n] 的, build 时使用 v[1.._n], _n=n=v.size()-1
    LazySegmentTree(const vector<S> &v): _n((int) v.size() - 1), a(_n << 2, e()), tag(_n << 2, id()) {
        auto build = [&] (auto &&build, int u, int l, int r) -> void {
            if(l == r) {
                a[u] = v[l];
                return;
            }
            int mid = (l + r) / 2;
            build(build, u << 1, l, mid);
            build(build, u << 1 | 1, mid + 1, r);
            a[u] = op(a[u << 1], a[u << 1 | 1]);

        };
        build(build, 1, 1, _n);
    }

    void pushdown(int u) {
        a[u << 1] = mapping(tag[u], a[u << 1]);
        a[u << 1 | 1] = mapping(tag[u], a[u << 1 | 1]);
        tag[u << 1] = composition(tag[u], tag[u << 1]);
        tag[u << 1 | 1] = composition(tag[u], tag[u << 1 | 1]);
        tag[u] = id();
    }

    // return op(a[l], .., a[r])
    S prod(int L, int R) {
        auto dfs = [&] (auto&& dfs, int u, int l, int r, int L, int R) -> S {
            if(l > R || L > r) return e();
            if(l >= L && r <= R) return a[u];
            pushdown(u);
            int mid = (l + r) / 2;
            return op(
                dfs(dfs, u << 1, l, mid, L, R),
                dfs(dfs, u << 1 | 1, mid + 1, r, L, R)
            );
        };
        return dfs(dfs, 1, 1, _n, L, R);
    }

    // apply f to a[l] .. a[r]
    void apply(int L, int R, F f) {
        auto dfs = [&] (auto&& dfs, int u, int l, int r, int L, int R, F f) -> void {
            if(l > R || L > r) return;
            if(l >= L && r <= R) {
                a[u] = mapping(f, a[u]);
                tag[u] = composition(f, tag[u]);
                return;
            }
            pushdown(u);
            int mid = (l + r) / 2;
            dfs(dfs, u << 1, l, mid, L, R, f);
            dfs(dfs, u << 1 | 1, mid + 1, r, L, R, f);
            a[u] = op(a[u << 1], a[u << 1 | 1]);
        };
        return dfs(dfs, 1, 1, _n, L, R, f);
    }

    // set v[p] := x
    void set(int p, S x) {
        auto dfs = [&] (auto&& dfs, int u, int l, int r, int p, S x) -> void {
            if(l == r) {
                a[u] = x;
                return;
            }
            pushdown(u);
            int mid = (l + r) / 2;
            if(p <= mid) dfs(dfs, u << 1, l, mid, p, x);
            else dfs(dfs, u << 1 | 1, mid + 1, r, p, x);
            a[u] = op(a[u << 1], a[u << 1 | 1]);
        };
        return dfs(dfs, 1, 1, _n, p, x);
    }
    
    // return a[1] .. a[n]
    vector<S> get_all() {
        vector<S> ret(_n + 1);
        auto dfs = [&] (auto&& dfs, int u, int l, int r) -> void {
            if(l == r) {
                ret[l] = a[u];
                return;
            }
            pushdown(u);
            int mid = (l + r) / 2;
            dfs(dfs, u << 1, l, mid);
            dfs(dfs, u << 1 | 1, mid + 1, r);
        };
        dfs(dfs, 1, 1, _n);
        return ret;
    }

    // binary search a `L<=r<=_n`: check(op(a[L]..a[r-1]))=true, check(op(a[L]..a[r]))=false
    // if not exist: return -1
    template<typename Check>
    int max_right(int L, Check check) {
        auto dfs = [&] (auto&& dfs, int u, S &s, int l, int r, int L, Check check) -> int {
            if(r < L) return -1;
            if(l >= L) {
                S ss = op(s, a[u]);
                if(check(ss)) {
                    s = ss;
                    return -1;
                }
                if(l == r) return l;
            }
            pushdown(u);
            int mid = (l + r) / 2;
            int pos = dfs(dfs, u << 1, s, l, mid, L, check);
            if(pos != -1) return pos;
            return dfs(dfs, u << 1 | 1, s, mid + 1, r, L, check);
        };
        S temp_s = e();
        return dfs(dfs, 1, temp_s, 1, _n, L, check);
    }

    // binary search a `1<=l<=R`: check(op(a[l+1]..a[R]))=true, check(op(a[l]..a[R]))=false
    // if not exist: return -1
    template<typename Check>
    int min_left(int R, Check check) {
        auto dfs = [&] (auto&& dfs, int u, S &s, int l, int r, int R, Check check) -> int {
            if(l > R) return -1;
            if(r <= R) {
                S ss = op(a[u], s);
                if(check(ss)) {
                    s = ss;
                    return -1;
                }
                if(l == r) return l;
            }
            pushdown(u);
            int mid = (l + r) / 2;
            int pos = dfs(dfs, u << 1 | 1, s, mid + 1, r, R, check);
            if(pos != -1) return pos;
            return dfs(dfs, u << 1, s, l, mid, R, check);
        };
        S temp_s = e();
        return dfs(dfs, 1, temp_s, 1, _n, R, check);
    }
};

signed main() {

    int n, m;

    cin >> n >> m;
    cin >> mod;

    vector<S> a(n + 1);
    for(int i = 1; i <= n; i++) cin >> a[i].sum, a[i].sz = 1;
    
    LazySegmentTree<S, op, e, F, mapping, comp, id> T(a);

    while(m--) {
        int op; cin >> op;
        if(op == 1) {
            int l, r; cin >> l >> r;
            int mul; cin >> mul;
            T.apply(l, r, {0, mul});
        } else if(op == 2) {
            int l, r, add; cin >> l >> r >> add;
            T.apply(l, r, {add, 1});
        } else if(op == 3) {
            int l, r; cin >> l >> r;
            S x = T.prod(l, r);
            cout << x.sum << "\n";
        } else if(op == 4) {
            int l, val; cin >> l >> val;
            int pos = T.min_left(l, [&] (S x) {
                return x.sum <= val;
            });
            cout << pos << "\n";
        }
    }

    return 0;
}