/*
https://www.luogu.com.cn/problem/P3373
区间修改、区间查询、动态开点
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
         auto id,
         typename P>
struct DynamicLazySegmentTree {
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

    int _n, node_cnt;
    vector<S> a;
    vector<F> tag;
    vector<int> lc, rc;

    // default_prod 是没有经过任何操作的原始数组建立的线段树 (l,r) 这个结点的值
    P default_prod;

    // 下标从 1 开始，v[1..n], node_id=1 是 root
    DynamicLazySegmentTree(int n, P default_prod): _n(n), node_cnt(1), a(2), tag(2), lc(2), rc(2), default_prod(default_prod) {
        static_assert(std::is_convertible_v<decltype(default_prod), std::function<S(int, int)>>,
                  "default_prod must work as S(int, int)");
        a[1] = default_prod(1, n);
        tag[1] = id();
        lc[1] = rc[1] = 0;
    }

    int new_node(int l, int r) {
        ++node_cnt;
        // if(node_cnt >= (int) a.size()) {
        //     a.resize(node_cnt + 1000, e());
        //     tag.resize(node_cnt + 1000, id());
        //     lc.resize(node_cnt + 1000);
        //     rc.resize(node_cnt + 1000);
        // }
        a.push_back(e());
        tag.push_back(id());
        lc.push_back(0);
        rc.push_back(0);
        a[node_cnt] = default_prod(l, r);
        tag[node_cnt] = id();
        lc[node_cnt] = rc[node_cnt] = 0;
        return node_cnt;
    }

    void pushdown(int u) {
        a[lc[u]] = mapping(tag[u], a[lc[u]]);
        a[rc[u]] = mapping(tag[u], a[rc[u]]);
        tag[lc[u]] = composition(tag[u], tag[lc[u]]);
        tag[rc[u]] = composition(tag[u], tag[rc[u]]);
        tag[u] = id();
    }

    // return op(a[l], .., a[r])
    S prod(int L, int R) {
        auto dfs = [&] (auto&& dfs, int u, int l, int r, int L, int R) -> S {
            if(l > R || L > r) return e();
            if(l >= L && r <= R) return a[u];
            int mid = (l + r) / 2;
            if(!lc[u]) lc[u] = new_node(l, mid);
            if(!rc[u]) rc[u] = new_node(mid + 1, r);
            pushdown(u);
            return op(
                dfs(dfs, lc[u], l, mid, L, R),
                dfs(dfs, rc[u], mid + 1, r, L, R)
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
            int mid = (l + r) / 2;
            if(!lc[u]) lc[u] = new_node(l, mid);
            if(!rc[u]) rc[u] = new_node(mid + 1, r);
            pushdown(u);
            dfs(dfs, lc[u], l, mid, L, R, f);
            dfs(dfs, rc[u], mid + 1, r, L, R, f);
            a[u] = op(a[lc[u]], a[rc[u]]);
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
            int mid = (l + r) / 2;
            if(!lc[u]) lc[u] = new_node(l, mid);
            if(!rc[u]) rc[u] = new_node(mid + 1, r);
            pushdown(u);
            if(p <= mid) dfs(dfs, lc[u], l, mid, p, x);
            else dfs(dfs, rc[u], mid + 1, r, p, x);
            a[u] = op(a[lc[u]], a[rc[u]]);
        };
        return dfs(dfs, 1, 1, _n, p, x);
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
            int mid = (l + r) / 2;
            if(!lc[u]) lc[u] = new_node(l, mid);
            if(!rc[u]) rc[u] = new_node(mid + 1, r);
            pushdown(u);
            int pos = dfs(dfs, lc[u], s, l, mid, L, check);
            if(pos != -1) return pos;
            return dfs(dfs, rc[u], s, mid + 1, r, L, check);
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
            int mid = (l + r) / 2;
            if(!lc[u]) lc[u] = new_node(l, mid);
            if(!rc[u]) rc[u] = new_node(mid + 1, r);
            pushdown(u);
            int pos = dfs(dfs, rc[u], s, mid + 1, r, R, check);
            if(pos != -1) return pos;
            return dfs(dfs, lc[u], s, l, mid, R, check);
        };
        S temp_s = e();
        return dfs(dfs, 1, temp_s, 1, _n, R, check);
    }
};

class P {
public:
    S operator () (int l, int r) {
        return e();
    }
};

signed main() {
    int n, m;

    cin >> n >> m;
    cin >> mod;

    vector<S> a(n + 1);
    for(int i = 1; i <= n; i++) cin >> a[i].sum, a[i].sz = 1;
    DynamicLazySegmentTree<S, op, e, F, mapping, comp, id, P> T(n, P());
    for(int i = 1; i <= n; i++) {
        T.set(i, a[i]);
    }
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