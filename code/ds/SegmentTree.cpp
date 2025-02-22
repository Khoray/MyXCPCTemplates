/*
https://ac.nowcoder.com/acm/contest/88888/F
单点修改、区间查询
*/

#include<bits/stdc++.h>
using namespace std;

int op(int a, int b) {
    return min(a, b);
}

int e() { return INT_MAX; }

template<class S,
         auto op,
         auto e>
struct SegmentTree {
    static_assert(std::is_convertible_v<decltype(op), std::function<S(S, S)>>,
                  "op must work as S(S, S)");
    static_assert(std::is_convertible_v<decltype(e), std::function<S()>>,
                  "e must work as S()");

    int _n;
    vector<S> a;

    // 下标从 1 开始，v[1..n]
    SegmentTree(int n): a(n << 2, e()), _n(n) {}

    // 传入的 v 是 [0..n] 的, build 时使用 v[1.._n], _n=n=v.size()-1
    SegmentTree(const vector<S> &v): _n((int) v.size() - 1), a(_n << 2, e()) {
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

    // return op(a[l], .., a[r])
    S prod(int L, int R) {
        auto dfs = [&] (auto&& dfs, int u, int l, int r, int L, int R) -> S {
            if(l > R || L > r) return e();
            if(l >= L && r <= R) return a[u];
            int mid = (l + r) / 2;
            return op(
                dfs(dfs, u << 1, l, mid, L, R),
                dfs(dfs, u << 1 | 1, mid + 1, r, L, R)
            );
        };
        return dfs(dfs, 1, 1, _n, L, R);
    }

    // set v[p] := x
    void set(int p, S x) {
        auto dfs = [&] (auto&& dfs, int u, int l, int r, int p, S x) -> void {
            if(l == r) {
                a[u] = x;
                return;
            }
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
            int mid = (l + r) / 2;
            int pos = dfs(dfs, u << 1 | 1, s, mid + 1, r, R, check);
            if(pos != -1) return pos;
            return dfs(dfs, u << 1, s, l, mid, R, check);
        };
        S temp_s = e();
        return dfs(dfs, 1, temp_s, 1, _n, R, check);
    }
};

void solve() {
    int n; cin >> n;
    vector<vector<int>> arr(n + 1);
    vector<int> sz(n + 1), sum_sz(n + 1);
    for(int i = 1; i <= n; i++) {
        cin >> sz[i];
        arr[i].resize(sz[i] + 1);
        for(int j = 1; j <= sz[i]; j++) {
            cin >> arr[i][j];
        }
    }
    for(int i = 1; i <= n; i++) {
        sum_sz[i] = sum_sz[i - 1] + sz[i];
    }
    vector<int> a(sum_sz[n] + 1);
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= sz[i]; j++) {
            a[sum_sz[i - 1] + j] = arr[i][j];
        }
    }
    SegmentTree<int, op, e> T(a);
    int q; cin >> q;

    while(q--) {
        int op; cin >> op;
        if(op == 1) {
            int u, v, x; cin >> u >> v >> x;
            T.set(sum_sz[u - 1] + v, x);
        } else {
            int u; cin >> u;
            cout << T.prod(1, sum_sz[u]) << '\n';
        }
    }

}

signed main() {
    #ifndef _DBG
        ios::sync_with_stdio(false);
        cin.tie(0);
    #endif
    int t = 1; // cin >> t;
    while(t--) {
        solve();
    }
}