#include<bits/stdc++.h>
using namespace std;
int mod;
struct info {
    int sum;
    int sz;
};
struct F {
    int add, mul;
};
inline info merge(info a, info b) {
    return {(a.sum + b.sum) % mod, (a.sz + b.sz) % mod};
}
inline info mapping(F f, info a) {
    return {(1ll * a.sum * f.mul % mod + 1ll * f.add * a.sz % mod) % mod, a.sz};
}
// f作用于g
inline F comp(F f, F g) {
    return {(1ll * g.add * f.mul % mod + f.add) % mod, 1ll * f.mul * g.mul % mod};
}
const int N = 1e5 + 1;
struct segmenttree {

    info tr[N << 2];
    F tag[N << 2];
    int treel[N << 2], treer[N << 2];
    void build(int u, int l, int r, vector<int> &a) {
        treel[u] = l;
        treer[u] = r;
        tag[u] = {0, 1};
        if(l == r) {
            tr[u].sum = a[l];
            tr[u].sz = 1;
            return;
        }
        int mid = l + r >> 1;
        build(u << 1, l, mid, a);
        build(u << 1 | 1, mid + 1, r, a);
        tr[u] = merge(tr[u << 1], tr[u << 1 | 1]);
    }

    void pushdown(int u) {
        tr[u << 1] = mapping(tag[u], tr[u << 1]);
        tr[u << 1 | 1] = mapping(tag[u], tr[u << 1 | 1]);
        tag[u << 1] = comp(tag[u], tag[u << 1]);
        tag[u << 1 | 1] = comp(tag[u], tag[u << 1 | 1]);
        tag[u] = {0, 1};
    }

    info query(int u, int L, int R) {
        if(treel[u] > R || treer[u] < L) return {0, 0};
        if(treel[u] >= L && treer[u] <= R) {
            return tr[u];
        }
        pushdown(u);
        int mid = treel[u] + treer[u] >> 1;
        return merge(
            query(u << 1, L, R),
            query(u << 1 | 1, L, R)
        );
    }

    void modify(int u, int L, int R, F f) {
        if(treel[u] > R || treer[u] < L) return;
        if(treel[u] >= L && treer[u] <= R) {
            tr[u] = mapping(f, tr[u]);
            tag[u] = comp(f, tag[u]);
            return;
        }
        pushdown(u);
        modify(u << 1, L, R, f);
        modify(u << 1 | 1, L, R, f);
        tr[u] = merge(tr[u << 1], tr[u << 1 | 1]);
    }

    template<typename Check>
    int max_right(int u, info &s, int l, int r, int L, Check check) {
        if(r < L) return -1;
        if(l >= L) {
            info ss = merge(s, tr[u]);
            if(check(ss)) {
                s = ss;
                return -1;
            }
            if(l == r) return l;
        }
        pushdown(u);
        int mid = l + r >> 1;
        int pos = max_right(u << 1, s, l, mid, L, check);
        if(pos != -1) return pos;
        return max_right(u << 1 | 1, s, mid + 1, r, L, check);
    }


};

segmenttree segt;
signed main() {
    int n, m;

    cin >> n >> m;
    cin >> mod;

    vector<int> a(n + 1);
    for(int i = 1; i <= n; i++) cin >> a[i];
    
    segt.build(1, 1, n, a);

    while(m--) {
        int op; cin >> op;
        if(op == 1) {
            int l, r; cin >> l >> r;
            int mul; cin >> mul;
            segt.modify(1, l, r, {0, mul});
        } else if(op == 2) {
            int l, r, add; cin >> l >> r >> add;
            segt.modify(1, l, r, {add, 1});
        } else if(op == 3) {
            int l, r; cin >> l >> r;
            info x = segt.query(1, l, r);
            cout << x.sum << "\n";
        } else if(op == 4) {
            int l, val; cin >> l >> val;
            info u = {0, 0};
            int pos = segt.max_right(1, u, 1, n, l, [&] (info x) {
                return x.sum <= val;
            });
            cout << pos << "\n";
        }
    }

    return 0;
}