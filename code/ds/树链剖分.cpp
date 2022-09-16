vector<int> dep(n + 1), parent(n + 1), siz(n + 1), top(n + 1), son(n + 1);
function<void(int, int)> dfs1 = [&] (int u, int pa) -> void {
    siz[u] = 1;
    parent[u] = pa;
    dep[u] = dep[pa] + 1;
    for(int &v : adj[u]) {
        if(v == parent[u]) continue;
        dfs1(v);
        siz[u] += siz[v];
        if(siz[v] > siz[son[u]]) son[u] = v;
    }
};
function<void(int, int)> dfs2 = [&] (int u, int tp) -> void {
    if(son[u]) dfs2(u, tp);
    for(int &v : adj[u]) {
        if(v == parent[u] || v == son[u]) continue;
        dfs2(v, v);
    }
};
dfs(1, 0);
dfs2(1, 1);
function<int(int, int)> lca = [&] (int u, int v) {
    while(top[u] != top[v]) {
        if(dep[top[u]] > dep[top[v]]) {
            u = parent[top[u]];
        } else {
            v = parent[top[v]];
        }
    }
    if(dep[u] < dep[v]) {
        return u;
    } else {
        return v;
    }
};
