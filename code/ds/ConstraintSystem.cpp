struct ConstraintSystem {
	using ll = long long;
	vector<vector<pair<int, int>>> adj;
	int _n;

	ConstraintSystem(int n): _n(n), adj(n + 1) {
		for(int i = 1; i <= n; i++) {
			adj[0].emplace_back(i, 0);
		}
	}

	void leq_than(int u, int v, int w) { // xu - xv <= w
		adj[v].emplace_back(u, w);
	}

	void geq_than(int u, int v, int w) { // xu - xv >= w
		adj[u].emplace_back(v, -w);
	}

	void less_than(int u, int v, int w) { // xu - xv < w
		leq_than(u, v, w - 1);
	}

	void greater_than(int u, int v, int w) { // xu - xv > w
		geq_than(u, v, w + 1);
	}

	void eq(int u, int v, int w) {
		leq_than(u, v, w);
		geq_than(u, v, w);
	}

	// get all result: solve(0, dis);
	// get max{x_u - x_v}: solve(v, dis); ans = dis[u];
	int solve(int start, vector<ll> &dis) {
		queue<int> q;
		dis.resize(_n + 1);
		fill(dis.begin(), dis.end(), (ll) 1e18);
		dis[start] = 0;
		q.push(start);
		vector<int> inq(_n + 1), tot(_n + 1);
		inq[start] = 1;
		while(q.size()) {
			int u = q.front();
			q.pop();
			inq[u] = 0;
			tot[u]++;
			if(tot[u] > _n) {
				return 0;
			}
			for(auto [v, w] : adj[u]) {
				if(dis[v] > dis[u] + w) {
					dis[v] = dis[u] + w;
					if(!inq[v]) {
						inq[v] = 1;
						q.push(v);
					}
				}
			}
		}
		return 1;
	}
};