void solve() {
    // d 是间隔输出（对题目无影响
	int n, d; cin >> n >> d;
	vector<int> pf; // 质因数分解
	int pn = phi[n];
	while(notpri[pn]) {
		int now = notpri[pn];
		pf.push_back(now);
		while(pn % now == 0) pn /= now;
	}
	if(pn != 1) {
		pf.push_back(pn);
	}
	int cnt = 0;
	int ming = -1;
	vector<int> ans, vis(n); // 记录答案
    // 找到最小的原根 min_g
	for(int i = 1; i < n; i++) {
		if(__gcd(i, n) != 1) continue;
		int judge = 1;
		for(auto &p : pf) {
			if(ksm(i, phi[n] / p, n) == 1) {
				judge = 0;
				break;
			}
		}
		if(judge) {
			ming = i;
			break;
		}
	}
    // 还原出所有原根 g
	if(ming > 0) {
		for(int i = 1; i < n; i++) {
			if(__gcd(i, phi[n]) == 1) {
				int cur = ksm(ming, i, n);
				if(!vis[cur]) {
					vis[cur] = 1;
					ans.push_back(cur);
				}
			}
		}
	}
    // 排序输出所有原根
	cout << ans.size() << '\n';
	sort(ans.begin(), ans.end());
	for(auto as : ans) {
		cnt++;
		if(cnt % d == 0) {
			cout << as << ' ';
		}
	}
	cout << '\n';
}