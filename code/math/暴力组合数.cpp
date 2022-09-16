int binom(int n, int k) {
    if(n < 0 || k < 0 || k > n) { return 0; }
    k = min(n - k, k);
    int u = 1, v = 1;
    for(int i = 0; i < k; i++) {
    	v = v * (i + 1) % mod;
    	u = u * (n - i) % mod;
	}
	return u * ksm(v, mod - 2) % mod;
}