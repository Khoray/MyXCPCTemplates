int LagrangeInterpolation(vector<int> &y, int l, int r, int n) {
    if(n <= r && n >= l) return y[n];
	vector<int> lg(r - l + 3), rg(r - l + 3);
    int ret = 0;
	lg[0] = 1;
	rg[r - l + 2] = 1;
    for(int i = l; i <= r; i++) {
        lg[i - l + 1] = (ll) lg[i - l] * (n - i) % mod;
    }
	for(int i = r; i >= l; i--) {
        rg[i - l + 1] = (ll) rg[i - l + 2] * (n - i) % mod;
    }
    for(int i = l; i <= r; i++) {
        if(r - i & 1) {
            ret = (ret - (ll) y[i] * lg[i - l] % mod * rg[i - l + 2] % mod * facinv[i - l] % mod * facinv[r - i] % mod + mod) % mod;
        } else {
            ret = (ret + (ll) y[i] * lg[i - l] % mod * rg[i - l + 2] % mod * facinv[i - l] % mod * facinv[r - i] % mod) % mod;
        }
    }
    return ret;
}
