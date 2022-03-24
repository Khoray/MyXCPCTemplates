const int N = 31;
int C[N][N];

void init_C() {
	C[0][0] = C[1][0] = C[1][1] = 1;
	for(int i = 2; i < N; i++) {
		C[i][0] = 1;
		for(int j = 1; j <= i; j++)
			C[i][j] = (C[i - 1][j] + C[i - 1][j - 1]) % mod;
	}
}

int binom(int n, int k) {
	if(n < 0 || k < 0 || n < k) return 0;
	return C[n][k];
}
