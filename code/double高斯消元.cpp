const int N = 1005;
template<int N>
struct Matrix {
    bitset<N> mat[N];
    Matrix() { }
    bitset<N> &operator [] (int idx) { return mat[idx]; }
};
template<int N>
int guass(int n, int equ, Matrix<N> a, vector<double> b, vector<double> &ans) {
    fill(ans.begin(), ans.end(), 0);
    vector<int> fre(n + 1);
	int row, col;
	for(row = 1, col = 1; col <= n; col++) {
		double mx = fabs(a[row][col]);
		int mxp = row;
		for(int i = row + 1; i <= equ; i++) {
			if(fabs(a[row][col]) > mx) {
				mx = fabs(a[row][col]);
				mxp = i;
			}
		}
		if(mxp != row) {
			for(int i = col; i <= n; i++) {
				swap(a[row][i], a[mxp][i]);
			}
			swap(b[row], b[mxp]);
		}
		if(fabs(a[row][col]) < eps) {
			fre[col] = 1;
		}
		for(int i = row + 1; i <= equ; i++) {
			if(fabs(a[i][col]) > eps) {
				double k = a[i][col] / a[row][col];
				for(int j = col; j <= n; j++) {
					a[i][j] -= a[row][j] * k;
				}
				b[i] -= b[row] * k;
			}
		}
		row++;
	}
	// 判断解是否存在
	if(row <= equ) {
		for(int i = row; i <= equ; i++) {
			if(fabs(b[i]) > eps) return -1;
		}
	}

	// 回代求解
	int all = 0;
	for(col = n; col >= 1; col--) {
		if(fre[col]) {
			ans[col] = 0;
			all++;
		} else {
			row--;
			ans[col] = b[row];
			for(int i = col + 1; i <= n; i++) {
				ans[col] -= ans[i] * a[row][i]; 
			}
			ans[col] /= a[row][col];
		}
	}
	return all;
}