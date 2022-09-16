#define mul(a, b) (ll(a) * (b) % mod)
#define add(a, b) (((a) += (b)) >= mod ? (a) -= mod : 0) // (a += b) %= P
#define dec(a, b) (((a) -= (b)) < 0 ? (a) += mod: 0)  // ((a -= b) += P) %= P

const int N = 1005;
const int mod = 1e9 + 7;
template<int N>
struct Matrix {
    int mat[N][N];
    Matrix() { }
    int* operator [] (int idx) { return mat[idx]; }
};
template<int N>
int guass(int n, int equ, Matrix<N> a, vector<int> b, vector<int> &ans) {
    fill(ans.begin(), ans.end(), 0);
    vector<int> fre(n + 1);
	int row, col;
	for(row = 1, col = 1; col <= n; col++) {
		if(!a[row][col]) {
			int sw = 0;
			for(int i = row + 1; i <= equ; i++) {
				if(a[i][col]) {
					for(int j = col; j <= n; j++) {
						swap(a[row][j], a[i][j]);
					}
					swap(b[row], b[i]);
					sw = 1;
					break;
				}
			}
			if(!sw) {
				fre[col] = 1;
				continue;
			}
		}
		for(int i = row + 1; i <= equ; i++) {
			if(a[i][col]) {
				int k = a[i][col] * ksm(a[row][col]) % mod;
				for(int j = col; j <= n; j++) {
					dec(a[i][j], a[row][j] * k % mod);
				}
				dec(b[i], b[row] * k % mod);
			}
		}
		row++;
	}
	if(row <= equ) {
		for(int i = row; i <= equ; i++) {
			if(b[i]) return -1;
		}
	}
	int all = 0;
	for(col = n; col >= 1; col--) {
		if(fre[col]) {
			ans[col] = 0;
			all++;
		} else {
			row--;
			ans[col] = b[row];
			for(int i = col + 1; i <= n; i++) {
				dec(ans[col], ans[i] * a[row][i] % mod); 
			}
			mul(ans[col], ksm(a[row][col]));
		}
	}
	return all;
}