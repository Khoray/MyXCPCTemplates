const int N = 1005;
template<int N>
struct Matrix {
    bitset<N> mat[N];
    Matrix() { }
    bitset<N> &operator [] (int idx) { return mat[idx]; }
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
					swap(a[row], a[i]);
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
				a[i] ^= a[row];
				b[i] ^= b[row];
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
			ans[col] = 1;
			all++;
		} else {
			row--;
			ans[col] = b[row];
			for(int i = col + 1; i <= n; i++) {
				if(a[row][i]) ans[col] ^= ans[i]; 
			}
		}
	}
	return all;
}