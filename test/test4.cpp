#include <bits/stdc++.h>
using namespace std;

template<typename T>
inline void Read(T &n){
	char ch; bool flag=false;
	while(!isdigit(ch=getchar()))if(ch=='-')flag=true;
	for(n=ch^48;isdigit(ch=getchar());n=(n<<1)+(n<<3)+(ch^48));
	if(flag)n=-n;
}

const int MAXN = 64005;
const int MOD = 998244353;
const int G = 3;

inline int inc(int a, int b){
	a += b;
	if(a>=MOD) a -= MOD;
	return a;
}

inline void iinc(int &a, int b){a = inc(a,b);}

inline int dec(int a, int b){
	a -= b;
	if(a<0) a += MOD;
	return a;
}

inline void ddec(int &a, int b){a = dec(a,b);}

inline int ksm(int base, int k=MOD-2){
	int res = 1;
	while(k){
		if(k&1)
			res = 1ll*res*base%MOD;
		base = 1ll*base*base%MOD;
		k >>= 1;
	}
	return res;
}

typedef vector<int> poly;

int tr[MAXN<<2], wn[MAXN<<2];
inline int prework(int n){
	int len = 1; while(len<=n) len<<=1;
	for(register int i=0; i<=len; i++) tr[i] = (tr[i>>1]>>1)|((i&1)?len>>1:0);
	wn[0] = 1; wn[1] = ksm(G,(MOD-1)/len); for(register int i=2; i<=len; i++) wn[i] = 1ll*wn[i-1]*wn[1]%MOD;
	return len;
}

poly f, g;
inline void ntt(poly &f, int n, int flag){
	for(register int i=0; i<n; i++) if(i<tr[i]) swap(f[i],f[tr[i]]);
	for(register int len=2; len<=n; len<<=1){
		int base = flag*n/len;
		for(register int l=0; l<n; l+=len){
			int now = flag==1?0:n;
			for(register int i=l; i<l+len/2; i++){
				int tmp = 1ll*f[i+len/2]*wn[now]%MOD;
				f[i+len/2] = dec(f[i],tmp);
				f[i] = inc(f[i],tmp);
				now += base;
			}
		}
	}
	if(flag==-1){
		int invn = ksm(n);
		for(register int i=0; i<n; i++) f[i] = 1ll*f[i]*invn%MOD;
	}
}

void poly_inv(poly f, poly &res, int n){
	static poly tmp;
	if(n==1) return res.resize(1), res[0] = ksm(f[0]), void();
	poly_inv(f,res,n+1>>1);
	int len = prework(n<<1);
	tmp.resize(len); res.resize(len);
	for(register int i=0; i<n; i++) tmp[i] = f[i];
	ntt(tmp,len,1); ntt(res,len,1);
	for(register int i=0; i<len; i++) res[i] = dec(inc(res[i],res[i]),1ll*tmp[i]*res[i]%MOD*res[i]%MOD);
	ntt(res,len,-1);
	tmp.resize(0); res.resize(n);
}

int a[MAXN];

int n, m;

poly A[MAXN<<2], B[MAXN<<2], invA, F;
inline int lc(int x){return x<<1;}
inline int rc(int x){return x<<1|1;}

void solve1(int x, int l, int r){
	if(l==r) return A[x].resize(2), A[x][0] = 1, A[x][1] = dec(0,a[l]), void();
	int mid = l+r >> 1;
	solve1(lc(x),l,mid); solve1(rc(x),mid+1,r);
	int num = A[lc(x)].size()+A[rc(x)].size()-2;
	int len = prework(r-l+1);
	f = A[lc(x)]; g = A[rc(x)];
	f.resize(len); g.resize(len); A[x].resize(len);
	ntt(f,len,1); ntt(g,len,1);
	for(register int i=0; i<len; i++) A[x][i] = 1ll*f[i]*g[i]%MOD;
	ntt(A[x],len,-1);
	f.resize(0); g.resize(0); A[x].resize(r-l+2);
}

int ans[MAXN];
void solve2(int x, int l, int r){
    cout << l << ' ' << r << '\n';
    for(auto tp:B[x]) cout << tp << ' ';
    cout << '\n';
	if(l==r) return ans[l] = B[x][0], void();
	int mid = l+r >> 1;
	int len = prework(r-l+1<<1);
    for(auto tp:A[rc(x)]) cout << tp << ' ';
    cout << '\n';
	reverse(A[lc(x)].begin(),A[lc(x)].end());
	reverse(A[rc(x)].begin(),A[rc(x)].end());
	A[lc(x)].resize(len); A[rc(x)].resize(len); B[x].resize(len); B[lc(x)].resize(len); B[rc(x)].resize(len);
	ntt(A[lc(x)],len,1); ntt(A[rc(x)],len,1), ntt(B[x],len,1);
	for(register int i=0; i<len; i++) B[lc(x)][i] = 1ll*A[rc(x)][i]*B[x][i]%MOD, B[rc(x)][i] = 1ll*A[lc(x)][i]*B[x][i]%MOD;
	ntt(B[lc(x)],len,-1); ntt(B[rc(x)],len,-1);
	for(register int i=0; i<=mid-l+1; i++) B[lc(x)][i] = B[lc(x)][i+r-mid]; B[lc(x)].resize(mid-l+2);
	for(register int i=0; i<=r-mid; i++) B[rc(x)][i] = B[rc(x)][i+mid-l+1]; B[rc(x)].resize(r-mid+1);
	solve2(lc(x),l,mid); solve2(rc(x),mid+1,r);
}

int main(){
	Read(n); Read(m);
	F.resize(n+1);
	for(register int i=0; i<=n; i++) Read(F[i]);
	for(register int i=0; i<m; i++) Read(a[i]);
	n = max(n,m-1);
	solve1(1,0,n); poly_inv(A[1],invA,n+1);
	reverse(invA.begin(),invA.end());
    // cout << "invA:"; for(auto tp:invA) cout << tp << ' ';
    // cout << n << '\n';
	int len = prework(n<<1);
	F.resize(len); invA.resize(len); B[1].resize(len);
	ntt(F,len,1); ntt(invA,len,1);
	for(register int i=0; i<len; i++) B[1][i] = 1ll*F[i]*invA[i]%MOD;
	ntt(B[1],len,-1); for(register int i=0; i<=n; i++) B[1][i] = B[1][i+n]; B[1].resize(n+1);
	solve2(1,0,n);
	for(register int i=0; i<m; i++) printf("%d\n",ans[i]);
	return 0;
}