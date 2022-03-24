// all in which f[pq] = f[p] * f[q] (gcd(p, q) = 1)
const int N = 1e7 + 5;
int pri[N / 5], notpri[N], prinum, minpri_cnt[N], f[N];

int calc_f(int val, int power) {
	
}

void init_pri() {
    for(int i = 2; i < N; i++) {
        if(!notpri[i]) pri[++prinum] = i, minpri_cnt[i] = 1, f[i] = calc_f(i, 1);
        for(int j = 1; j <= prinum && pri[j] * i < N; j++) {
            notpri[pri[j] * i] = pri[j];
            if(i % pri[j] == 0) {
            	minpri_cnt[pri[j] * i] = minpri_cnt[i] + 1;
            	f[pri[j] * i] = f[i] / calc_f(pri[j], minpri_cnt[i]) * calc_f(pri[j], minpri_cnt[i] + 1);
                break;
            }
            minpri_cnt[pri[j] * i] = 1;
            f[pri[j] * i] = f[i] * calc_f(pri[j], 1);
        }
    }
}