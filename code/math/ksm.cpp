int ksm(int a, int b = mod - 2, int MOD_KSM = mod) {
    int ret = 1;
    while(b) {
        if(b & 1) {
            ret = (ll) ret * a % MOD_KSM;
        }
        a = (ll) a * a % MOD_KSM;
        b >>= 1;
    }
    return ret;
}