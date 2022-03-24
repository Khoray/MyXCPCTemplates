int inv(int x, int p) {
    int y, k;
    int gcd = exgcd(y, k, x, p);
    int moder = p / gcd;
    return (y % moder + moder) % moder;
}