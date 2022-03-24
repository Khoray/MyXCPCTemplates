// ax + by = gcd(a, b)
// return gcd(a, b)
// all sol: x = x_0 + k[a, b] / a, y = y_0 - k[a, b] / b;
int exgcd(int &x, int &y, int a, int b) {
    if(!b) {
        x = 1;
        y = 0;
        return a;
    }
    int ret = exgcd(x, y, b, a % b);
    int t = x;
    x = y;
    y = t - a / b * y;
    return ret;
}
