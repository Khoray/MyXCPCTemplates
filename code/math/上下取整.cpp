// b must be positive integer

int updiv(int a, int b) {
    return a > 0 ? (a + b - 1) / b : a / b;
}

int downdiv(int a, int b) {
    return a > 0 ? a / b : (a - b + 1) / b;
}
