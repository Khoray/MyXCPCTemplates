#include<bits/stdc++.h>
using namespace std;
// 运算符
template<typename T>
/// @brief 二元运算符
/// @tparam a, b, sg
/// @return a (sg) b
function<T(T, T, char)> ca = [&] (T a, T b, char sg) {
    // 操作取模在这里改
    if(sg == '+') return a + b;
    if(sg == '-') return a - b;
    if(sg == '*') return a * b;
    if(sg == '/') return a / b;
    return (T)0;
};
// 读入取模在这里改
/// @param 要读入的字符串s 和当前指针pt，都要是引用类型
// 最后pt停在数字的最后一位
function<int(string&, int&)> readint = [&] (string &s, int &pt) {
    int ret = 0, fl = 1;
    if(s[pt] == '-') fl = -1, pt++;
    while(pt < s.length() && isdigit(s[pt])) {
        ret = ret * 10 + s[pt++] - '0';
    }
    pt--;
    return ret * fl;
};
function<double(string&, int&)> readdouble = [&] (string &s, int &pt) {
    int fl = 1;
    string ret;
    if(s[pt] == '-') fl = -1, pt++;
    while(pt < s.length() && (isdigit(s[pt]) || s[pt] == '.')) {
        ret += s[pt++];
    }
    pt--;
    return stod(ret) * fl;
};

// priority['('] = 0;
template<typename T>
T calc(string s, map<char, int> priority, function<T(T, T, char)> calculate, function<T(string&, int&)> readnum) {
    int pt = 0;
    string tmp;
    for(int i = 0; i < s.length(); i++) if(s[i] != ' ') tmp += s[i];
    s = tmp;
    vector<T> stk;
    vector<char> sgn;
    auto calc_single = [&] () {
        char sg = sgn.back();
        sgn.pop_back();
        T v = stk.back();
        stk.pop_back();
        T u = stk.back();
        stk.pop_back();
        stk.push_back(calculate(u, v, sg));
    };
    int mode = 0;
    // 0 is number, 1 is sign
    while(pt < s.length()) {
        if(mode == 0) {
            if(isdigit(s[pt]) || s[pt] == '-') {
                stk.push_back(readnum(s, pt));// 读数字
                mode = 1;
            } else if(s[pt] == '(') {
                sgn.push_back(s[pt]);
            }
        } else {
            if(s[pt] != ')') {
                while(sgn.size() && priority[sgn.back()] >= priority[s[pt]]) {
                    calc_single();
                }
                sgn.push_back(s[pt]);
                mode = 0;
            } else if(s[pt] == ')') {
                while(sgn.back() != '(') {
                    calc_single();
                }
                sgn.pop_back();
            }
        }
        pt++;
    }
    while(sgn.size()) {
        calc_single();
    }
    assert(stk.size() == 1);
    return stk[0];
};

int main() {
    // 设定优先级
    map<char, int> prio;
    prio['+'] = prio['-'] = 1;
    prio['*'] = 2;
    prio['/'] = 2;
    string s; cin >> s;
    cout << calc<double>(s, prio, ca<double>, readdouble);
    return 0;
}