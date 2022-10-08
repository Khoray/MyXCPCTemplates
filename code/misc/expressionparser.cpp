// 运算符
function<int(int, int, char)> ca = [&] (int a, int b, char sg) {
    if(sg == '+') return (a + b) % 10000;
    if(sg == '-') return a - b;
    if(sg == '*') return (a * b) % 10000;
    return 0;
};
// 设定优先级
map<char, int> prio;
prio['+'] = prio['-'] = 1;
prio['*'] = 2;
// calculate 传入前两个是操作数，第三个是操作符   priority['('] = 0;
int parse_expression(string s, map<char, int> priority, function<int(int, int, char)> calculate) {
    int pt = 0;
    string tmp;
    for(int i = 0; i < s.length(); i++) if(s[i] != ' ') tmp += s[i];
    s = tmp;
    vector<int> stk;
    vector<char> sgn;
    // 读入取模在这里改
    auto readint = [&]() {
        int ret = 0, fl = 1;
        if(s[pt] == '-') fl = -1, pt++;
        while(pt < s.length() && isdigit(s[pt])) {
            ret = ret * 10 + s[pt] - '0';
            pt++;
        }
        return ret * fl;
    };
    auto calc_single = [&] () {
        char sg = sgn.back();
        sgn.pop_back();
        int v = stk.back();
        stk.pop_back();
        int u = stk.back();
        stk.pop_back();
        stk.push_back(calculate(u, v, sg));
    };
    int mode = 0;
    // 0 is number, 1 is sign
    while(pt < s.length()) {
        if(mode == 0) {
            if(isdigit(s[pt]) || s[pt] == '-') {
                int num = readint();
                stk.push_back(num);
                mode = 1;
            } else if(s[pt] == '(') {
                sgn.push_back(s[pt]);
                pt++;
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
            pt++;
        }
    }
    while(sgn.size()) {
        calc_single();
    }
    assert(stk.size() == 1);
    return stk[0];
};