#include<bits/stdc++.h>
using namespace std;

using ll = long long;

template<class T>
struct MediumStructure {
    T sumup, sumdown;
    multiset<T> up, down;

    MediumStructure(): sumup(0), sumdown(0) {}

    void balance() {
        while(down.size() > up.size() + 1) {
            T x = *down.begin();
            sumdown -= x;
            sumup += x;
            down.erase(down.find(x));
            up.insert(x);
        }
        while(up.size() > down.size()) {
            T x = *up.rbegin();
            sumup -= x;
            sumdown += x;
            up.erase(prev(up.end()));
            down.insert(x);
        }
    }

    void add(T x) {
        if(down.empty() || x >= *down.begin()) {
            down.insert(x);
            sumdown += x;
        } else {
            up.insert(x);
            sumup += x;
        }
        balance();
    }

    void remove(T x) {
        assert(!down.empty());
        if(x >= *down.begin()) {
            auto it = down.find(x);
            assert(it != down.end());
            down.erase(it);
            sumdown -= x;
        } else {
            auto it = up.find(x);
            assert(it != up.end());
            up.erase(up.find(x));
            sumup -= x;
        }
        balance();
    }
};

template<class T>
struct StaticKthMaxStructure {
    int k;
    T sumup, sumdown;
    multiset<T> up, down;

    StaticKthMaxStructure(int k): k(k), sumup(0), sumdown(0) {}

    void balance() {
        while(down.size() > k) {
            T x = *down.begin();
            sumdown -= x;
            sumup += x;
            down.erase(down.find(x));
            up.insert(x);
        }
        while(up.size() && down.size() < k) {
            T x = *up.rbegin();
            sumup -= x;
            sumdown += x;
            up.erase(prev(up.end()));
            down.insert(x);
        }
    }

    void add(T x) {
        if(down.empty() || x >= *down.begin()) {
            down.insert(x);
            sumdown += x;
        } else {
            up.insert(x);
            sumup += x;
        }
        balance();
    }

    void remove(T x) {
        assert(!down.empty());
        if(x >= *down.begin()) {
            auto it = down.find(x);
            assert(it != down.end());
            down.erase(it);
            sumdown -= x;
        } else {
            auto it = up.find(x);
            assert(it != up.end());
            up.erase(up.find(x));
            sumup -= x;
        }
        balance();
    }
};

void solve() {
    StaticKthMaxStructure<int> ss(3);
    ss.add(1);
    ss.add(2);
    ss.add(3);
    ss.add(4);
    ss.add(5);
    ss.add(6);
    ss.add(7);
    cout << *ss.down.begin() << " " << ss.sumup << '\n';
    ss.remove(7);
    cout << *ss.down.begin() << '\n';
    


}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
    int t = 1; // cin >> t;
    while(t--) solve();
    return 0;
}