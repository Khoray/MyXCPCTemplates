#include<bits/stdc++.h>
using namespace std;

#define out(args...) { cout << "Line " << __LINE__ << ": [" << #args << "] = ["; debug(args); cout << "]\n"; }

template<typename T> void debug(T a) { cout << a; }

template<typename T, typename...args> void debug(T a, args...b) {
	cout << a << ", ";
	debug(b...);
}

template<typename T>
ostream& operator << (ostream &os, const vector<T> &a) {
	os << "[";
	int f = 0;
	for(auto &x : a) os << (f++ ? ", " : "") << x;
	os << "]";
	return os;
}

template<typename T>
ostream& operator << (ostream &os, const set<T> &a) {
	os << "{";
	int f = 0;
	for(auto &x : a) os << (f++ ? ", " : "") << x;
	os << "}";
	return os;
}

template<typename T>
ostream& operator << (ostream &os, const multiset<T> &a) {
	os << "{";
	int f = 0;
	for(auto &x : a) os << (f++ ? ", " : "") << x;
	os << "}";
	return os;
}

template<typename A, typename B>
ostream& operator << (ostream &os, const map<A, B> &a) {
	os << "{";
	int f = 0;
	for(auto &x : a) os << (f++ ? ", " : "") << x;
	os << "}";
	return os;
}

template<typename A, typename B>
ostream& operator << (ostream &os, const pair<A, B> &a) {
	os << "(" << a.first << ", " << a.second << ")";
	return os;
}

template<typename A, size_t N>
ostream& operator << (ostream &os, const array<A, N> &a) {
	os << "{";
	int f = 0; 
	for(int i = 0; i < N; i++) {
		os << (f++ ? ", " : "") << a[i];
	}
	os << "}";
	return os;
}