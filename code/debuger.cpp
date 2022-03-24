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

template<typename A, typename B, typename C>
ostream& operator << (ostream &os, const tuple<A, B, C> &a) {
	os << "(" << get<0>(a) << ", " << get<1>(a) << ", " << get<2>(a) << ")";
	return os;
}

template<typename A, typename B, typename C, typename D>
ostream& operator << (ostream &os, const tuple<A, B, C, D> &a) {
	os << "(" << get<0>(a) << ", " << get<1>(a) << ", " << get<2>(a) << ", " << get<3>(a) << ")";
	return os;
}