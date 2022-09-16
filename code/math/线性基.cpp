struct linear_basis {
	vector<int> base, kth;
	int size, max_size, builded;
	linear_basis(int n) : base(n), size(0), max_size(n), builded(0) {}

	void insert(int x) {
		builded = 0;
		for(int i = max_size - 1; i >= 0; i--) {
			if((x >> i) & 1) {
				if(!base[i]) {
					base[i] = x, size++; 
					break;
				}
				else x ^= base[i];
			}
		}
	}

	int get_max() {
		int ret = 0;
		for(int i = max_size - 1; i >= 0; i--) {
			if(!((ret >> i) & 1) && base[i]) ret ^= base[i];
		}
		return ret;
	}

	bool can_eq(int x) {
		int now = 0;
		for(int i = max_size - 1; i >= 0; i--) {
			if(((now >> i) & 1) != ((x >> i) & 1)) {
				if(!base[i]) return false;
				else now ^= base[i];
			}
		}
		return true;
	}


	int get_kth(int k) {
        if(k >= 1ll << size) return -1;
		if(!builded) buildk();
		int ret = 0;
		for(int i = size - 1; ~i; i--) {
			if(k >> i & 1) {
				ret ^= kth[i];
			}
		}
		return ret;
	}

	int get_rank(int x) { // return the number of values less than x TODO
		int tmpsz = size, ret = 0, now = 0;
		for(int i = max_size - 1; i >= 0; i--) {
			if(base[i]) tmpsz--;
			if((x >> i) & 1) {
				if(!((now >> i) & 1)) {
					ret += tmpsz * tmpsz;
					if(!base[i]) break;
					else now ^= base[i];
				} else {
					if(base[i]) ret += tmpsz * tmpsz;
				}
			} else {
				if((now >> i) & 1) {
					if(!base[i]) break;
					else now ^= base[i];
				}
			}
		}
		return ret;
	}
private:
	void buildk() {
		builded = 1;
		kth.resize(size);
		int cnt = size;
		for(int i = max_size - 1; ~i; i--) {
			if(base[i]) {
				for(int j = i - 1; ~j; j--) {
					if(base[i] >> j & 1) {
						base[i] ^= base[j];
					}
				}
			}
		}
		for(int i = max_size - 1; ~i; i--) {
			if(base[i]) {
				kth[--cnt] = base[i];
			}
		}
	}
};