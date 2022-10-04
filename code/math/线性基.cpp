struct linear_basis {
	vector<int> base, kth;
	int size, max_size, builded;
	/// @brief 构造线性基，向量长度是n
	/// @param n 
	linear_basis(int n) : base(n), size(0), max_size(n), builded(0) {}
	/// @brief 插入一个数
	/// @param x 
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
	/// @brief 获取最大值
	/// @return 最大值
	int get_max() {
		int ret = 0;
		for(int i = max_size - 1; i >= 0; i--) {
			if(!((ret >> i) & 1) && base[i]) ret ^= base[i];
		}
		return ret;
	}
	/// @brief 查询是否能等于x
	/// @param x 
	/// @return 查询是否能等于x
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

	/// @brief 询问第k大的值，insert后需要先buildk
	/// @param k 
	/// @return 第k大的值
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
	/// @brief 找到小于等于x的数的个数
	/// @param x 
	/// @return 小于等于x的数的个数
	int get_rank(int x) { 
		int tmpsz = size, ret = 0, now = 0;
		for(int i = max_size - 1; i >= 0; i--) {
			if(base[i]) tmpsz--;
			if((x >> i) & 1) {
				if(!((now >> i) & 1)) {
					ret += 1ll << tmpsz;
					if(!base[i]) return -1;
					else now ^= base[i];
				} else {
					if(base[i]) ret += 1ll << tmpsz;
				}
			} else {
				if((now >> i) & 1) {
					if(!base[i]) return -1;
					else now ^= base[i];
				}
			}
		}
		return ret;
	}
private:
	/// @brief 消成上阶梯矩阵
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