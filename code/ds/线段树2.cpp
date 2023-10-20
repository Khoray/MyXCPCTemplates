template<class S, S (*merge)(S, S), S (*s_id)(), class F, S (*mapping)(F, S), F (*composition)(F, F), F (*f_id)()>
struct lazy_segment_tree {
    int sz;
    vector<S> a;
    struct node {
        int l, r;
        S val;
        F tag;
        node *lc, *rc;
    };
    
    node *root;
    lazy_segment_tree(int n) : sz(n), a(n + 1), root(new node()) {}
    lazy_segment_tree(vector<S> &x) : sz(x.size() - 1), a(x), root(new node()) { build(root, 1, sz); }

    void build(node *now, int L, int R) {
        now->l = L, now->r = R, now->val = s_id(), now->tag = f_id();
        if(L == R) {
            now->val = a[L];
            return;
        }
        int mid = L + R >> 1;
        build(now->lc = new node(), L, mid);
        build(now->rc = new node(), mid + 1, R);
        now->val = merge(now->lc->val, now->rc->val);
    }
    
    void build(int L, int R) {
    	build(root, L, R);
	}

    S query(node *now, int L, int R) {
        if(now->l > R || now->r < L) return s_id();
        if(now->l >= L && now->r <= R) return now->val;
        push_down(now);
        return merge(query(now->lc, L, R), query(now->rc, L, R));
    }
    
    S query(int L, int R) {
    	return query(root, L, R);
	}

    void update(node *now, int pos, S val) {
        if(now->l == now->r && now->l == pos) {
            now->val = val;
            return;
        }
        int mid = now->l + now->r >> 1;
        push_down(now);
        update(pos <= mid ? now->lc : now->rc, pos, val);
        now->val = merge(now->lc->val, now->rc->val);
    }
    
    void update(int pos, S val) {
    	update(root, pos, val);
	}

    void push_down(node *now) {
        now->lc->val = mapping(now->tag, now->lc->val);
        now->rc->val = mapping(now->tag, now->rc->val);
        now->lc->tag = composition(now->tag, now->lc->tag);
        now->rc->tag = composition(now->tag, now->rc->tag);
        now->tag = f_id();
    }

    void update_range(node *now, int L, int R, F f) {
        if(now->l > R || now->r < L) return;
        if(now->l >= L && now->r <= R) {
            now->val = mapping(f, now->val);
            now->tag = composition(f, now->tag);
            return;
        }
        push_down(now);
        update_range(now->lc, L, R, f);
        update_range(now->rc, L, R, f);
        now->val = merge(now->lc->val, now->rc->val);
    }
    
    void update_range(int L, int R, F f) {
    	update_range(root, L, R, f);
	}

	template<class T>
    pair<int, S> find_r(node *now, int pos, S now_val, T check_val) {
        // 如果线段树的区间完全小于要查询的点
        if(now->r < pos) return {sz + 1, s_id()};
        // 如果线段树的区间完全大于要查询的点        
        if(now->l >= pos) {
            S all_val = merge(now_val, now->val);
            if(check_val(all_val)) return {sz + 1, all_val};
            if(now->l == now->r) return {now->l, all_val};
        }
        // 如果不满足条件，在这个区间内二分
        auto [lp, lval] = find_r(now->lc, pos, now_val, check_val);
        if(lp != sz + 1) {
            return {lp, lval};
        } else {
            return find_r(now->rc, pos, lval, check_val);
        }
    }
    
    template<class T>
    pair<int, S> find_r(int pos, S now_val, T check_val) {
    	return find_r(root, pos, now_val, check_val);
	}
};