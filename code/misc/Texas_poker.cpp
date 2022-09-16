#include <bits/stdc++.h>
using namespace std;
struct card{
  char suit;
  int rank;
  card(){
  }
  bool operator <(card C){
    return rank < C.rank || rank == C.rank && suit < C.suit;
  }
};
istream& operator >>(istream& is, card& C){
  string S;
  is >> S;
  if (S[0] == 'A'){
    C.rank = 14;
  } else if (S[0] == 'K'){
    C.rank = 13;
  } else if (S[0] == 'Q'){
    C.rank = 12;
  } else if (S[0] == 'J'){
    C.rank = 11;
  } else if (S[0] == 'T'){
    C.rank = 10;
  } else {
    C.rank = S[0] - '0';
  }
  C.suit = S[1];
  return is;
}
vector<int> hand(vector<card> C){
  sort(C.begin(), C.end());
  set<char> suits;
  for (int i = 0; i < 5; i++){
    suits.insert(C[i].suit);
  }
  if (suits.size() == 1 && C[4].rank - C[0].rank == 4){
    if (C[4].rank == 14){
      return vector<int>{9};
    } else {
      return vector<int>{8, C[4].rank};
    }
  }
  if (suits.size() == 1 && C[3].rank == 5 && C[4].rank == 14){
    return vector<int>{8, 5};
  }
  if (C[0].rank == C[3].rank){
    return vector<int>{7, C[0].rank, C[4].rank};
  }
  if (C[1].rank == C[4].rank){
    return vector<int>{7, C[1].rank, C[0].rank};
  }
  if (C[0].rank == C[2].rank && C[3].rank == C[4].rank){
    return vector<int>{6, C[0].rank, C[3].rank};
  }
  if (C[2].rank == C[4].rank && C[0].rank == C[1].rank){
    return vector<int>{6, C[2].rank, C[0].rank};
  }
  if (suits.size() == 1){
    return vector<int>{5, C[4].rank, C[3].rank, C[2].rank, C[1].rank, C[0].rank};
  }
  if (C[1].rank - C[0].rank == 1 && C[2].rank - C[1].rank == 1 && C[3].rank - C[2].rank == 1 && C[4].rank - C[3].rank == 1){
    return vector<int>{4, C[4].rank};
  }
  if (C[0].rank == 2 && C[1].rank == 3 && C[2].rank == 4 && C[3].rank == 5 && C[4].rank == 14){
    return vector<int>{4, 5};
  }
  if (C[0].rank == C[2].rank){
    return vector<int>{3, C[0].rank, C[4].rank, C[3].rank};
  }
  if (C[1].rank == C[3].rank){
    return vector<int>{3, C[1].rank, C[4].rank, C[0].rank};
  }
  if (C[2].rank == C[4].rank){
    return vector<int>{3, C[2].rank, C[1].rank, C[0].rank};
  }
  if (C[0].rank == C[1].rank && C[2].rank == C[3].rank){
    return vector<int>{2, C[2].rank, C[0].rank, C[4].rank};
  }
  if (C[0].rank == C[1].rank && C[3].rank == C[4].rank){
    return vector<int>{2, C[3].rank, C[0].rank, C[2].rank};
  }
  if (C[1].rank == C[2].rank && C[3].rank == C[4].rank){
    return vector<int>{2, C[3].rank, C[1].rank, C[0].rank};
  }
  if (C[0].rank == C[1].rank){
    return vector<int>{1, C[0].rank, C[4].rank, C[3].rank, C[2].rank};
  }
  if (C[1].rank == C[2].rank){
    return vector<int>{1, C[1].rank, C[4].rank, C[3].rank, C[0].rank};
  }
  if (C[2].rank == C[3].rank){
    return vector<int>{1, C[2].rank, C[4].rank, C[1].rank, C[0].rank};
  }
  if (C[3].rank == C[4].rank){
    return vector<int>{1, C[3].rank, C[2].rank, C[1].rank, C[0].rank};
  }
  return vector<int>{0, C[4].rank, C[3].rank, C[2].rank, C[1].rank, C[0].rank};
}
int dfs(vector<vector<int>> &alice, vector<vector<int>> &bob, int a, int b, int p){
  if (p == 6){
    if (alice[a] > bob[b]){
      return 1;
    } else if (alice[a] < bob[b]){
      return -1;
    } else {
      return 0;
    }
  } else {
    int mx = -1;
    for (int i = 0; i < 6; i++){
      if ((a >> i & 1) == 0 && (b >> i & 1) == 0){
        if (p % 2 == 0){
          mx = max(mx, -dfs(alice, bob, a | (1 << i), b, p + 1));
        } else {
          mx = max(mx, -dfs(alice, bob, a, b | (1 << i), p + 1));
        }
      }
    }
    return mx;
  }
}
int main(){
  int T;
  cin >> T;
  for (int i = 0; i < T; i++){
    vector<card> a(2);
    cin >> a[0] >> a[1];
    vector<card> b(2);
    cin >> b[0] >> b[1];
    vector<card> c(6);
    cin >> c[0] >> c[1] >> c[2] >> c[3] >> c[4] >> c[5];
    vector<vector<int>> alice(1 << 6), bob(1 << 6);
    for (int j = 0; j < (1 << 6); j++){
      if (__builtin_popcount(j) == 3){
        vector<card> ha = {a[0], a[1]};
        vector<card> hb = {b[0], b[1]};
        for (int k = 0; k < 6; k++){
          if ((j >> k & 1) == 1){
            ha.push_back(c[k]);
            hb.push_back(c[k]);
          }
        }
        alice[j] = hand(ha);
        bob[j] = hand(hb);
      }
    }
    int ans = dfs(alice, bob, 0, 0, 0);
    if (ans == 1){
      cout << "Alice" << endl;
    } else if (ans == -1){
      cout << "Bob" << endl;
    } else {
      cout << "Draw" << endl;
    }
  }
}