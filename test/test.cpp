 
#include<bits/stdc++.h>
using namespace std;
typedef long long LL;
 
LL pri[16] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53};
 
LL ans,n,MAXY;
    //求n以内因子最多的数,多个答案输出小的
void dfs(LL Y,LL dep,LL num,LL limit){
        //因子数 深度 当前的数
    if(dep==16)return;
    //if(Y>=n)return;
    if(Y>MAXY)MAXY=Y,ans=num;
    if(Y==MAXY)ans=min(ans,num);
 
    for(int i=1;i<=limit;i++){
        //因为后面的次数不会超过前面,所以这里添加上限就是一个极大的剪枝
        if(n/pri[dep]<num)break;//保证不超过n
        dfs(Y*(i+1),dep+1,num*=pri[dep],i);
    }
}
 
int main(){
    int t;cin>>t;
    while(t--){
        scanf("%lld",&n);
        ans=1e18;
        MAXY=0;
        dfs(1,0,1,60);
        printf("%lld %lld\n",ans,MAXY);
    }
}
 
 