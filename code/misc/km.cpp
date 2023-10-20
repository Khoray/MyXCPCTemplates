#include<iostream>
#include<cstdio>
#include<cmath>
#include<cstring>
using namespace std;
typedef long long ll;
const ll Maxn=505;
const ll inf=1e18;
ll n,m,map[Maxn][Maxn],matched[Maxn];
ll slack[Maxn],pre[Maxn],ex[Maxn],ey[Maxn];//ex,ey顶标
bool visx[Maxn],visy[Maxn];
void match(ll u)
{	
    ll x,y=0,yy=0,delta;
    memset(pre,0,sizeof(pre));
    for(ll i=1;i<=n;i++)slack[i]=inf;
    matched[y]=u;
    while(1)
    {	
        x=matched[y];delta=inf;visy[y]=1;
        for(ll i=1;i<=n;i++)
        {	
            if(visy[i])continue;
            if(slack[i]>ex[x]+ey[i]-map[x][i])
            {	
                slack[i]=ex[x]+ey[i]-map[x][i];
                pre[i]=y;
            }
            if(slack[i]<delta){delta=slack[i];yy=i;}
        }
        for(ll i=0;i<=n;i++)
        {	
            if(visy[i])ex[matched[i]]-=delta,ey[i]+=delta;
            else slack[i]-=delta;
        }
        y=yy;
        if(matched[y]==-1)break;
    }
    while(y){matched[y]=matched[pre[y]];y=pre[y];}
}
ll KM()
{	
    memset(matched,-1,sizeof(matched));
    memset(ex,0,sizeof(ex));
    memset(ey,0,sizeof(ey));
    for(ll i=1;i<=n;i++)
    {	
        memset(visy,0,sizeof(visy));
        match(i);
    }
    ll res=0;
    for(ll i=1;i<=n;i++)
        if(matched[i]!=-1)res+=map[matched[i]][i];
    return res;
}
int main()
{	
    scanf("%lld%lld",&n,&m);
    for(int i=1;i<=n;i++)
        for(int j=1;j<=n;j++)
            map[i][j]=-inf;
    for(ll i=1;i<=m;i++)
    {	
        ll u,v,w;
        scanf("%lld%lld%lld",&u,&v,&w);
        map[u][v]=w;
    }
    printf("%lld\n",KM());
    for(ll i=1;i<=n;i++)
        printf("%lld ",matched[i]);
    printf("\n");
    return 0;
}