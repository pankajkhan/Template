// MAIN TEMPLATE 

#include<bits/stdc++.h>
using namespace std;

#define ll long long
#define ull unsigned long long
#define N 1000006
#define mod 1000000007
#define boost ios_base::sync_with_stdio(false);cin.tie(0)
#define prec(n) fixed<<setprecision(n)

#define pii pair<int,int> 
#define pll pair<ll,ll> 
#define fi first
#define se second
#define pb push_back
#define vi vector<int> 
#define vll vector<ll> 

ll modulo(ll num){ return ((num%mod)+mod)%mod;} // for negative integer
ll power(ll b,ll e,ll MOD=mod){ll ans=1; while(e){if(e%2) ans=(ans*b)%MOD; b=(b*b)%MOD; e/=2;} return ans;}

int main()
{
    boost;
    
    
    return 0;
}

// Implementation of extended euclid's gcd
ll x,y,xx,yy ;
void egcd(ll a, ll b)
{
if(b==0)
{
x=1;
y=0;
return;
}
egcd(b,a%b);
xx=y;
yy=x-(a/b)*y;
x=xx;y=yy;
}

// Normal Sieve

int prime[N];
bool isprime[N];

void init()
{
	fill(isprime,isprime+N,1);
	int i,j;
	for(i=2;i*i<N;i++)
	{
		if(!isprime[i])  continue;
		
		for(j=i*i;j<N;j+=i)
		isprime[j]=0;
	}
	
	j=0;
	for(i=2;i<N;i++) if(isprime[i]) prime[j++] = i ; 
}


// Linear Sieve . Also useful for Factorization.
const int N = 10000000;
int lp[N+1];
vector<int> pr;
for (int i=2; i<=N; ++i) {
    if (lp[i] == 0) {
        lp[i] = i;
        pr.push_back (i);
    }
    for (int j=0; j<(int)pr.size() && pr[j]<=lp[i] && i*pr[j]<=N; ++j)
        lp[i * pr[j]] = pr[j];
}  

//******************************************************

#define N 400005
 
class segtree
{
    ll lo[N],hi[N],sum[N],delta[N];
    public:
    segtree(int n,int def)
    {
        init(1,0,n-1,def);
    }
    void init(int i,int a,int b,int def)
    {
        lo[i]=a;
        hi[i]=b;
        delta[i]=0;
        if(a==b)
        sum[i]=def;
        if(a==b)
        return ;
        int mid=(a+b)/2;
        init(2*i,a,mid,def);
        init(2*i+1,mid+1,b,def);
    }
    void prop(int i)
    {
        delta[2*i]=(delta[2*i]+delta[i])%mod;
        delta[2*i+1]=(delta[2*i+1]+delta[i])%mod;
        
 
        delta[i]=0;
    }
    void update(int i)
    {
        sum[i]=sum[2*i]%mod+delta[2*i]*(hi[2*i]-lo[2*i]+1)%mod+sum[2*i+1]%mod+delta[2*i+1]*(hi[2*i+1]-lo[2*i+1]+1)%mod;
    }
    void incre(int a,int b,int val)
    {
        increment(1,a,b,val);
    }
    int query(int a,int b)
    {
     return summation(1,a,b);
    }
    void increment(int i,int a,int b,int val)
    {
        if(lo[i]>b || hi[i]<a)
        return ;
        if(lo[i]>=a && hi[i]<=b)
        {
            delta[i]=(delta[i]+val)%mod;
            return;
        }
        prop(i);
        increment(2*i,a,b,val);
        increment(2*i+1,a,b,val);
        update(i);
    }
    
    int summation(int i,int a, int b)
    {
        if(lo[i]>b || hi[i]<a)
        return 0;
        if(lo[i]>=a && hi[i]<=b)
        return (sum[i]%mod+delta[i]*(hi[i]-lo[i]+1)%mod)%mod ;
        
        prop(i);
        ll sumL=summation(2*i,a,b);
        ll sumR=summation(2*i+1,a,b);
        
        update(i);
        
        return (sumL+sumR)%mod;
    }
};
 
segtree t1(n,0),t2(m,1); // Initialization 

//********************************************* FENWICK TREE 
//Warning : 1-based Indexing must be used everywhere.
int BIT[N]={},a[N],size ;
 
void update(int i, int val)
{
	i=i;
	while(i<=size)
	{
		BIT[i]+=val;
		i+=(i & -i);
	}
}
 
int sum(int i)
{
	int sum=0;
	
	while(i>0)
	{
		sum+=BIT[i];
		i-=(i & -i)  ;
	}
	
	return sum;
}

//Djikstra's Algorithm 

int dis[15],b,i;
fill(dis,dis+15,2*mod);

dis[1]=0;

priority_queue< pii , vector< pii > , greater<pii > > pq;  // careful about '>' spacing

pq.push({0,1});

while(!pq.empty())
{
    b=pq.top().second;
    pq.pop();
    
    for(i=0;i<adj[b].size();i++)
    {
        if(dis[ adj[b][i].first ] >  dis[b]+ adj[b][i].second )
        {
            dis[adj[b][i].first ]=dis[b]+adj[b][i].second ;
            pq.push(make_pair(dis[ adj[b][i].first ],adj[b][i].first));
        }
    }
}

//Prim's Algorithm. 

bool tv[3005]={};
int d[3005];
fill(d,d+3005,inf);
d[1]=0;
int sum=0;

priority_queue< pii , vector< pii > , greater< pii > > pq ;
pq.push({0,1});

while(!pq.empty())
{
    int dis=0;
    u=pq.top().second;
    dis=pq.top().first;
    pq.pop();
   
    if(tv[u])   // its already in tree
    continue;

    sum+=dis;
    tv[u]=1;
   
    for(i=0;i<adj[u].size();i++)
    {
       if( adj[u][i].second  < d[adj[u][i].first]  )
       {
           d[adj[u][i].first]=adj[u][i].second ;
           pq.push({adj[u][i].second,adj[u][i].first});
       }
       
    }
}

// Floyd Warshall's Algorithm Sample Solution

cin>>n>>m;
    
int W[n+1][n+1][n+1],i,j,u,v,k,w ;

for(i=1;i<=n;i++)
{
for(j=1;j<=n;j++)
W[i][j][0]=inf;
W[i][i][0]=0;
}
while(m--)
{
    cin>>u>>v>>w ;   W[u][v][0]=w;W[v][u][0]=w;
}
int maxx=0;

for(k=1;k<=n;k++)
{
    for(i=1;i<=n;i++)
    {
        for(j=1;j<=n;j++)
        {
            W[i][j][k]=min(W[i][j][k-1],W[i][k][k-1]+W[k][j][k-1]);
        }
    }
}

for(i=1;i<=n;i++)
for(j=1;j<=n;j++)
{
    maxx=max(maxx,W[i][j][n]);
}

    
// (ncr) mod m when m is prime. 

ll invf[N],f[2*N];

ll inv(ll num)
{
    return power(num,mod-2);
}
 
void pre()
{
    invf[1]=inv(1);
    for(ll i=2;i<N;i++)
    {
        invf[i]=(invf[i-1] * inv(i))%mod;
    }
    f[0] = f[1]=1;
    for(ll i=2;i<2*N;i++)
    f[i]=(f[i-1]*i)%mod;
   
}
 
ll ncr(ll n,ll r)
{
	if(r>n) return 0;
    if(r==0 || n==r)
        return 1;
    if(r==1)
        return n;
    return (f[n]*invf[r]%mod*invf[n-r]%mod)%mod;
}


// Lucas Theorem for ncr mod p when n and r are too large.. p is prime

int nCrModpLucas(int n, int r, int p)
{
   // Base case
   if (r==0)
      return 1;
  // Compute last digits of n and r in base p
   int ni = n%p, ri = r%p;
 
   return (nCrModpLucas(n/p, r/p, p) * // Last digits of n and r
           nCr(ni, ri, p)) % p;  // Remaining digits              
}

// Sample coin change code. 

ll m[1005][500];
ll coin_change(ll sum,ll ind)
{
    if(sum==0)
    return 1;
    
    if(a[ind]>sum)
    return 0;
    
    if(m[sum][ind]!=-1)
    return m[sum][ind];
    
    ll ways=0,tot=0;
    
    while(tot<=sum)
    {
        ways+=coin_change(sum-tot,ind+1);
        tot+=a[ind];
    }
    return m[sum][ind]=ways;
}

coin_change(money,0); // start with index 0 


