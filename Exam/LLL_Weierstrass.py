import math
def inv(u, v):
    u3, v3 = u, v
    u1, v1 = 1, 0
    while v3 > 0:
        q = u3 // v3
        u1, v1 = v1, u1 - v1*q
        u3, v3 = v3, u3 - v3*q
    while u1<0:
        u1 = u1 + v
    return u1
def div(x,y,n):
    return (x%n)*inv(y,n)%n
def po3(x,y,n):
    if y<0:
        return po3(inv(x,n),-y,n)
    r=1
    for i in range(y):
        r=r*x%n
    return r
def dot(al,bl):
    l=len(al)
    ret=0
    for i in range(l):
        ret=ret+al[i]*bl[i]
    return ret
def rou(x):
    return math.floor(x+0.5)
def gs(M,n):
    M_ret=[]
    for i in range(n):
        vi=M[i].copy()
        if i==0:
            M_ret.append(vi)
            continue
        for j in range(i):
            v=M_ret[j]
            uij=dot(vi,v)/dot(v,v)
            for k in range(n):
                vi[k]=vi[k]-uij*v[k]
        M_ret.append(vi)
    return M_ret    
def lll(M,n,delta=3/4):
    M_gs=gs(M,n)
    k=1
    while k<n:
        for j in range(k-1,-1,-1):
            ukj=dot(M[k],M_gs[j])/dot(M_gs[j],M_gs[j])
            if abs(ukj)>0.5:
                for h in range(n):
                    M[k][h]=M[k][h]-rou(ukj)*M[j][h]
                M_gs=gs(M,n)
        u=dot(M[k],M_gs[k-1])/dot(M_gs[k-1],M_gs[k-1])
        if dot(M_gs[k],M_gs[k])>(delta-u*u)*dot(M_gs[k-1],M_gs[k-1]):
            k=k+1
        else:
            # Swap M[k] and M[k-1]
            arr=M[k].copy()
            M[k]=M[k-1]
            M[k-1]=arr
            M_gs=gs(M,n)
            k=max(k-1,1)
    return M
def gef(M,n,X):
    arr=[]
    for i in range(n):
        arr.append(M[0][i]/pow(X,i))
    return arr
def gev(arr,n,x):
    res=arr[n-1]
    for k in range(1,n):
        res=res*x+arr[n-k-1]
    return res
def sol(arr,n,sol):
    arr_d=[]
    for i in range(1,n):
        arr_d.append(i*arr[i])
    x=sol
    while abs(gev(arr,n,x))>0.5:
        x=x-gev(arr,n,x)/gev(arr_d,n-1,x)
    return rou(x)
def douw(x,y,a,b,modn):
    if math.isinf(x) or math.isinf(y):
        return math.inf,math.inf
    l=(3*x*x+a)*inv(2*y,modn)%modn
    x_r=(l*l-2*x)%modn
    y_r=(l*(x-x_r)-y)%modn
    return x_r,y_r
def addw(x_p,y_p,x_q,y_q,a,b,modn):
    if math.isinf(x_p):
        return x_q,y_q
    if math.isinf(x_q):
        return x_p,y_p
    if x_p==x_q and y_p!=y_q:
        return math.inf,math.inf
    if x_p==x_q and y_p==y_q:
        return douw(x_p,y_p,a,b,modn)
    l=(y_p-y_q)*inv((x_p-x_q),modn)%modn
    x_r=(l*l-x_p-x_q)%modn
    y_r=(l*(x_p-x_r)-y_p)%modn
    return x_r,y_r
def mulw(x,y,a,b,exp,modn):
    n=1
    x_r=x
    y_r=y
    while n<exp:
        x_r,y_r=addw(x_r,y_r,x,y,a,b,modn)
        n=n+1
    return x_r,y_r
def bsgsw(x,y,x_A,y_A,a,b,order,modn):
    t=math.floor(math.sqrt(order))
    _xt,_yt=mulw(x,y,a,b,t,modn)
    k=1
    _xkt=_xt
    _ykt=_yt
    g_map=dict()
    g_map[(_xt,_yt)]=k
    while k<math.floor(order/t):
        k=k+1
        _xkt,_ykt=addw(_xkt,_ykt,_xt,_yt,a,b,modn)
        g_map[(_xkt,_ykt)]=k
    i=1
    while i<t+1:
        _xi,_yi=mulw(x,y,a,b,i,modn)
        _xAi,_yAi=addw(x_A,y_A,_xi,_yi,a,b,modn)
        if (_xAi,_yAi) in g_map:
            return (g_map[(_xAi,_yAi)]*t-i)%modn
        i=i+1
    return -1