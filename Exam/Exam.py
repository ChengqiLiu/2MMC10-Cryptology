import math
def f(x,a,b,G,H,modn):
    if(x%3==0):
        x=x*G%modn
        a=a+1
    elif(x%3==1):
        x=x*H%modn
        b=b+1
    elif(x%3==2):
        x=x*x%modn
        a=2*a
        b=2*b
    return x,a,b
def frho(x,n):
    x=(x*x+23)%n
    return x
def gcd(*a):
    l=len(a)
    c,d=a[0],a[1]
    if(l==2):
        while d > 0:
            c, d = d, c % d
        return c
    for i in range(2,l):
        c,d=gcd(c,d),a[i]
    return gcd(c,d) 
def lcm(*a):
    l=len(a)
    c,d=a[0],a[1]
    if(l==2):
        return c * d // gcd(c, d)
    for i in range(2,l):
        c,d=lcm(c,d),a[i]
    return lcm(c,d)
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
def ord(x,n):
    num,r=0,1
    while True:
        r=(r*x)%n
        num=num+1
        if r==1:
            break
    return num  
def rot(x,sqn,n):
    a=1
    while(po3(a,sqn,n)!=x):
        a=a+1
    return a
def crt(a_l,m_l):
    l=len(a_l)
    if l != len(m_l):
        return -1
    M=1
    for i in range(l):
        M=M*m_l[i]
    tM_l=[]
    for i in range(l):
        Mi=int(M/m_l[i])
        t=inv(Mi,m_l[i])
        tM_l.append(t*Mi)
    res=0
    for i in range(l):
        res=res+a_l[i]*tM_l[i]
    return res%M,M
def bsgs(base,power,order,modn):
    t=math.floor(math.sqrt(order))
    base_map=dict()
    k=0
    base_t=1
    base_map[base_t]=k
    k=k+1
    while k<math.floor(order/t):
        base_t=base_t*po3(base,t,modn)%modn
        base_map[base_t]=k
        k=k+1
    for i in range(t+1):
        _n=power*po3(base,i,modn)%modn
        if _n in base_map:
            return (base_map[_n]*t-i)%order
    return -1   
def ph(base,power,order,p_l,e_l,modn):
    l=len(p_l)
    if l != len(e_l):
        return -1
    a_l=[]
    for i in range(l):
        ai_l=[0 for j in range(e_l[i])]
        expj=order
        for j in range(e_l[i]):
            expj=int(expj/p_l[i])
            basej=po3(base,int(order/p_l[i]),modn)
            powerj=po3(power,expj,modn)
            cj=1
            for k in range(j):
                cj=cj*po3(base,expj*pow(p_l[i],k)*ai_l[k],modn)%modn
            powerj=powerj*inv(cj,modn)%modn
            res=bsgs(basej,powerj,p_l[i],modn)
            ai_l[j]=res
            if i==0:
                print('a',j," is ",res,sep='')
            else:
                print("log_",basej,'(',powerj,')'," is ",res,sep='')
        res=0
        for j in range(e_l[i]):
            res=res*p_l[i]+ai_l[e_l[i]-1-j]
        a_l.append(res)
    m_l=[pow(p_l[i],e_l[i]) for i in range(l)]
    return crt(a_l,m_l)
def pr(s0,f0,a0,b0,G,H,modn):
    i=0
    _s,_as,_bs=s0,a0,b0
    _f,_af,_bf=f0,a0,b0
    start=True
    while(_s!=_f or start):
        _s,_as,_bs=f(_s,_as,_bs,G,H,modn)
        _f,_af,_bf=f(_f,_af,_bf,G,H,modn)
        _f,_af,_bf=f(_f,_af,_bf,G,H,modn)
        i=i+1
        if(i==1): start=False
        print(i,_s,_as,_bs,_f,_af,_bf)
def prho(n,x1):
    x=x1
    _x=frho(x,n)
    p=gcd((x-_x)%n,n)
    while(p==1):
        print(x,_x,end='; ')
        x=frho(x,n)
        _x=frho(_x,n)
        _x=frho(_x,n)
        p=gcd((x-_x)%n,n)
    if(p==n): return -1
    print(x,_x,end='; ')
    return p