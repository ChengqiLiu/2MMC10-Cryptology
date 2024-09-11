import math
def ps(modn):
    for i in range(modn):
        print(i*i%modn,end=' ')
def pt(typ,arg,modn):
    if typ not in [0,1,2,3]:
        print("Error!")
        return
    for x in range(modn):
        for y in range(modn):
            sqx=x*x%modn
            sqy=y*y%modn
            if typ==0 and (sqx+sqy)%modn==(1+arg*sqx*sqy)%modn \
                or typ==1 and (arg[0]*sqx+sqy)%modn==(1+arg[1]*sqx*sqy)%modn\
                    or typ==2 and sqy==((sqx+arg[0])*x+arg[1])%modn\
                        or typ==3 and (arg[1]*sqy)%modn==x*(sqx+arg[0]*x+1)%modn:
                        print((x,y),end=' ')
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
def adde(x1,y1,x2,y2,d,modn):
    x_res=(x1*y2+x2*y1)*inv(1+d*x1*y1*x2*y2,modn)%modn
    y_res=(y1*y2-x1*x2)*inv(1-d*x1*y1*x2*y2,modn)%modn
    return x_res,y_res
def mule(x,y,times,d,modn):
    x_res,y_res=0,1
    for i in range(times):
        x_res,y_res=adde(x_res,y_res,x,y,d,modn)
    return x_res,y_res
def orde(x,y,d,modn):
    x_res,y_res=0,1
    order=0
    while True:
        x_res,y_res=adde(x_res,y_res,x,y,d,modn)
        order=order+1
        if x_res==0 and y_res==1:
            break
    return order
def loge(x1,y1,x2,y2,d,modn):
    x,y=x1,y1
    n=1
    while True:
        x,y=adde(x1,y1,x,y,d,modn)
        n=n+1
        if x_res==x2 and y_res==y2:
            break
    return n
def addm(x_p,y_p,x_q,y_q,A,B,modn):
    if x_p==-1:
        return x_q,y_q
    if x_q==-1:
        return x_p,y_p
    if x_p==x_q and y_p!=y_q:
        return -1,-1
    if x_p==x_q and y_p==y_q:
        return doum(x_p,y_p,A,B,modn)
    l=(y_p-y_q)*inv((x_p-x_q),modn)%modn
    x_r=(B*l*l%modn-A-x_p-x_q)%modn
    y_r=(l*(x_p-x_r)-y_p)%modn
    return x_r,y_r
def mulm(x,y,A,B,exp,modn):
    n=1
    x_r=x
    y_r=y
    while n<exp:
        x_r,y_r=addm(x_r,y_r,x,y,A,B,modn)
        n=n+1
    return x_r,y_r
def doum(x,y,A,B,modn):
    if x==-1 or y==-1:
        return -1,-1
    l=((3*x*x+2*A*x+1)%modn)*inv(2*B*y,modn)%modn
    x_r=(B*l*l%modn-A-2*x)%modn
    y_r=(l*(x-x_r)-y)%modn
    return x_r,y_r
def logm(x1,y1,x2,y2,A,B,modn):
    x,y=x1,y1
    n=1
    while True:
        x,y=addm(x1,y1,x,y,A,B,modn)
        n=n+1
        if x==x2 and y==y2:
            break
    return n
def prt(n):
    for k in range(2,int(math.sqrt(n))+1):
        if (n%k == 0):return False
    return True