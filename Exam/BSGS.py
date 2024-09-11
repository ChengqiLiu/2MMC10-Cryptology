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

def douw(x,y,a,b,modn):
    '''Double the point (x,y) modulo modn. Result is (x_r,y_r).'''
    if math.isinf(x) or math.isinf(y):
        return math.inf,math.inf
    l=(3*x*x+a)*inv(2*y,modn)%modn
    x_r=(l*l-2*x)%modn
    y_r=(l*(x-x_r)-y)%modn
    return x_r,y_r

def addw(x_p,y_p,x_q,y_q,a,b,modn):
    '''Add two poins (x_p,y_p) and (x_q,y_q) modulp modn. Result is (x_r,y_r).'''
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
    '''Calculate the exp times of the point (x,y) modulo modn. Result is (x_r,y_r).'''
    n=1
    x_r=x
    y_r=y
    while n<exp:
        x_r,y_r=addw(x_r,y_r,x,y,a,b,modn)
        n=n+1
    return x_r,y_r

def bsgsw(x,y,x_A,y_A,a,b,order,modn):
    '''
    Solve DLP using Baby Step Giant Step. 
    Return log_{(x,y)}(x_A,x_A) modulo modn. 
    The point (x,y) has order "order".
    '''
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
        
if __name__=="__main__":
    a=1
    b=20
    modn=41
    p=53
    P_x=3
    P_y=38
    PA_x=25 
    PA_y=34
    print("The answer is:",bsgsw(P_x,P_y,PA_x,PA_y,a,b,p,modn))