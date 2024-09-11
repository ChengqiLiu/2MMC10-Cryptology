import math
from Crypto.Util.number import inverse
a=1
b=20

def Double(x,y,modn):
    '''Double the point (x,y) modulo modn. Result is (x_r,y_r).'''
    if math.isinf(x) or math.isinf(y):
        return math.inf,math.inf
    l=(3*x*x+a)*inverse(2*y,modn)%modn
    x_r=(l*l-2*x)%modn
    y_r=(l*(x-x_r)-y)%modn
    return x_r,y_r

def Add(x_p,y_p,x_q,y_q,modn):
    '''Add two poins (x_p,y_p) and (x_q,y_q) modulp modn. Result is (x_r,y_r).'''
    if math.isinf(x_p):
        return x_q,y_q
    if math.isinf(x_q):
        return x_p,y_p
    if x_p==x_q and y_p!=y_q:
        return math.inf,math.inf
    if x_p==x_q and y_p==y_q:
        return Double(x_p,y_p,modn)
    l=(y_p-y_q)*inverse((x_p-x_q),modn)%modn
    x_r=(l*l-x_p-x_q)%modn
    y_r=(l*(x_p-x_r)-y_p)%modn
    return x_r,y_r

def Power(x,y,exp,modn):
    '''Calculate the exp times of the point (x,y) modulo modn. Result is (x_r,y_r).'''
    n=1
    x_r=x
    y_r=y
    while n<exp:
        x_r,y_r=Add(x_r,y_r,x,y,modn)
        n=n+1
    return x_r,y_r

def BSGS(x_A,y_A,x,y,order,modn):
    '''
    Solve DLP using Baby Step Giant Step. 
    Return log_{(x,y)}(x_A,x_A) modulo modn. 
    The point (x,y) has order "order".
    '''
    t=math.floor(math.sqrt(order))
    _xt,_yt=Power(x,y,t,modn)
    k=1
    _xkt=_xt
    _ykt=_yt
    g_map=dict()
    g_map[(_xt,_yt)]=k
    while k<math.floor(order/t):
        k=k+1
        _xkt,_ykt=Add(_xkt,_ykt,_xt,_yt,modn)
        g_map[(_xkt,_ykt)]=k
    i=1
    while i<t+1:
        _xi,_yi=Power(x,y,i,modn)
        _xAi,_yAi=Add(x_A,y_A,_xi,_yi,modn)
        if (_xAi,_yAi) in g_map:
            return (g_map[(_xAi,_yAi)]*t-i)%modn
        i=i+1
    return -1
        
if __name__=="__main__":
    modn=41
    p=53
    P_x=3
    P_y=38
    PA_x=25 
    PA_y=34
    print("The answer is:",BSGS(PA_x,PA_y,P_x,P_y,p,modn))