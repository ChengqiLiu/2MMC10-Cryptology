import math
from Crypto.Util.number import inverse
a=1
b=20
modn=41

def Double(x,y):
    '''Double the point (x,y) modulo modn. Result is (x_r,y_r).'''
    if math.isinf(x) or math.isinf(y):
        return math.inf,math.inf
    l=(3*x*x+a)*inverse(2*y,modn)%modn
    x_r=(l*l-2*x)%modn
    y_r=(l*(x-x_r)-y)%modn
    return x_r,y_r

def Add(x_p,y_p,x_q,y_q):
    '''Add two poins (x_p,y_p) and (x_q,y_q) modulp modn. Result is (x_r,y_r).'''
    if math.isinf(x_p):
        return x_q,y_q
    if math.isinf(x_q):
        return x_p,y_p
    if x_p==x_q and y_p!=y_q:
        return math.inf,math.inf
    if x_p==x_q and y_p==y_q:
        return Double(x_p,y_p)
    l=(y_p-y_q)*inverse((x_p-x_q),modn)%modn
    x_r=(l*l-x_p-x_q)%modn
    y_r=(l*(x_p-x_r)-y_p)%modn
    return x_r,y_r

def Power(x,y,exp):
    '''Calculate the exp times of the point (x,y) modulo modn. Result is (x_r,y_r).'''
    n=1
    x_r=x
    y_r=y
    while n<exp:
        x_r,y_r=Add(x_r,y_r,x,y)
        n=n+1
    return x_r,y_r

def BSGS(x_A,y_A,x,y,order):
    '''
    Solve DLP using Baby Step Giant Step. 
    Return log_{(x,y)}(x_A,y_A) modulo modn. 
    The point (x,y) has order "order".
    '''
    t=math.floor(math.sqrt(order))
    _xt,_yt=Power(x,y,t)
    k=1
    _xkt=_xt
    _ykt=_yt
    g_map=dict()
    g_map[(_xt,_yt)]=k
    while k<math.floor(order/t):
        k=k+1
        _xkt,_ykt=Add(_xkt,_ykt,_xt,_yt)
        g_map[(_xkt,_ykt)]=k
    i=1
    while i<t+1:
        _xi,_yi=Power(x,y,i)
        _xAi,_yAi=Add(x_A,y_A,_xi,_yi)
        if (_xAi,_yAi) in g_map:
            return (g_map[(_xAi,_yAi)]*t-i)%modn
        i=i+1
    return -1

def kP_plus_mPA(k,m,P_x,P_y,PA_x,PA_y):
    '''Calculate the point k*P+m*PA, result is (x,y).'''
    x1,y1=Power(P_x,P_y,k)
    x2,y2=Power(PA_x,PA_y,m)
    return Add(x1,y1,x2,y2)

def f(x,y,l):
    '''Step function f for Pollard rho algorithm. '''
    R0x,R0y=23,22
    R0a,R0b=23,13
    R1x,R1y=30,20
    R1a,R1b=19,11
    R2x,R2y=0,26
    R2a,R2b=2,41
    R3x,R3y=27,38
    R3a,R3b=25,37
    if x%4==0:
        x_r,y_r=Add(x,y,R0x,R0y)
        l[0]=l[0]+R0a
        l[1]=l[1]+R0b
        return x_r,y_r
    elif x%4==1:
        x_r,y_r=Add(x,y,R1x,R1y)
        l[0]=l[0]+R1a
        l[1]=l[1]+R1b
        return x_r,y_r
    elif x%4==2:
        x_r,y_r=Add(x,y,R2x,R2y)
        l[0]=l[0]+R2a
        l[1]=l[1]+R2b
        return x_r,y_r
    else:
        x_r,y_r=Add(x,y,R3x,R3y)
        l[0]=l[0]+R3a
        l[1]=l[1]+R3b
        return x_r,y_r
        
def f2(x,y,l):
    '''Step function f for Pollard rho algorithm. '''
    P_x,P_y=3,38
    PA_x,PA_y=25,34
    if x%3==0:
        x_r,y_r=Add(x,y,P_x,P_y)
        l[0]=l[0]+1
        return x_r,y_r
    elif x%3==1:
        x_r,y_r=Double(x,y)
        l[0]=2*l[0]
        l[1]=2*l[1]
        return x_r,y_r
    else:
        x_r,y_r=Add(x,y,PA_x,PA_y)
        l[1]=l[1]+1
        return x_r,y_r

def Pollard_rho(x_A,y_A,x,y,order):
    '''
    Solve DLP using Pollard rho algorithm. 
    Return log_{(x,y)}(x_A,y_A) modulo modn. 
    The point (x,y) has order "order".
    '''
    #(Wx,Wy)=w_l[0]*(x,y)+w_l[1]*(x_A,y_A)
    Wx,Wy=20,39
    W_l=[2,3]
    i=0
    Sx,Sy=Wx,Wy
    S_l=W_l.copy()
    Fx,Fy=Wx,Wy
    F_l=W_l.copy()
    print('S'+str(i)+":(",Sx,',',Sy,') ','F'+str(i)+":(",Fx,',',Fy,') ',sep='')
    while True:
        Sx,Sy=f2(Sx,Sy,S_l)
        Fx,Fy=f2(Fx,Fy,F_l)
        Fx,Fy=f2(Fx,Fy,F_l)
        i=i+1
        print('S'+str(i)+":(",Sx,',',Sy,') ','F'+str(i)+":(",Fx,',',Fy,') ',sep='')
        if Sx==Fx and Sy==Fy:
            break
    if math.gcd(F_l[0]-S_l[0],p)==1:
        return (F_l[1]-S_l[1])*inverse(F_l[0]-S_l[0],p)%p
    else:
        return -1

if __name__=="__main__":
    # p=53
    # P_x,P_y=3,38
    # PA_x,PA_y=25,34
    # res=Pollard_rho(PA_x,PA_y,P_x,P_y,p)
    # print("The answer is:",res)
    # x_r,y_r=Power(P_x,P_y,res)
    # if x_r==PA_x and y_r==PA_y:
    #     print("The result is verified to be correct.")
    # else:
    #     print("The result is verified to be wrong!")
    print(pow(3,6,31))