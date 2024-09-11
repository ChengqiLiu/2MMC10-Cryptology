import math

def inverse(u, v):
    """The inverse of :data:`u` *mod* :data:`v`."""

    u3, v3 = u, v
    u1, v1 = 1, 0
    while v3 > 0:
        q = u3 // v3
        u1, v1 = v1, u1 - v1*q
        u3, v3 = v3, u3 - v3*q
    while u1<0:
        u1 = u1 + v
    return u1

def Add_Edw(x1,y1,x2,y2,d,modn):
    x_res=(x1*y2+x2*y1)*inverse(1+d*x1*y1*x2*y2,modn)%modn
    y_res=(y1*y2-x1*x2)*inverse(1-d*x1*y1*x2*y2,modn)%modn
    return x_res,y_res

def Mul_Edw(x,y,times,d,modn):
    x_res,y_res=0,1
    for i in range(times):
        x_res,y_res=Add_Edw(x_res,y_res,x,y,d,modn)
    return x_res,y_res

def Find_Order_Edw(x,y,d,modn):
    x_res,y_res=0,1
    order=0
    while True:
        x_res,y_res=Add_Edw(x_res,y_res,x,y,d,modn)
        order=order+1
        if x_res==0 and y_res==1:
            break
    return order

def Log_BF_Edw(x1,y1,x2,y2,d,modn):
    x,y=x1,y1
    n=1
    while True:
        x,y=Add_Edw(x1,y1,x,y,d,modn)
        n=n+1
        if x_res==x2 and y_res==y2:
            break
    return n

if __name__=="__main__":
    x,y=3,7
    d=11
    modn=17
    print(Log_BF_Edw(x,y,14,7,d,modn))