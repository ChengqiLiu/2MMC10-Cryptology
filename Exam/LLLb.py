import math
# import numpy as np

def dot(al,bl):
    l=len(al)
    ret=0
    for i in range(l):
        ret=ret+al[i]*bl[i]
    return ret

def rou(x):
    '''rou x to nearest integer.'''
    return math.floor(x+0.5)

def gs(M,n):
    '''
    Input a n*n matrix M. Each row is a base.
    Return a n*n marix M_ret, with Gram–Schmidt process.
    '''
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
    '''
    Input a n*n matrix M. Each row is a base.
    Return a n*n marix M_ret, with lll lattice basis reduction.
    '''
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
    '''
    Return the function array arr from n*n Matrix M.
    The function is F(x)=arr[0]+arr[1]x+...+arr[n-1]x^n-1.
    '''
    arr=[]
    for i in range(n):
        arr.append(M[0][i]/pow(X,i))
    return arr

def gev(arr,n,x):
    ''' Return the value F(x).'''
    res=arr[n-1]
    for k in range(1,n):
        res=res*x+arr[n-k-1]
    return res

def sol(arr,n,sol):
    ''' sol F(x) by Newton's method. sol is estimated value of x.'''
    arr_d=[]
    for i in range(1,n):
        arr_d.append(i*arr[i])
    x=sol
    while abs(gev(arr,n,x))>0.5:
        x=x-gev(arr,n,x)/gev(arr_d,n-1,x)
    return rou(x)

if __name__=="__main__":    
    n=7684607040813031964568123727442263397500506224420545927139570285341
    a=3855587076697238083701498334674944
    X=pow(2,36)
    B=[[n,0,0,0],[a,X,0,0],[0,a*X,pow(X,2),0],[0,0,a*pow(X,2),pow(X,3)]]
    length=4
    matrix_l=lll(B,length)
    print("The matrix after lll reduction is:\n",matrix_l)
    
    F=gef(matrix_l,length,X)
    print("The function of a is:",F)

    r=sol(F,length,X)
    print("r is:",r)
    if n%(a+r)==0:
        print("The answer is correct.")
    else:
        print("The answer is wrong!")

    print("p is:",int(a+r))
    print("q is:",int(n/(a+r)))