import math
import numpy as np

def Round(x):
    '''Round x to nearest integer.'''
    return math.floor(x+0.5)

def GramSchmidt(M,n):
    '''
    Input a n*n matrix M. Each row is a base.
    Return a n*n marix M_ret, with Gram–Schmidt process.
    '''
    M_ret=np.array([])
    for i in range(n):
        vi=M[i].copy()
        if i==0:
            M_ret=np.append(M_ret,vi)
            M_ret=np.resize(M_ret,(1,n))
            continue
        for j in range(i):
            v=M_ret[j]
            uij=np.dot(vi,v)/np.dot(v,v)
            vi=vi-uij*v
            vi=np.resize(vi,(1,n))
        M_ret=np.append(M_ret,vi,axis=0)
    return M_ret
        
def LLL(M,n,delta=3/4):
    '''
    Input a n*n matrix M. Each row is a base.
    Return a n*n marix M_ret, with LLL lattice basis reduction.
    '''
    M_gs=GramSchmidt(M,n)
    k=1
    while k<n:
        for j in range(k-1,-1,-1):
            ukj=np.dot(M[k],M_gs[j])/np.dot(M_gs[j],M_gs[j])
            if abs(ukj)>0.5:
                M[k]=M[k]-Round(ukj)*M[j]
                M_gs=GramSchmidt(M,n)
        u=np.dot(M[k],M_gs[k-1])/np.dot(M_gs[k-1],M_gs[k-1])
        if np.dot(M_gs[k],M_gs[k])>(delta-u*u)*np.dot(M_gs[k-1],M_gs[k-1]):
            k=k+1
        else:
            # Swap M[k] and M[k-1]
            arr=np.copy(M[k])
            M[k]=M[k-1]
            M[k-1]=arr
            M_gs=GramSchmidt(M,n)
            k=max(k-1,1)
    return M

def GetFunction(M,n,X):
    '''
    Return the function array arr from n*n Matrix M.
    The function is F(x)=arr[0]+arr[1]x+...+arr[n-1]x^n-1.
    '''
    arr=np.array([])
    for i in range(n):
        arr=np.append(arr,M[0][i]/pow(X,i))
    return arr

def GetValue(arr,n,x):
    ''' Return the value F(x).'''
    res=arr[n-1]
    for k in range(1,n):
        res=res*x+arr[n-k-1]
    return res

def Solve(arr,n,sol):
    ''' Solve F(x) by Newton's method. sol is estimated value of x.'''
    arr_d=np.array([])
    for i in range(1,n):
        arr_d=np.append(arr_d,i*arr[i])
    x=sol
    while abs(GetValue(arr,n,x))>0.5:
        x=x-GetValue(arr,n,x)/GetValue(arr_d,n-1,x)
    return Round(x)

if __name__=="__main__":    
    n=7684607040813031964568123727442263397500506224420545927139570285341
    a=3855587076697238083701498334674944
    X=pow(2,36)
    B=np.array([[n,0,0,0],[a,X,0,0],[0,a*X,pow(X,2),0],[0,0,a*pow(X,2),pow(X,3)]])
    length=4
    matrix_l=LLL(B,length)
    
    
    print("The matrix after LLL reduction is:\n",matrix_l)
    
    F=GetFunction(matrix_l,length,X)
    print("The function of a is:",F.tolist())

    r=Solve(F,length,X)
    print("r is:",r)
    if n%(a+r)==0:
        print("The answer is correct.")
    else:
        print("The answer is wrong!")

    print("p is:",int(a+r))
    print("q is:",int(n/(a+r)))