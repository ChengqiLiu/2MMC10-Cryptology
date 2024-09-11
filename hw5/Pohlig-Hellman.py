import math
from Crypto.Util.number import inverse

def BSGS(base,power,order,modn):
    '''
    Solve log_base(power) mod modn using BSGS algorithm.
    Return an integer, or -1 if error.
    '''
    t=math.floor(math.sqrt(order))
    base_map=dict()
    k=0
    base_t=1
    base_map[base_t]=k
    k=k+1
    while k<math.floor(order/t):
        base_t=base_t*pow(base,t,modn)%modn
        base_map[base_t]=k
        k=k+1
    for i in range(t+1):
        _n=power*pow(base,i,modn)%modn
        if _n in base_map:
            return (base_map[_n]*t-i)%order
    return -1

def CRT(a_l,m_l):
    '''
    Use Chinese Remainder Theorem to solve the euation system:
    x=a_l[0] mod m_l[0]
    x=a_l[1] mod m_l[1]
    ...
    x=a_l[l] mod m_l[l]
    Output (res,modn), where res is the minimum positive integer solution,
    modn is the multiple. General solution is res+k*modn.
    '''
    l=len(a_l)
    if l != len(m_l):
        return -1
    
    M=1
    for i in range(l):
        M=M*m_l[i]

    tM_l=[]
    for i in range(l):
        Mi=int(M/m_l[i])
        t=inverse(Mi,m_l[i])
        tM_l.append(t*Mi)

    res=0
    for i in range(l):
        res=res+a_l[i]*tM_l[i]
    return res%M,M

def Pohlig_Hellman(base,power,order,p_l,e_l,modn):
    '''
    Solve log_base(power) mod modn using Pohlig-Hellman algorithm. 
    Return an integer, or -1 if error.
    The order of base is order=p_l[0]^e_l[0]*p_l[1]^e_l[1]...*p_l[l]^e_l[l].

    '''
    l=len(p_l)
    if l != len(e_l):
        return -1

    a_l=[]
    for i in range(l):
        ai_l=[0 for j in range(e_l[i])]
        expj=order
        for j in range(e_l[i]):
            expj=int(expj/p_l[i])
            basej=pow(base,int(order/p_l[i]),modn)
            powerj=pow(power,expj,modn)
            cj=1
            for k in range(j):
                cj=cj*pow(base,expj*pow(p_l[i],k)*ai_l[k],modn)%modn
            powerj=powerj*inverse(cj,modn)%modn
            res=BSGS(basej,powerj,p_l[i],modn)
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
    return CRT(a_l,m_l)

    
if __name__=="__main__":
    p=1249
    p_l=[2,3,13]
    e_l=[5,1,1]
    g=7
    hb=1195
    ans=Pohlig_Hellman(g,hb,p-1,p_l,e_l,p)[0]
    print("The answer is",ans)
    if pow(g,ans,p)==hb:
        print("The answer is correct.")
    else:
        print("The answer is wrong!")
    