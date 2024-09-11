import math
from Exam import *
A=82
p=1019
x=401
for i in range(1,10):
    x=div(pow(x*x-1,2),4*x*(x*x+A*x+1),p)
    print(i,x)
