Px=1000
Py=2
PAx=837670
PAy=538535
p=1000003
a=1
x=Px
y=Py
while not (x ==PAx and y ==PAy):
    x_2=(x*Py+y*Px)%p
    y_2=(y*Py-x*Px)%p
    a=a+1
    x=x_2
    y=y_2

print("a is: "+str(a))


