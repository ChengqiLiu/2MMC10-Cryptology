x=0
while x<13:
    y=0
    while y<13: 
        if (x*x+y*y-1+5*x*x*y*y)%13==0 :
            print("("+str(x)+","+str(y)+")")
        y=y+1
    x=x+1
