from Pollard_rho import kP_plus_mPA

if __name__=="__main__":
    p=53
    P_x=3
    P_y=38
    PA_x=25 
    PA_y=34
    print("W0: ",kP_plus_mPA(2,3,P_x,P_y,PA_x,PA_y))
    print("R0: ",kP_plus_mPA(23,13,P_x,P_y,PA_x,PA_y))
    print("R1: ",kP_plus_mPA(19,11,P_x,P_y,PA_x,PA_y))
    print("R2: ",kP_plus_mPA(2,41,P_x,P_y,PA_x,PA_y))
    print("R3: ",kP_plus_mPA(25,37,P_x,P_y,PA_x,PA_y))