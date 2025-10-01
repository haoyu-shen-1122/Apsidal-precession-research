#Calculating semimajoraxis
def semimajoraxis(M1, M2, Ps):
    a = (c.G*(M1+M2)*Ps**2/(4*np.pi**2))**(1/3)
    return a
#Calculate the dimensionless P/U_rel, make sure to decompose when using astropy
def P_Urel(M1,M2,a,e):
    P_Urel = (3*c.G*(M1+M2)) / (c.c**2 * a * (1 - e**2))
    return P_Urel
#Calculate the dimensionless P/U_tide, make sure to decompose when using astropy
#When using it, calculate by P_Utide(star1) + P_Utide(star2) to get P_Utide for the whole system
def P_Utide(M1, M2, e, R, logk2, a):
    k2 = 10**logk2
    f_e = (1 + 1.5*e**2 + 0.125*e**4) / (1 - e**2)**5
    Ct = (R/a)**5 * (M2/M1) * 15 * f_e
    return k2*Ct
#Calculate the dimensionless P/U_rot, make sure to decompose when using astropy
#When using it, calculate by P_Urot(star1) + P_Urot(star2) to get P_Urot for the whole system
def P_Urot(e, P, k2, R, a, M1, M2, w_r):
    g = (1 - e**2)**(-2)
    w_k = 2*np.pi/P
    return k2 * (R/a)**5*(1+M2/M1)*g*(w_r/w_k)**2
def getU(P_Urel, P_Utide, P_Urot, P):
    num = P
    denom = P_Urel + P_Utide + P_Urot
    return num/denom