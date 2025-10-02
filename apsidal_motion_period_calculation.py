from astropy import units as u
from astropy import constants as c

def semimajoraxis(M1, M2, Ps):
    """Function of semimajor calculation

    Args:
        M1: The mass of one object.
        M2: The mass of other object.

    Returns:
        Semimajor axis of the system.

    """
    a = (c.G*(M1+M2)*Ps**2/(4*np.pi**2))**(1/3)
    return a

def P_Urel(M1,M2,a,e):
     """Function to calculate P/U_rel

    Args:
        M1: The mass of one object.
        M2: The mass of other object.
        a: semimajor axis.
        e: eccentricity

    Returns:
        P/U_rel, make sure to decompose when using astropy.

    """
    P_Urel = (3*c.G*(M1+M2)) / (c.c**2 * a * (1 - e**2))
    return P_Urel
#Calculate the dimensionless P/U_tide, make sure to decompose when using astropy
#When using it, calculate by P_Utide(star1) + P_Utide(star2) to get P_Utide for the whole system
def P_Utide(M1, M2, e, R, logk2, a):
    """Function to calculate one term of P/U_tide

    Args:
        M1: The mass of one object.
        M2: The mass of other object.
        e: eccentricity
        R: Radius of the object corresponding to M1

    Returns:
        One term of P/U_tide, make sure to decompose when using astropy.
        
    Example:
        P_Utide(star1_parameter) + P_Utide(star2_parameter)

    """
    k2 = 10**logk2
    f_e = (1 + 1.5*e**2 + 0.125*e**4) / (1 - e**2)**5
    Ct = (R/a)**5 * (M2/M1) * 15 * f_e
    return k2*Ct
#Calculate the dimensionless P/U_rot, make sure to decompose when using astropy
#When using it, calculate by P_Urot(star1) + P_Urot(star2) to get P_Urot for the whole system
def P_Urot(e, P, k2, R, a, M1, M2, w_r):
    """Function to calculate one term of P/U_rot

    Args:
        M1: The mass of one object.
        M2: The mass of other object.
        e: eccentricity
        R: Radius of the object corresponds to M1
        k2: love number of the object correspnds to M1, which indicate the shape of the star/planet
        w_r: self rotational frequency
        P: Sidereal period

    Returns:
        One term of P/U_rot, make sure to decompose when using astropy.
        
    Example:
        P_Urot(star1_parameter) + P_Urot(star2_parameter)

    """
    g = (1 - e**2)**(-2)
    w_k = 2*np.pi/P
    return k2 * (R/a)**5*(1+M2/M1)*g*(w_r/w_k)**2
def getU(P_Urel, P_Utide, P_Urot, P):
    """Function to calcualte apsidal precession period. Make sure every parameter is '.decompes()' after using astropy
    Args:
        P_Urel: The calculated P/U_rel. 
        P_Utide: The exact calculated P/U_tide, which is P_Utide(star1_parameter) + P_Utide(star2_parameter).
        P_Urot: The exact calculated P/U_rot, which is P_Urot(star1_parameter) + P_Urot(star2_parameter).
    Rturns:
        Apsidal precession

    """
    num = P
    denom = P_Urel + P_Utide + P_Urot
    return num/denom