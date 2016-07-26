#!/usr/bin/env python
"""
This file contains a function computing the Bethe-Bloch formula.
"""

def get_I(Z, approx = 'power'):
    """
    Returns the mean excitation energy in Joules for a given atomic charge 
    number.
    
    Parameters
    ----------
    
    Z : int
        atomic charge number
    approx : str
        Approximation method to use when Z is not found in the lookup table.
        Options are:
        
            * '10eV' for I(Z) = 10 eV * Z
            * 'power' for I(Z) = 16 eV * Z^0.9
        
    Returns
    -------
    
    I : float
        excitation energy in Joule
    """
    from scipy.constants import eV
    mean_excitation_energies = {1: 15.6, 2: 43.0, 3: 37.4, 4: 61.1, 5: 63.2,
                           6: 79.1, 7: 89.3, 8: 101.5, 9: 120.3, 10: 126.8,
                           11: 149.4, 12: 157.0, 13: 163.4, 14: 172.3, 
                           15: 181.6, 16: 192.0, 17: 176.0, 18: 198.0, 
                           19: 242.5, 20: 188.3, 21: 205.2, 22: 225.5, 
                           23: 240.5, 24: 244.2, 25: 258.5, 26: 282.9, 
                           27: 285.7, 28: 302.6, 29: 317.2, 30: 316.5, 
                           31: 329.4, 32: 337.2, 33: 347.1, 34: 356.8, 
                           35: 368.6, 36: 362.6, 37: 392.0, 38: 396.2,
                           39: 401.0, 40: 373.4, 41: 406.2, 42: 420.4,
                           43: 429.3, 44: 436.6, 45: 436.1, 46: 454.4,
                           47: 475.3, 48: 460.5, 49: 478.4, 50: 496.7,
                           51: 490.8, 52: 499.8, 53: 510.2, 54: 517.5}
    try:
        I_mean = mean_excitation_energies[int(Z)]*eV
    except KeyError:
        if approx == '10eV':
            I_mean = 10*eV*float(Z)
        elif approx == 'power':
            I_mean = 16*eV*float(Z)**0.9
    return I_mean
    
def bethebloch(v, Z_proj, Z_target, A_target, rho, I = None, approx = 'power'):
    """
    Computes the dE/dx for given parameters according to the Bethe- Bloch
    formula.
    
    Parameters
    ----------
    
    v : float
        Projectile speed in m/s
    Z_proj : int
        Projectile atomic charge number
    Z_target : int
        Target material atomic charge number
    rho : float
        Target material density in kg/m^3
    A_target : float
        Target atomic mass number
    I : float
        Target mean ionization potential. Taken from hardcoded list if not
        passed as a parameter.
    approx : str
        Approximation method to use when Z is not found in the lookup table.
        Options are:
        
            * '10eV' for I(Z) = 10 eV * Z
            * 'power' for I(Z) = 16 eV * Z^0.9
    
    Returns
    -------
    
    dEdx : float
        dE/dx for the given parameters in J/m
    """
    from scipy.constants import pi, c, epsilon_0, e, m_u, m_e
    from numpy import log
    if I is None:
        I = get_I(Z_target, approx = approx)
    beta = v/c
    n = Z_target*rho/(A_target*m_u)
    dEdx = - 1.*(4*pi/(m_e*c**2))*(n*Z_proj**2/beta**2)*\
        ((e**2/(4*pi*epsilon_0))**2)*\
        (log(2*m_e*c**2*beta**2/(I*(1 - beta**2))) - beta**2)
    return dEdx
    