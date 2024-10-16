import numpy as np

def melt_calc(z, u, T, S, cw=3974,ci=2009,L=335000,lambda1=-0.0573,lambda2=0.0832,\
              lambda3=0.000761,iceT=-10,Cd=2.5e-3,gammaT=0.022,gammaS=0.00062):

    """
    Solves for the submarine melt rate of ice given ambient water conditions

    3-Equation Melt Parameterization (e.g. Jenkins 2011) 
    This function is called by BPTmodel, line_plume, and point_plume
    Can also be used as a stand-alone estimate of submarine melt neglecting plume dynamics

    Returns ice melt rate (m/s) and temperature (C) and salinity (psu) at the ice-ocean boundary.

    See also BPTmodel, line_plume, point_plume

    Parameters
    ----------
    z : array_like
        the depth (m) coordinates at which melt is being calculated. should be <0, increasing towards the free surface.
    u : array_like or float
        velocity that drives melting (m/s), same length as z or a single float that is turned into a uniform profile.
    T : array_like
        Conservative temperature outside the BL (e.g. in plume) (C), same length as z
    S : array_like
        Absolute salinity outside the BL (e.g. in plume) (psu), same length as z
    cw : float, optional
        heat capacity of seawater (J kg^-1 C^-1)
    ci : float, optional
        heat capacity of ice (J kg^-1 C^-1)
    L : float, optional
        latent heat of fusion (J kg^-1)
    lambda1 : float, optional
        variation of freezing point with salinity (C psu^-1)
    lambda2 : float, optional
        freezing point offset (C)
    lambda3 : float, optional
        variation of freezing point with depth (C m^-1)
    iceT : float, optional
        ice temperature (C)
    Cd : float, optional
        drag coefficient
    gammaT : float, optional
        thermal transfer coefficient
    gammaS : float, optional
        haline transfer coefficient

    Returns
    -------
    melt : ndarray
        submarine melt rate. same size as z
    Tb : ndarray
        boundary layer temperature. same size as z
    Sb : ndarray
        boundary layer salinity. same size as z

    """

    if np.size(T) != np.size(S) or np.size(T) != np.size(z) or np.size(S) != np.size(z):
        raise ValueError('Vectors for ambient properties must all have the same dimesnions')
    
    if np.size(u)!=1:
        if np.size(u) != np.size(z):
            raise ValueError('Vector for along-ice velocity must be scalar or have same dimensions as depth profile')
        
    if np.sum(z>0)>0:
        raise ValueError('depths must be less than 0')

    
    # if single value for u is entered, make it a uniform profile
    if np.size(u)==1:
        u=np.ones(np.size(z))*u


    # define a, b, c to solve quadratic equation for Sb
    a = lambda1 * (gammaT*cw - gammaS*ci)

    b = gammaS*ci*(lambda1 * S - lambda2-lambda3*z + iceT - (L/ci)) - gammaT*cw*(T-lambda2-lambda3*z)

    c = gammaS * S * (ci * (lambda2 + lambda3*z-iceT)+L)

    # calculate Sb from a, b, c
    Sb = (1/(2*a)) * (-b-np.sqrt(b**2 - 4*a*c))
    # calculate Tb from linearized freezing temp. equation
    Tb = lambda1*Sb+lambda2+lambda3*z

    # calculate melt
    melt = gammaS * np.sqrt(Cd) * np.abs(u) * (S-Sb)/Sb

    return melt, Tb, Sb