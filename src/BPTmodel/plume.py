import numpy as np
import gsw
from melt import melt_calc

def line_plume(z,X,experiment):

    """
    This code solves BPT assuming a LINE plume geometry. This function is called by BPTmodel.
    It is not intended to be used as a stand-alone function.

    See also BPTmodel

    Parameters
    ----------
    z : array_like
        depth
    X : array_like
        Variables for ode solver. X[0] = the plume radius (m), X[0] = the plume velocity (m s^-1),
        X[1] = the plume conservative temperature (C), and X[2] = the plume absolute salinity (g/kg).
    experiment : xarray dataset
        experiment parameters. includes g = gravity acceleration, rho0 = reference density,
        ci / cw = heat capacity of ice / water, lambda* = coefficients for linearized freezing temp equation,
        L = latent heat of fusion, alpha = entrainment coefficient, Cd = drag coefficient,
        gammaT/S = turbulent transfer coefficient for temp/salt, iceT = temperature of ice,
        ambient_temperature/salinity/density = properties of ambient water,
        u_horizontal = horizontal velocity that enhances melt rates

    Returns
    -------
    yDot : ndarray, shape(4,)
        Rate of change of variables. yDot[0]] = d/dz of plume radius (m),
        yDot[1] = d/dz of plume velocity (m s^-1), yDot[2] = d/dz of plume temperature (C),
        yDot[3] = d/dz of plume salinity (psu)

    Functions Called
    ----------------
    melt_calc(...)

    """

    yDot = np.zeros(4)

    # calculate melt rate with calc_melt()
    wmelt = np.sqrt(X[1]**2 + experiment.horizontal_velocity**2)

    iSort = np.argsort(experiment.depth.values)
    ambientDensity = np.interp(z,experiment.depth.values[iSort],experiment.ambient_density.values[iSort])
    ambientTemperature = np.interp(z,experiment.depth.values[iSort],experiment.ambient_temperature.values[iSort])
    ambientSalinity = np.interp(z,experiment.depth.values[iSort],experiment.ambient_salinity.values[iSort])

    [melt, Tb, Sb] = melt_calc(z, wmelt, X[2], X[3], lambda1=experiment.lambda1, lambda2=experiment.lambda2, lambda3=experiment.lambda3, gammaT=experiment.gammaT, gammaS=experiment.gammaS,\
                                Cd=experiment.Cd, iceT=experiment.ice_temperature, L=experiment.latent_heat_of_fusion, ci=experiment.ice_heat_capacity, cw=experiment.seawater_heat_capacity)
    
    # ODEs for radius, velocity, temperature, and salinity. From Jenkins 2011.  

    yDot[0]=2*experiment.alpha+2*melt/(X[1])-X[0]*experiment.g*(ambientDensity-gsw.density.sigma0(X[3],X[2]))/(X[1]*X[1]*experiment.reference_density)+experiment.Cd

    yDot[1]=-experiment.alpha*X[1]/X[0]-melt/(X[0])+experiment.g*(ambientDensity-gsw.density.sigma0(X[3],X[2]))/(X[1]*experiment.reference_density)-experiment.Cd*X[1]/(X[0])

    yDot[2]=experiment.alpha*(ambientTemperature-X[2])/X[0]+melt*(Tb-X[2])/(X[0]*X[1])-experiment.gammaT*np.sqrt(experiment.Cd)*(X[2]-Tb)/(X[0])*wmelt/X[1]

    yDot[3]=experiment.alpha*(ambientSalinity-X[3])/X[0]+melt*(Sb-X[3])/(X[0]*X[1])-experiment.gammaS*np.sqrt(experiment.Cd)*(X[3]-Sb)/(X[0])*wmelt/X[1]

    return yDot


def point_plume(z,X,experiment):

    """

    this code solves BPT assuming a POINT plume geometry.
    This function is called by BPTmodel. It is not intended to be used as a stand-alone function.

    See also BPTmodel

    Parameters
    ----------
    z : array_like
        depth
    X : array_like
        Variables for ode solver. X[0] = the plume radius (m), X[1] = the plume velocity (m s^-1),
        X[2] = the plume conservative temperature (C), and X[3] = the plume absolute salinity (g/kg).
    experiment : xarray dataset
        experiment parameters. includes g = gravity acceleration, rho0 = reference density,
        ci / cw = heat capacity of ice / water, lambda* = coefficients for linearized freezing temp equation,
        L = latent heat of fusion, alpha = entrainment coefficient, Cd = drag coefficient,
        gammaT/S = turbulent transfer coefficient for temp/salt, iceT = temperature of ice,
        ambient_temperature/salinity/density = properties of ambient water,
        u_horizontal = horizontal velocity that enhances melt rates

    Returns
    -------
    yDot : ndarray, shape(4,)
        Rate of change of variables. yDot(1) = d/dz of plume radius (m),
        yDot(2) = d/dz of plume velocity (m s^-1), yDot(3) = d/dz of plume temperature (C),
        yDot(4) = d/dz of plume salinity (psu)

    Functions Called
    ----------------
    melt_calc(...)  

    """

    yDot = np.zeros(4)

    # calculate melt rate with calc_melt()
    wmelt = np.sqrt(X[1]**2 + experiment.horizontal_velocity**2)

    iSort = np.argsort(experiment.depth.values)
    ambientDensity = np.interp(z,experiment.depth.values[iSort],experiment.ambient_density.values[iSort])
    ambientTemperature = np.interp(z,experiment.depth.values[iSort],experiment.ambient_temperature.values[iSort])
    ambientSalinity = np.interp(z,experiment.depth.values[iSort],experiment.ambient_salinity.values[iSort])

    [melt, Tb, Sb] = melt_calc(z, wmelt, X[2], X[3], lambda1=experiment.lambda1, lambda2=experiment.lambda2, lambda3=experiment.lambda3, gammaT=experiment.gammaT, gammaS=experiment.gammaS,\
                                Cd=experiment.Cd, iceT=experiment.ice_temperature, L=experiment.latent_heat_of_fusion, ci=experiment.ice_heat_capacity, cw=experiment.seawater_heat_capacity)

    # ODEs for radius, velocity, temperature, and salinity. From Cowton et al., 2015
    yDot[0]=2*experiment.alpha+4*melt/(np.pi*X[1])-X[0]*experiment.g*(ambientDensity-gsw.density.sigma0(X[3],X[2]))/(2*X[1]*X[1]*experiment.reference_density)+2*experiment.Cd/np.pi

    yDot[1]=-2*experiment.alpha*X[1]/X[0]-4*melt/(np.pi*X[0])+experiment.g*(ambientDensity-gsw.density.sigma0(X[3],X[2]))/(X[1]*experiment.reference_density)-4*experiment.Cd*X[1]/(np.pi*X[0])

    yDot[2]=2*experiment.alpha*(ambientTemperature-X[2])/X[0]+4*melt*(Tb-X[2])/(np.pi*X[0]*X[1])-4*experiment.gammaT*np.sqrt(experiment.Cd)*(X[2]-Tb)/(np.pi*X[0])*wmelt/X[1]

    yDot[3]=2*experiment.alpha*(ambientSalinity-X[3])/X[0]+4*melt*(Sb-X[3])/(np.pi*X[0]*X[1])-4*experiment.gammaS*np.sqrt(experiment.Cd)*(X[3]-Sb)/(np.pi*X[0])*wmelt/X[1]

    return yDot

def corner_plume(z,X,experiment):
    """
    QUARTER_PLUME

    This code solves BPT assuming a QUARTER CONICAL plume geometry with two faces in contact with ice.
    This function is called by BPTmodel. It is not intended to be used as a stand-alone function.

    See also BPTmodel

    Parameters
    ----------
    z : array_like
        depth
    X : array_like
        Variables for ode solver. X[0] = the plume radius (m), X[1] = the plume velocity (m s^-1),
        X[2] = the plume conservative temperature (C), and X[3] = the plume absolute salinity (g/kg).
    experiment : xarray dataset
        experiment parameters. includes g = gravity acceleration, rho0 = reference density,
        ci / cw = heat capacity of ice / water, lambda* = coefficients for linearized freezing temp equation,
        L = latent heat of fusion, alpha = entrainment coefficient, Cd = drag coefficient,
        gammaT/S = turbulent transfer coefficient for temp/salt, iceT = temperature of ice,
        ambient_temperature/salinity/density = properties of ambient water,
        u_horizontal = horizontal velocity that enhances melt rates

    Returns
    -------
    yDot : ndarray, shape(4,)
        Rate of change of variables. yDot(1) = d/dz of plume radius (m),
        yDot(2) = d/dz of plume velocity (m s^-1), yDot(3) = d/dz of plume temperature (C),
        yDot(4) = d/dz of plume salinity (psu)

    FUNCTIONS called:
        melt_calc(...)  
        
    """

    yDot = np.zeros(4)

    # calculate melt rate with calc_melt()
    wmelt = np.sqrt(X[1]**2 + experiment.horizontal_velocity**2)

    iSort = np.argsort(experiment.depth.values)
    ambientDensity = np.interp(z,experiment.depth.values[iSort],experiment.ambient_density.values[iSort])
    ambientTemperature = np.interp(z,experiment.depth.values[iSort],experiment.ambient_temperature.values[iSort])
    ambientSalinity = np.interp(z,experiment.depth.values[iSort],experiment.ambient_salinity.values[iSort])

    [melt, Tb, Sb] = melt_calc(z, wmelt, X[2], X[3], lambda1=experiment.lambda1, lambda2=experiment.lambda2, lambda3=experiment.lambda3, gammaT=experiment.gammaT, gammaS=experiment.gammaS,\
                                Cd=experiment.Cd, iceT=experiment.ice_temperature, L=experiment.latent_heat_of_fusion, ci=experiment.ice_heat_capacity, cw=experiment.sewater_heat_capacity)

    # ODEs for radius, velocity, temperature, and salinity. 
    yDot[0]=2*experiment.alpha+8*melt/(np.pi*X[1])-X[0]*experiment.g*(ambientDensity-gsw.density.sigma0(X[3],X[2]))/(2*X[1]*X[1]*experiment.reference_density)+4*experiment.Cd/np.pi

    yDot[1]=-2*experiment.alpha*X[1]/X[0]-8*melt/(np.pi*X[0])+experiment.g*(ambientDensity-gsw.density.sigma0(X[3],X[2]))/(X[1]*experiment.reference_density)-8*experiment.Cd*X[1]/(np.pi*X[0])

    yDot[2]=2*experiment.alpha*(ambientTemperature-X[2])/X[0]+8*melt*(Tb-X[2])/(np.pi*X[0]*X[1])-8*experiment.gammaT*np.sqrt(experiment.Cd)*(X[2]-Tb)/(np.pi*X[0])*wmelt/X[1]

    yDot[3]=2*experiment.alpha*(ambientSalinity-X[3])/X[0]+8*melt*(Sb-X[3])/(np.pi*X[0]*X[1])-8*experiment.gammaS*np.sqrt(experiment.Cd)*(X[3]-Sb)/(np.pi*X[0])*wmelt/X[1]

    return yDot
