
#region import statements

import numpy as np
from scipy import interpolate as interp
from scipy.integrate import solve_ivp
import gsw
import xarray as xr
from BPTmodel import plume
from BPTmodel.melt import melt_calc

#endregion

def BPTmodel(aD,aT,aS,GL,Q, Ttype = 'PT', Stype = 'SA', lat = 71.5, lon=-51.5, iceT = -10, W = 100, uh = 0,\
             ui = None, alpha = 0.1, Cd = 2.5e-3, gammaT = 0.022, gammaS = 0.00062,\
             intMethod = 'linear', extValue_T = None, extValue_S = None, type = 'line'):
    """
    1D model of buoyant plumes in the context of marine-terminating glaciers, solved in depth space

    Parameters
    ----------
    aD : 1 dimensional array of floats
        ambient ocean depth. Values should be negative and aproach 0 at the free surface
    aT : 1 dimensional array of floats
        ambient ocean temperature. potential temperature referenced to 0 dbar if Ttype = 'PT', conservative temperature if Ttype = 'CT'
    aS : 1 dimensional array of floats
        ambient ocean salinity. practical salinity units if Stype = 'PSU', absolute salinity if Stype = 'SA'
    GL : float
        Grounding line depth. should be negative
    Q : float
        subglacial discharge flux. Should be greater than zero. Use a small value (eg. 1e-10) for 'stacked' melt plumes.
    Ttype : string, optional
        specify what type of temperature is used in aT. 'PT' or 'CT'
    Stype : string, optional
        specify what type of salinity is used in aS. 'PSU' or 'SA'
    lat : float, optional
        latitude for calculating density with GSW. defaults to western coast of greenland
    lon : float, optional
        longitude for calculating density with GSW. defaults to western coast of greenland
    iceT : int, optional
        temperature of ice. Must be less than or equal to 0.
    W : int, optional
        discharge outlet with for line plumes (m). Must be greater than or equal to 1.
    uh: int, optional
        horizontal velocity to be used for calculating melt (m/s).
    ui: float, optional
        initial upwelling velocity. If None, will be calculated based on a balance of momentum and buoyancy.
    alpha: float, optional
        entrainment coefficient
    Cd: float, optional
        drag coefficient
    gammaT: float, optional
        thermal transfer coefficient
    gammaS: float, optional
        haline transfer coefficient
    intMethod: str, optional
        interpolation method for ambient profile (methods of interpolation for ??????????)
    extValue_T: float, optional
        value to fill ambient temperature profile with to extend to full water depth. If None or 'extrapolate', will extrapolate.
    extValue_S: float, optional
        value to fill ambient salinity profile with to extend to full water depth. If None or 'extrapolate', will extrapolate.
    type: str, optional
        sepcify what typ of plume model to run. 'line', 'point', 'stacked', or 'corner'. 'stacked' runs a series of line-plumes to simulate no discharge. 'corner' is experimental.

    Returns
    -------
    experiment: xarray dataset
        Contains inputs given to set up model run, constants used in model run, and resulting plume variables, including velocity, temperature, salinity, and density

    Functions Called
    ----------------
    line_plume(...), point_plume(...), or corner_plume(...)
    melt_calc(...)

    Notes
    -----
    For the most current version, see <a href="matlab: 
    web('https://github.com/BridgetOvall/BPTModel.git')">BPTmodel on GitHub</a>.

"""

    # supress integration failure warning that occurs for any plume that does not reach the surface
    #warning('off','MATLAB:ode45:IntegrationTolNotMet')

    #region handle input errors

    if sum(aD>0)>0:
        raise ValueError('aD invalid. depths muss be less than or equal to 0')
    
    if GL>=0:
        raise ValueError('GL invalid. must be less than 0')
    
    if Q<=0:
        raise ValueError('Q invalid. must be greater than 0')
    
    if iceT>0:
        raise ValueError('iceT invalid. must be less than or equal to 0')
    
    if W<1:
        raise ValueError('W invalid. must be greater than or equal to 1')

    # for stacked plumes, always use a unit width
    if type == 'stacked':
        if W!=1:
            print('Stacked plumes are always calculated per meter terminus width. Setting W=1')
        W=1
        
    # simple validations of input profiles
    # --- This is designed to catch the most obvious problems, but should not be used as a validation of adequate data ---
    validTS = np.logical_and(~np.isnan(aT),~np.isnan(aS))
    top=np.max(aD[validTS])
    bottom=np.min(aD[validTS])
    getinput=0
    if np.sum(~np.isnan(aT))<10 or np.sum(~np.isnan(aS))<10:
        raise ValueError('Ambient ocean profiles must contain more than 10 good data points')
    
    if extValue_T is None or extValue_T == 'extrapolate':
        print("Missing temperature data will be extrapolated per the methods set by the interpolation method. ")
        extValue_T = 'extrapolate'
    else:
        print("Missing temperature data will be filled with a value of "+extValue_T+". ")

    if extValue_S is None or extValue_S == 'extrapolate':
        print("Missing salinity data will be extrapolated per the methods set by the interpolation method. ")
        extValue_S = 'extrapolate'
    else:
        print("Missing salinity data will be filled with a value of "+extValue_S+". ")

    if top<-10:
        print("Ambient ocean profiles are missing data over the top "+ str(abs(top)) + " m. " + TintMsg + SintMsg +" \n")
        getinput=1

    if bottom-GL>5:
        print("Ambient ocean profiles are missing data over the bottom " + (bottom-GL) + " m (above the grounding line). " + TintMsg + SintMsg +" \n")
        getinput=1

    if np.sum(np.isnan(aT))/aT.size>.1 or np.sum(np.isnan(aS))/aS.size>.1:
        print("Ambient ocean profiles contain greater than 10% NaNs. \n")
        getinput=1

    if abs(np.mean(np.diff(aD)))>5:
        print("Ambient ocean profiles have a mean depth interval greater than 5 m\n")
        getinput=1

    if getinput==1:
        userinfo=input("Would you like to continue? (y/n)\n")
        if not userinfo:
            userinfo='n'
        if userinfo!='y':
            print("\nModel run aborted\n")
            return None

    #endregion

    #region get inputs into correct form for calculation and raise any remaining input errors if needed

    #calculate ambient pressure
    aP = gsw.conversions.p_from_z(aD,lat)

    # convert aT and aS to absolute salinity and conservative temperature
    if Stype == 'PSU':
        aS = gsw.conversions.SA_from_SP(aS,aP,lon,lat)
    elif Stype != 'SA':
        raise ValueError('Stype invalid. Must be \'PSU\' or \'SA\'.')
    
    if Ttype == 'PT':
        print('assuming potential temperature is referenced to 0 dbar')
        aT = gsw.conversions.CT_from_pt(aS,aT)
    elif Ttype != 'CT':
        raise ValueError('Ttype invalid. Must be \'PT\' or \'CT\'.')
    
    # convert discharge to m2/s for line plume or keep as m3/s for point plume
    match type:
        case 'line' | 'stacked':
            qsg = Q/W
        case 'point' | 'corner':
            qsg = Q
        case _:
            raise ValueError('type invalid. must be \'line\', \'point\', \'stacked\', or \'corner\'')
        

    # interpolate ambient profile to 0.1 m depth increments, extrapolate/fill to full depth range
    aD_good = aD[validTS]
    aT_good = aT[validTS]
    aS_good = aS[validTS]

    depth = np.round(np.arange(0,GL-.1,-.1),1)
    pressure = gsw.conversions.p_from_z(depth,lat)

    if intMethod in ['linear','nearest','next','previous']:
        Tinterpolator = interp.interp1d(aD_good,aT_good, kind=intMethod,fill_value=extValue_T)
        Sinterpolator = interp.interp1d(aD_good,aS_good, kind=intMethod,fill_value=extValue_S)
    else:
        if extValue_T == 'extrapolate':
            exT = None
        else:
            exT = False
        if extValue_S == 'extrapolate':
            exS = None
        else:
            exS = False
        if intMethod == 'cubic':
            Tinterpolator = interp.CubicSpline(aD_good,aT_good,extrapolate=exT)
            Sinterpolator = interp.CubicSpline(aD_good,aS_good,extrapolate=exS)
        elif intMethod == 'pchip':
            Tinterpolator = interp.PchipInterpolator(aD_good,aT_good,extrapolate=exT)
            Sinterpolator = interp.PchipInterpolator(aD_good,aS_good,extrapolate=exS)
        elif intMethod in ['akima','makima']:
            Tinterpolator = interp.Akima1DInterpolator(aD_good,aT_good,method=intMethod,extrapolate=exT)
            Sinterpolator = interp.Akima1DInterpolator(aD_good,aS_good,method=intMethod,extrapolate=exS)
        else:
            raise ValueError('intMethod invalid. must be \'linear\', \'nearest\', \'next\', \'previous\', \'cubic\', \'pchip\', \'akima\', or \'makima\'.')
    
    ambientTemp = Tinterpolator(depth)
    ambientSalt = Sinterpolator(depth)
    if not intMethod in ['linear','nearest','next','previous']:
        if extValue_T != 'extrapolate':
            ambientTemp = np.nan_to_num(ambientTemp,nan=extValue_T)
        if extValue_S != 'extrapolate':
            ambientSalt = np.nan_to_num(ambientSalt,nan=extValue_S)

    #calculate ambient potential density profile 
    ambientRho = gsw.density.sigma0(ambientSalt,ambientTemp) #potential density referenced to 0

    #endregion

    #region Define experiment structure with inputs and constants

    experiment = xr.Dataset(
        data_vars = dict(
            ambient_temperature = ('depth',ambientTemp,{'units':'deg C', 'notes':'Conservative Temperature'}),
            ambient_salinity = ('depth',ambientSalt, {'units':'g/kg','notes':'Absolute Salinity'}),
            ambient_density = ('depth',ambientRho, {'units':'kg/m^3','notes':'potential density referenced to 0 dbar'})
        ),
        coords = dict(
            depth = ('depth', depth, {'units':'meters'})
        ),
        attrs = dict(
            groundingLine_depth = GL,
            discharge_volumeFlux = Q,
            lat = lat,
            lon = lon,
            ice_temperature = iceT,
            outlet_width = W,
            horizontal_velocity = uh,
            initial_velocity = ui,
            alpha = alpha,
            Cd = Cd,
            gammaT = gammaT,
            gammaS = gammaS,
            intMethod = intMethod,
            type = type,
            g = 9.81, #gravitational acceleration (m s^-2)
            seawater_heat_capacity = 3974, #heat capacity of seawater (J kg^-1 C^-1)
            ice_heat_capacity = 2009.0, #heat capacity of ice  (J kg^-1 C^-1)
            latent_heat_of_fusion = 335000, #latent heat of fusion (J kg^-1)
            reference_density = 1028, #reference density (kg m^-3) 
            lambda1 = -0.0573, #variation of freezing point with salinity (C psu^-1)
            lambda2 = 0.0832, #freezing point offset (C)
            lambda3 = 0.000761, #variation of freezing point with depth (C m^-1)
            surface = 0, #depth of surface where integration stops
            attribute_units = ['distances, depths, and lengths (m)','temperatures (Conservative Temperature, C)','salinities (Absolute Salinity g/kg)','velocities (m s^{-1})','densities (kg m^{-3})',\
                                'volume_flux (m^3 s^{-1})', 'g (m s^{-2})', 'heat_capacity (J kg^-1 C^-1)', 'latent_heat (J kg^{-1})','lambda1 (variation of freezing point with salinity, C psu^{-1})',\
                                'lambda2 (freezing point offset, C)', 'lambda3 (variation of freezing point with depth, C m^{-1})']
        )
    )

    ## calculate N2 of ambient profile
    experiment = experiment.assign(ambient_buoyancy_frequency = ('depth', np.real(np.sqrt(-experiment.g/experiment.reference_density * np.gradient(experiment.ambient_density.values)/np.gradient(experiment.depth.values))), {'units' : 's^{-1}'}))  ################ Bridget had just diff of density, but depth is in .1 m intervals?
    # experiment = experiment.assign(ambient_buoyancy_frequency = ('depth', np.real(np.sqrt(-experiment.g/experiment.reference_density * np.gradient(experiment.ambient_density.values))), {'units' : 's^{-1}'})) 


    #endregion

    #region Initial conditions of plume at grounding line

    ti = gsw.CT_freezing(0,gsw.p_from_z(GL,lat),0) #initial plume starts at freezing point
    si = 0.0001 # initial plume salinity  (note: needs > 0 for integration to converge)

    # find ambient rho and plume rho at grounding line --> calculate g' (gp)
    iGL = np.where(np.abs(depth)==np.abs(GL))[0][0]
    ambientRho_GL = ambientRho[iGL]
    plumeRho_GL = gsw.density.sigma0(si,ti)
    gp = (ambientRho_GL-plumeRho_GL)*experiment.g/experiment.reference_density

    # initial velocity of plume
    match type:
        case 'line' | 'stacked':
            if ui is None:
                ui = (gp*qsg/alpha)**(1/3)
            bi = qsg / ui; #plume thickness in cross-terminus direction (m)
        case 'point':
            if ui is None:
                ui = 2/np.pi*(np.pi**2*gp/(8*alpha))**(2/5)*qsg**(1/5)
            bi = np.sqrt(2*qsg / (np.pi*ui)) #plume radius (m)
        case 'corner':
            if ui is None:
                ui = 4/np.pi*(np.pi**2*gp/(32*alpha))**(2/5)*qsg**(1/5)
            bi = np.sqrt(4*qsg / (np.pi*ui)) #plume radius (m)
    
    #endregion

    #region Solve coupled system of ODEs 

    # using Explicit Runge-Kutta method of order 5(4).
    #    z = depth
    #    X is m x 4 with columns of plume: width, velocity, T , S

    match type:
        case 'line':
            sol = solve_ivp(plume.line_plume,(GL,experiment.surface),[bi,ui,ti,si],args=(experiment,),method='RK45',max_step=1)
            z = sol.t
            X = sol.y
        case 'point':
            sol = solve_ivp(plume.point_plume,(GL,experiment.surface),[bi,ui,ti,si],args=(experiment,),method='RK45',max_step=1)
            z = sol.t
            X = sol.y
        case 'corner':
            sol = solve_ivp(plume.corner_plume,(GL,experiment.surface),[bi,ui,ti,si],args=(experiment,),method='RK45',max_step=1)
            z = sol.t
            X = sol.y
        case 'stacked':
            di=GL # initial depth of bottom plume is grounding line depth
            count=0
            maxH=[]
            az = []
            aX = []
            while di<0:
                count+=1
                tsol = solve_ivp(plume.line_plume,(di,experiment.surface),[bi,ui,ti,si],args=(experiment,),method='RK45',max_step=1)
                tz = tsol.t
                tX = tsol.y
                [tz,ia]=np.unique(tz,return_index=True) #occassionally, intervals are too small to be distinguishable
                tX=tX[:,ia]
                az=az + [tz]
                aX=aX + [tX]
                #recalculate initial conditions for next plume
                di=round(tz[-1]+.1,1)
                ti = gsw.CT_freezing(0,gsw.p_from_z(di,lat),0)
                gp=(ambientRho[np.abs(depth)==np.abs(di)]-gsw.density.sigma0(si,ti))*experiment.g/experiment.reference_density
                ui=(gp[0]*qsg/experiment.alpha)**(1/3)
                bi=qsg/ui
                # store max plume height for each plume
                maxH = maxH + [di-.1]
                # return some progress status to user
                print('Plume {0:2.0f} outflows at {1:4.1f} m\n '.format(count,di-.1))
                if di<0:
                    print('    Solving for next plume ... \n')
                else:
                    print('    Plume integration complete. Preparing final output...\n')

#endregion
    
    #region format results

    ## Interpolate to specified grid 
    if type == 'stacked':
        radius = np.empty(np.shape(depth))*np.nan
        w = np.empty(np.shape(depth))*np.nan
        temperature = np.empty(np.shape(depth))*np.nan
        salinity = np.empty(np.shape(depth))*np.nan
        density = np.empty(np.shape(depth))*np.nan
        nD=np.empty(count)*np.nan
        for i in np.arange(count):
            dind = np.where(np.logical_and(depth>az[i][0],depth<az[i][-1]))[0]
            radius[dind] = np.interp(depth[dind],az[i],aX[i][0,:])
            w[dind] = np.interp(depth[dind],az[i],aX[i][1,:])
            temperature[dind] = np.interp(depth[dind],az[i],aX[i][2,:])
            salinity[dind] = np.interp(depth[dind],az[i],aX[i][3,:])
            density[dind] = gsw.density.sigma0(salinity[dind],temperature[dind])
            # find neutral density for each plume
            nDi=np.where(density[dind]<=ambientRho[dind])[0][0]
            nD[i] = depth[dind][nDi]
            ################################################################# Bridget was doing something with sum in matlab - I'm not super sure why
    else:
        [z,ia]=np.unique(z,return_index=True) #occassionally, intervals are too small to be distinguishable
        X=X[:,ia]
        radius=np.interp(depth,z,X[0,:])
        w=np.interp(depth,z,X[1,:])
        temperature=np.interp(depth,z,X[2,:])
        salinity=np.interp(depth,z,X[3,:])
        density = gsw.density.sigma0(salinity,temperature)

    ## Compute additional variables from solution of ODEs

    # define a geometric factor to account for differing plume areas
    match type:
        case 'line' | 'stacked':
            geom = experiment.outlet_width
        case 'point':
            geom = np.pi*radius/2
        case 'corner':
            geom = np.pi*radius/4

    # area (m^2)
    area = geom*radius

    # volume flux (m^3/s)
    Fv = area*w

    # calculate total velocity that drives melting: upwelling + horizontal
    uTot=np.sqrt(w**2 + experiment.horizontal_velocity**2)

    # melt rate (m/s --> m/day)
    [melt_ms,_,_] = melt_calc(depth, uTot, temperature, salinity, lambda1=experiment.lambda1, lambda2=experiment.lambda2, lambda3=experiment.lambda3, gammaS=experiment.gammaS, gammaT=experiment.gammaT,\
                                Cd=experiment.Cd, iceT=experiment.ice_temperature, L=experiment.latent_heat_of_fusion, ci=experiment.ice_heat_capacity, cw=experiment.seawater_heat_capacity)
    melt=24*60*60*melt_ms

    #momentum flux (kg m/s^2)
    Fm=w*w*(density+1000)*area

    ## Put plume properties in output structure

    experiment = experiment.assign(dict(
        plume_radius = ('depth',radius,{'units':'m'}),
        plume_velocity = ('depth',w,{'units':'m/s'}),
        plume_temperature = ('depth',temperature,{'units':'deg C', 'notes':'Conservative Temperature'}),
        plume_salinity = ('depth',salinity, {'units':'g/kg','notes':'Absolute Salinity'}),
        plume_density = ('depth',density, {'units':'kg/m^3','notes':'potential density referenced to 0 dbar'}),
        plume_area = ('depth',area,{'units':'m^2'}),
        plume_volumeFlux = ('depth',Fv,{'units':'m^3/s'}),
        plume_momentumFlux = ('depth',Fm,{'units':'kg m/s^2'}),
        plume_meltRate = ('depth',melt,{'units':'m/day'})
        )
    )

    #find depth of maximum melt
    mi = np.argmax(experiment.plume_meltRate.values)
    experiment = experiment.assign_attrs(dict(plume_maximum_melt_depth=depth[mi]))

    match type:
        case 'line' | 'point' | 'corner':
            # find terminal level based on density
            tl = np.where(experiment.plume_density.values <= ambientRho)[0][0]

            #find where integration stops
            wIndex = np.where(np.isnan(experiment.plume_velocity.values))[0]

            #if integration did not stop, then plume reached surface
            if np.size(wIndex) == 0:
                wIndex = 0
            else:
                wIndex = wIndex[-1]

            experiment = experiment.assign_attrs(dict(
                plume_neutral_density_depth = depth[tl],
                plume_maximum_height = depth[wIndex]        #find maximum height (which is above terminal depth)
                )
            )
        
        case 'stacked':
            # the above values have already been calculated separately for each plume
            experiment = experiment.assign_attrs(dict(
                plume_neutral_density_depth = nD,
                plume_maximum_height = maxH
                )
            )

    ## Return initial conditions
    experiment = experiment.assign_attrs(dict(
        initial_temperature = ti,
        initial_salinity = si
        )
    )

    #endregion

    return experiment