import os
import numpy as np
import params
import scipy.io as sio
import math
import datetime
import pandas as pd


def CalculateSpatialDerivatives(bathy,x,y,sx,sy,sfac,pifac,ntopo,toporesi):
    h = np.zeros((sx,sy))
    hx = np.zeros((sx,sy))
    hy = np.zeros((sx,sy))
    bathy = np.maximum(bathy,-1)
    dhdx = sfac*(np.roll(bathy,-1,0)-np.roll(bathy,0,0))
    dhdy = sfac*(np.roll(bathy,-1,1)-np.roll(bathy,0,1))        
        
    for i in range(sx):
        plong = x[i]
        if plong>360:
            plong -= 360
        for j in range(sy):
            plat = y[j]
            #
            ib = np.floor(toporesi * plong).astype(int)
            ibp = ib+1
            if ibp >= ntopo:
                ibp = ibp-ntopo

            ibs = np.floor(toporesi * plong-0.5).astype(int)
            ibsp = ibs+1
            plongs = plong
            if ibs < -1:
                ibs = ntopo - 1
                plongs = plong + 360

            if ibsp >= ntopo:
                ibsp = ibsp-ntopo

            jb = np.floor(toporesi * (plat+90)).astype(int)
            jbs = np.floor(toporesi * (plat+90)-0.5).astype(int)
            b1 = bathy[ib,jb]
            b2 = bathy[ib,jb+1]
            b3 = bathy[ibp,jb]
            b4 = bathy[ibp,jb+1]
            dely = toporesi * (plat+90)-jb
            delx = toporesi * plong-ib
            d1 = (1 - delx) * (1 - dely)
            d2 = dely * (1 - delx)
            d3 = delx * (1 - dely)
            d4 = delx * dely
            h[i,j] = np.exp(d1*np.log(b1+11) + d2*np.log(b2+11) + d3*np.log(b3+11) + d4*np.log(b4+11)) - 11
                
            b1 = dhdx[ibs,jbs]
            b2 = dhdx[ibs,jbs+1]
            b3 = dhdx[ibsp,jbs]
            b4 = dhdx[ibsp,jbs+1]
            dely = -0.5 + toporesi * (plat+90) - jbs
            delx = -0.5 + toporesi * plongs - ibs
            d1 = (1 - delx) * (1 - dely)
            d2 = dely * (1 - delx)
            d3 = delx * (1 - dely)
            d4 = delx * dely
            hx[i,j] = (b1 * d1+b2 * d2+b3 * d3+b4 * d4) / np.cos(pifac * plat)
            
            b1 = dhdy[ibs,jbs]
            b2 = dhdy[ibs,jbs+1]
            b3 = dhdy[ibsp,jbs]
            b4 = dhdy[ibsp,jbs+1]
            hy[i,j] = (b1*d1 + b2*d2 + b3*d3 + b4*d4)
    return h, hx, hy


def EstimateTopographicHeight(bxmin,bxmax,bymin,bymax,dellatlong):
    mat = sio.loadmat('data/bathymetry_high.mat')
    bathy = mat['bathy']
    ntopo,_ = np.shape(bathy)
    topores = 360/ntopo
    toporesi = 1/topores
    sfac = 1./(topores*60.*1852)
    pifac = math.acos(-1)/180
    
    x = np.round(np.arange(bxmin,bxmax+1e-8,dellatlong),4)
    y = np.round(np.arange(bymin,bymax+1e-8,dellatlong),4)
    sx = np.max(np.shape(x))
    sy = np.max(np.shape(y))

    if sx == sy:
        sx += 1
        x = np.concatenate([x,[bxmax+dellatlong]])        
        
    h, hx, hy = CalculateSpatialDerivatives(bathy,x,y,sx,sy,sfac,pifac,ntopo,toporesi)
    
    return h,hx,hy,x,y
            

def EstimateDragCoefficients(plat,plong,sfac):
    
    pifac = math.acos(-1)/180         # pi number    
    
    # Load neutral drag coefficients
    mat = sio.loadmat('data/C_Drag500.mat')
    cd = mat['cd']
    cd = 0.9 * cd / (1+50*cd)       # This corrects the drag coefficient to be better applied to gradient wind
                                    # see Esau et al. (2004)
    cd = np.maximum(cd,1e-3)        # Align over-water values with Fortran program (including some wave drag effect)

    # Interpolate drag coefficient and its gradients to POI
    sy = np.max(np.shape(plat))
    sx = np.max(np.shape(plong))
    dcddx = sfac * ( np.roll(cd,-1,0) - cd )
    dcddy = sfac * ( np.roll(cd,-1,1) - cd )    
    
    cdrag = np.zeros((sx,sy))
    cdx = np.zeros((sx,sy))
    cdy = np.zeros((sx,sy))
    
    for i in range(sx):
        for j in range(sy):
            ib = np.floor(4*plong[i]).astype(int)
            if ib >= 1440:
                ib -= 1440
            
            ibp = ib + 1
            if ibp > 1440-1:
                ibp = 0
            
            ibs = np.floor(4*plong[i]-0.5).astype(int)
            plongs = plong
            if ibs < 0:
                ibs += 1440
                plongs = plong+360
            
            if ibs >= 1440-1:
                ibs -= 1440
                plongs = plong-360 
                
            ibsp = ibs+1
            jb = np.floor(4*(plat[j]+90)).astype(int)
            jbs = np.floor(4*(plat[j]+90)-0.5).astype(int)
            b1 = cd[ib,jb]
            b2 = cd[ib,jb+1]
            b3 = cd[ibp,jb]
            b4 = cd[ibp,jb+1]
            b1x = dcddx[ibs,jbs]
            b2x = dcddx[ibs,jbs+1]
            b3x = dcddx[ibsp,jbs]
            b4x = dcddx[ibsp,jbs+1]
            b1y = dcddy[ibs,jbs]
            b2y = dcddy[ibs,jbs+1]
            b3y = dcddy[ibsp,jbs]
            b4y = dcddy[ibsp,jbs+1]
            dely = 4*(plat[j]+90) - jb
            delx = 4*plong[i] - ib
            d1 = (1-delx) * (1-dely)
            d2 = dely * (1.-delx)
            d3 = delx * (1.-dely)
            d4 = delx * dely
            cdrag[i,j] = d1*b1+d2*b2+d3*b3+d4*b4
            dely = -0.5 + 4 * (plat[j]+90) - jbs
            delx = -0.5 + 4 * plongs[i] - ibs
            d1 = (1.-delx) * (1.-dely)
            d2 = dely * (1-delx)
            d3 = delx * (1-dely)
            d4 = delx * dely
            cdx[i,j] = (d1*b1x + d2*b2x + d3*b3x + d4*b4x) / np.cos(pifac *plat[j])
            cdy[i,j] = d1*b1y + d2*b2y + d3*b3y + d4*b4y
            
    return cdrag,cdx,cdy


def CalculateDistancePOTfromTrack(plat,plong,latstore,longstore,nn,m,sx,sy,ngrid,dfac):
    pifac = math.acos(-1)/180         # pi number        
    dx = np.zeros((nn,m,sx,sy))
    dy = np.zeros((nn,m,sx,sy))
    for i in range(sx):
        for jj in range(sy):
            if ngrid == 1:
                j=i
            else:
                j=jj
            
            if (np.ndim(longstore)==1) & (np.ndim(latstore)==1):
                dx[:,:,i,j] = dfac * np.cos(pifac*plat[j]) * (plong[i]-longstore)
                dy[:,:,i,j] = dfac * (plat[j]-latstore)
            elif (np.ndim(longstore)==2) & (np.ndim(latstore)==2):
                dx[:,:,i,j] = dfac * np.cos(pifac*plat[j]) * (plong[i]-longstore)
                dy[:,:,i,j] = dfac * (plat[j]-latstore)                
            elif (np.ndim(longstore)==4) & (np.ndim(latstore)==4):
                dx[:,:,i,j] = dfac * np.cos(pifac*plat[j]) * (plong[i]-longstore[:,:,i,jj])
                dy[:,:,i,j] = dfac * (plat[j]-latstore[:,:,i,jj])

    radius = np.sqrt(dx*dx+dy*dy)   
    return radius,dx,dy


def CalculateWindPrimary(utf, vtf, vf, rf, rmf, vsef, rmsef, latf, longf, w, 
                         nn, jf, sx, sy, h, hx, hy, cdrag, hf, hxf, 
                         hyf, cdf, cdx, cdxf, cdy, cdyf, Hi, Htrop, Omega, pifac, 
                         deltar, deltari, knotfac, latfac, timereswi, costhetaf, 
                         sinthetaf, thresM, wprofile, adj_water=False):
    V = np.zeros((nn,jf,sx,sy))
    Vd = np.zeros((nn,jf,sx,sy))
    Vpp = np.zeros((nn,jf,sx,sy))
    Vmp = np.zeros((nn,jf,sx,sy))
    Vpm = np.zeros((nn,jf,sx,sy))
    Vmm = np.zeros((nn,jf,sx,sy))
    Vrp = np.zeros((nn,jf,sx,sy))
    Vrm = np.zeros((nn,jf,sx,sy))
    up = np.zeros((nn,jf,sx,sy))
    um = np.zeros((nn,jf,sx,sy))
    
    V = windprofiles( vf, rmf, rf, wprofile)
    Vd = windprofilem( vf, rmf, vsef, rmsef, rf, wprofile)
    Vpp[:,1:jf-1,:,:] = windprofiles( vf[:,2:jf,:,:], rmf[:,2:jf,:,:], rf[:,1:jf-1,:,:]+deltar, wprofile)
    Vpm[:,1:jf-1,:,:] = windprofiles( vf[:,2:jf,:,:], rmf[:,2:jf,:,:], np.maximum(rf[:,1:jf-1,:,:]-deltar,0), wprofile)
    Vmp[:,1:jf-1,:,:] = windprofiles( vf[:,0:jf-2,:,:], rmf[:,0:jf-2,:,:], rf[:,1:jf-1,:,:]+deltar, wprofile)
    Vmm[:,1:jf-1,:,:] = windprofiles( vf[:,0:jf-2,:,:], rmf[:,0:jf-2,:,:], np.maximum(rf[:,1:jf-1,:,:]-deltar,0), wprofile)
    Vrp[:,:,:,:] = windprofiles( vf[:,:,:,:], rmf[:,:,:,:], rf[:,:,:,:]+deltar, wprofile)
    Vrm[:,:,:,:] = windprofiles( vf[:,:,:,:], rmf[:,:,:,:], np.maximum(rf[:,:,:,:]-deltar,0), wprofile)
    
    # Convert to meters per second. Done now because windprofile expects knots.
    V = knotfac * V
    Vd = knotfac * Vd
    Vpp = knotfac * Vpp
    Vpm = knotfac * Vpm
    Vmp = knotfac * Vmp
    Vmm = knotfac * Vmm
    Vrp = knotfac * Vrp
    Vrm = knotfac * Vrm
    
    vph = 0.5 * (V+Vrp)
    vmh = 0.5 * (V+Vrm)
    u1temp = vtf*costhetaf - utf*sinthetaf
    u2temp = vtf**2 + utf**2
    vnetp = np.sqrt(vph**2 + 2*vph*latfac*u1temp + u2temp)
    vnetm = np.sqrt(vmh**2 + 2*vmh*latfac*u1temp + u2temp)   
    
    for n in range(nn):
        for j in range(jf):
            hf[n,j,:,:] = h[:,:]
            hyf[n,j,:,:] = hy[:,:]
            hxf[n,j,:,:] = hx[:,:]
            cdf[n,j,:,:] = cdrag[:,:]
            cdyf[n,j,:,:] = cdy[:,:]
            cdxf[n,j,:,:] = cdx[:,:]

    cdfac = 1e3 * 0.5 * deltar * (cdxf*costhetaf + cdyf*sinthetaf)
    Cdp = cdf+cdfac
    Cdm = cdf-cdfac
    Cdp = np.maximum(Cdp,0)
    Cdp = np.minimum(Cdp,0.005)
    Cdm = np.maximum(Cdm,0)
    Cdm = np.minimum(Cdm,0.005)
    
    #  These lines of code added/modified  December 2018 to account for
    #  roughness changes over water
    if adj_water:
        facp = 1 + 0.0193*vnetp
        facm = 1 + 0.0193*vnetm
        facp[hf>0] = 1     #  Do not let roughness change wind wind over land
        facm[hf>0] = 1
        uekp = -Hi * Cdp * facp * vph * vnetp
        uekm = -Hi * Cdm * facm * vmh * vnetm
    else: 
        uekp = -Hi * Cdp * vph * vnetp
        uekm = -Hi * Cdm * vmh * vnetm
            
    cp = 1000 * 2 * Omega * np.sin(pifac*abs(latf))
    dMdrp = cp * (rf+0.5*deltar) + (rf+0.5*deltar)*deltari*(Vrp-V) + 0.5*(Vrp+V)
    dMdrp = np.maximum(dMdrp,thresM)
    dMdrm = cp * (rf-0.5*deltar) + (rf-0.5*deltar)*deltari*(V-Vrm) + 0.5*(Vrm+V)
    dMdrm = np.maximum(dMdrm,thresM)
    efacp = np.minimum((-1+2*((rf[:,1:jf-1,:,:]+deltar)/rmf[:,1:jf-1,:,:])**2),1)
    efacm = np.minimum((-1+2*((rf[:,1:jf-1,:,:]-deltar)/rmf[:,1:jf-1,:,:])**2),1)
    
    up[:,1:jf-1,:,:] = (rf[:,1:jf-1,:,:]+deltar)*(-0.5*timereswi*efacp*(Vpp[:,1:jf-1,:,:]-Vmp[:,1:jf-1,:,:])+uekp[:,1:jf-1,:,:])/dMdrp[:,1:jf-1,:,:]
    um[:,1:jf-1,:,:] = (rf[:,1:jf-1,:,:]-deltar)*(-0.5*timereswi*efacm*(Vpm[:,1:jf-1,:,:]-Vmm[:,1:jf-1,:,:])+uekm[:,1:jf-1,:,:])/dMdrm[:,1:jf-1,:,:]
    w[:,1:jf-1,:,:] = -Htrop*deltari*((rf[:,1:jf-1,:,:]+deltar)*up[:,1:jf-1,:,:]-(rf[:,1:jf-1,:,:]-deltar)*um[:,1:jf-1,:,:])/np.maximum(rf[:,1:jf-1,:,:],1)    
    return w, cp, V, Vd, Vrp, Vrm, u1temp, u2temp, Cdp, Cdm


def CalculateWindSecondary(utf, vtf, vf, rf, rmf, vsef, rmsef, latf, longf, w, nn, jf, 
                           sx, sy,h, hx, hy, cdrag, hf, hxf, hyf, cdf, cdx, cdxf, cdy, cdyf,
                           Hi, Htrop, Omega, pifac, deltar, deltari, knotfac, latfac, 
                           timereswi, costhetaf, sinthetaf, u1temp, u2temp, Cdp, Cdm, wprofile):
    V = np.zeros((nn,jf,sx,sy))
    Vd = np.zeros((nn,jf,sx,sy))
    Vpp = np.zeros((nn,jf,sx,sy))
    Vmp = np.zeros((nn,jf,sx,sy))
    Vpm = np.zeros((nn,jf,sx,sy))
    Vmm = np.zeros((nn,jf,sx,sy))
    Vrp = np.zeros((nn,jf,sx,sy))
    Vrm = np.zeros((nn,jf,sx,sy))
    up = np.zeros((nn,jf,sx,sy))
    um = np.zeros((nn,jf,sx,sy))
    V =  windprofiles( vsef, rmsef, rf, wprofile)
    Vpp[:,1:jf-1,:,:] = windprofiles(vsef[:,2:jf,:,:], rmsef[:,2:jf,:,:], rf[:,1:jf-1,:,:]+deltar, wprofile)
    Vpm[:,1:jf-1,:,:] = windprofiles(vsef[:,2:jf,:,:], rmsef[:,2:jf,:,:], np.maximum(rf[:,1:jf-1,:,:]-deltar,0), wprofile)
    Vmp[:,1:jf-1,:,:] = windprofiles(vsef[:,0:jf-2,:,:], rmsef[:,0:jf-2,:,:], rf[:,1:jf-1,:,:]+deltar, wprofile)
    Vmm[:,1:jf-1,:,:] = windprofiles(vsef[:,0:jf-2,:,:], rmsef[:,0:jf-2,:,:], np.maximum(rf[:,1:jf-1,:,:]-deltar,0), wprofile)
    Vrp = windprofiles(vsef, rmsef, rf+deltar, wprofile)
    Vrm = windprofiles(vsef, rmsef, np.maximum(rf-deltar,0), wprofile) 
    
    # Convert to meters per second. Done now because windprofile expects knots.
    V = knotfac * V
    Vd = knotfac * Vd
    Vpp = knotfac * Vpp
    Vpm = knotfac * Vpm
    Vmp = knotfac * Vmp
    Vmm = knotfac * Vmm
    Vrp = knotfac * Vrp
    Vrm = knotfac * Vrm
    
    vph = 0.5*(V+Vrp)
    vmh = 0.5*(V+Vrm)
    vnetp = np.sqrt(vph**2+2*vph*latfac*u1temp+u2temp)
    vnetm = np.sqrt(vmh**2+2*vmh*latfac*u1temp+u2temp)
    uekp = -Hi*Cdp*vph*vnetp
    uekm = -Hi*Cdm*vmh*vnetm
    
    cp = 1000*2*Omega*np.sin(pifac*abs(latf))
    dMdrp = cp*(rf+0.5*deltar)+(rf+0.5*deltar)*deltari*(Vrp-V)+0.5*(Vrp+V)
    dMdrp = np.maximum(dMdrp,2)
    dMdrm = cp*(rf-0.5*deltar)+(rf-0.5*deltar)*deltari*(V-Vrm)+0.5*(Vrm+V)
    dMdrm = np.maximum(dMdrm,2)
    efacp = np.minimum((-1+2*((rf[:,1,:,:]+deltar)/rmsef[:,1,:,:])**2),1)
    efacm = np.minimum((-1+2*((rf[:,1,:,:]-deltar)/rmsef[:,1,:,:])**2),1)
    
    #  Do not consider time rate of change when either velocity is zero
    efacp = efacp*np.minimum(Vpp[:,1:jf-1,:,:],1)*np.minimum(Vmp[:,1:jf-1,:,:],1)
    efacm = efacm*np.minimum(Vpm[:,1:jf-1,:,:],1)*np.minimum(Vmm[:,1:jf-1,:,:],1)
    up[:,1:jf-1,:,:] = (rf[:,1:jf-1,:,:]+deltar)*(-0.5*timereswi*efacp*(Vpp[:,1:jf-1,:,:]-Vmp[:,1:jf-1,:,:])+uekp[:,1:jf-1,:,:])/dMdrp[:,1:jf-1,:,:]
    um[:,1:jf-1,:,:] = (rf[:,1:jf-1,:,:]-deltar)*(-0.5*timereswi*efacm*(Vpm[:,1:jf-1,:,:]-Vmm[:,1:jf-1,:,:])+uekm[:,1:jf-1,:,:])/dMdrm[:,1:jf-1,:,:]
    
    #  Do not calculate vertical velocities if either radial velocity is zero.
    #  This is necessary because secondary wind maximum can vanish from one time 
    #  step to the next, so that time rate of change blows up
    ufac = np.minimum(np.abs(30*up[:,1:jf-1,:,:]),1)*np.minimum(np.abs(30*um[:,1:jf-1,:,:]),1)
    w[:,1:jf-1,:,:] = -Htrop*ufac*deltari*((rf[:,1:jf-1,:,:]+deltar)*up[:,1:jf-1,:,:]-(rf[:,1:jf-1,:,:]-deltar)*um[:,1:jf-1,:,:])/np.maximum(rf[:,1:jf-1,:,:],1)
        
    return w, cp, V, Vd, Vrp, Vrm,


def rainfieldx(nt,monthplot,dayplot,hourplot):
    """
    This script calculates the distribution of surface rain rate, in mm/hr, 
    for a given storm at a specified time.
    nt is the track number of the storm, monthplot is the month (0-12), 
    dayplot is the day, and hourplot is the GMT time.
    x and y are vectors containing the longitudes and latitudes of the grid
    rainrate is the surface rainfall rate, in mm/hr, at each grid point.
    """    
    mat = sio.loadmat('data/temp.mat')
    vstore = mat['vstore']
    vsestore = mat['vsestore']
    latstore = mat['latstore']
    longstore = mat['longstore']
    rmstore = mat['rmstore']
    rmsestore = mat['rmsestore']
    ut, vt, _ = utrans(latstore, longstore)
    
    # T600store = mat['T600store']
    u850store = mat['u850store']
    v850store = mat['v850store']
    monthstore = mat['monthstore']
    daystore = mat['daystore']
    hourstore = mat['hourstore']
    
    magfac = params.magfac
    randfac = params.randfac
    bound = params.bound
    deltax = params.deltax
    deltay = params.deltay
    bxmin = params.bxmin
    bxmax = params.bxmax
    bymin = params.bymin
    bymax = params.bymax
    dellatlong = params.dellatlong
    q900 = params.q900
    eprecip = params.eprecip
    wrad = params.wrad
    
    m_to_mm = 1000
    rowa_over_rowl = 0.00117
       
    nrm, mrm = np.shape(rmstore)
    rfac = magfac * (1+np.zeros((nrm,mrm)))
    pifac = math.acos(-1)/180
    knotfac = 1852./3600    
    jmaxd = np.argmin(vstore, axis=1)   # Get the length of each event by finding the first 0 element

    dum = np.zeros((jmaxd[nt]))
    V = np.column_stack((dum+2000, monthstore[nt,0:jmaxd[nt]].T, daystore[nt,0:jmaxd[nt]].T, hourstore[nt,0:jmaxd[nt]].T, dum, dum))
    timev = pd.to_datetime(dict(year=V[:,0], month=V[:,1], day=V[:,2], hour=V[:,3], minute=V[:,4], second=V[:,5]))
    timeplot = datetime.datetime(2000, monthplot, dayplot, hourplot, 0, 0)
    diff0 = timeplot - timeplot
    timediff = timeplot - timev
    timediff[timediff<diff0] = diff0
    jplot = np.argmin(timediff)
    jplot = np.maximum(jplot,1)
    jplot = np.minimum(jplot,398)

    latstorm = latstore[nt,jplot-1:jplot+2]
    longstorm = longstore[nt,jplot-1:jplot+2]
    vstorm = vstore[nt,jplot-1:jplot+2]
    rmstorm = rfac[nt,jplot-1:jplot+2] * [rmstore[nt,jplot-1:jplot+2]]
    vsestorm = vsestore[nt,jplot-1:jplot+2]
    rmsestorm = rfac[nt,jplot-1:jplot+2] * [rmsestore[nt,jplot-1:jplot+2]]
    utstorm = ut[nt,jplot-1:jplot+2]
    vtstorm = vt[nt,jplot-1:jplot+2]
    ush = np.zeros((3))
    vsh = np.zeros((3))

    if 'u850store' in locals():
        vdrift = 1.5*3600/1852
        vdrift = vdrift*latstore[1,1]/(abs(latstore[1,1])+1e-8)
        u850storm = u850store[nt,jplot-1:jplot+2]
        v850storm = v850store[nt,jplot-1:jplot+2]
        ush = 5 * knotfac * (utstorm-u850storm)
        vsh = 5 * knotfac * (vtstorm-vdrift * np.cos(pifac*latstorm) - v850storm)

    if bound == 'auto':
        bxmin = np.floor(longstorm[1]-deltax)
        bxmax = np.ceil(longstorm[1]+deltax)
        bymin = np.floor(latstorm[1]-deltay)
        bymax = np.ceil(latstorm[1]+deltay)

    h, hx, hy, x ,y = EstimateTopographicHeight(bxmin, bxmax, bymin, bymax, dellatlong)
    w = pointwfield(latstorm, longstorm, vstorm, rmstorm, vsestorm, rmsestorm, 
                    utstorm, vtstorm, ush, vsh, y, x, h, hx, hy)
    temp = eprecip * m_to_mm * 3600 * rowa_over_rowl * q900 * np.maximum(w[0,1,:,:]-wrad, 0)
    rainrate = temp
    rainrate = rainrate.transpose()
    return rainrate        
                

def windprofiles(vm, rm, r, wp):
    """
    This function return radial profiles of azimuthal wind, V, given matrices
    containing the maximum circular wind speed, vm (knots), the radii of maximum
    wind, rm (km), distances r (km) of each event from the points of interest, 
    and the wind profile wp (see below must be equal to 1, 2 or 3). 
    """    
    wprofile = wp                   # Use holland (1) or emanuel (2) or er2011 (3) wind profile
    wprofile = np.min([wprofile,3])   # Forbid use of experimental profiles    
    vm = vm*1852/3600             # Convert maximum wind speed to m/s
    
    if wprofile == 1:
        # Holland 2010 wind profile model        
        bs = 1.8
        rn = 300/20     # Holland parameters
        xn = 0.8
        
        rh = r/rm
        x = 0.5+(rh-1)*(xn-0.5)/(rn-1)
        # x = np.maximum(x,0.5)
        x[x<0.5] = 0.5
        V = vm*(rh^-bs*np.exp(1-rh^-bs))**x  
    
    elif wprofile == 2:
        # Outer radius (km)        
        r0=1000
        
        #  Re-scale radius of maximum winds by random number drawn from log-normal distribution
        # Shape parameters
        b = 0.25
        nb = 0.9
        mb = 1.6          
        mb2 = 2*mb
        fac1 = (1-b)*(mb+nb)
        fac2 = b*(1+2*mb)
        fac3 = 2*(mb+nb)
        fac4 = 2*mb+1  
        rat = r/np.maximum(rm,1)
        V = vm*(np.maximum((r0-r),0)/(r0-rm))*np.sqrt(rat**mb2*(fac1/(nb+mb*rat**fac3)+fac2/(1+mb2*rat**fac4)))
        
    elif wprofile == 3:
        crat = 1
        f = 5.0e-5
        f = f*1000  # effectively converts radii from kms to meters

        Mm = rm*vm+0.5*f*rm**2
        rn = r/rm
        if crat == 1:
            M = Mm*(2*rn**2/(1+rn**2))
        else:
            M = Mm*(2*rn**2/(2-crat+crat*rn**2))**(1/(2-crat))
        V = M / (r+1e-8)  # (Add long tail to V to avoid discontinuity in vorticity at outer radius) 3/2013
        V[V<0] = 0
        
    V = V * 3600 / 1852    # Convert wind speed to knots
    return V


def windprofilem(vm, rm, vm2, rm2, r, wp):
    """
    This function return radial profiles of azimuthal wind, V, given matrices
    containing the maximum circular wind speeds, vm and vm2 (knots), the radii of maximum
    wind, rm  and rm2 (km), associated with primary and secondary wind maxima
    distances r (km) of each event from the points of interest, 
    and the wind profile wp (see below must be equal to 1, 2, 3, or 4).
    """    
    wprofile = wp           # Use holland (1) or emanuel (2) or er2011 (3) or ec2013 (4) wind profile
    vm = vm * 1852/3600     # Convert maximum wind speed to m/s
    vm2 = vm2 * 1852/3600   # Convert maximum wind speed to m/s
    se = np.sum(rm2)        # Test if there are any secondary eyewalls
    # se = 0
    
    if wprofile == 1:
    # Holland 2010 wind profile model
        bs = 1.8
        rn = 300/20     # Holland parameters
        xn = 0.8
        #
        rh = r/rm
        x = 0.5+(rh-1)*(xn-0.5)/(rn-1)
        x[x<0.5] = 0.5
        V = vm * (rh**-bs * np.exp(1 - rh**-bs))**x   
        #
        if se != 0:
            rh = r/rm2
            x = 0.5+(rh-1)*(xn-0.5)/(rn-1)
            x[x<0.5] = 0.5
            V2 = vm2 * (rh**-bs * np.exp(1-rh**-bs))**x   
    
    elif wprofile == 2:
        # Outer radius (km)
        r0 = 1000
        # Shape parameters
        b = 0.25
        nb = 0.9
        mb = 1.6  
        mb2 = 2*mb
        fac1 = (1-b)*(mb+nb)
        fac2 = b*(1+2*mb)
        fac3 = 2*(mb+nb)
        fac4 = 2*mb+1 
        rat = r/np.maximum(rm,1)
        V = vm * (np.maximum((r0-r),0)/(r0-rm))*np.sqrt(rat**mb2 * (fac1/(nb+mb*rat**fac3)+fac2/(1+mb2*rat**fac4)))

        if se != 0:
            rat = r/np.maximum(rm2,1)
            V2 = vm2 * (np.maximum((r0-r),0)/(r0-rm2))*np.sqrt(rat**mb2 * (fac1/(nb+mb*rat**fac3)+fac2/(1+mb2*rat**fac4)))
    
    elif wprofile == 3:
        crat = 1
        f = 5.0e-5
        f = f*1000      # effectively converts radii from kms to meters
        Mm = rm * vm
        rn = r / rm
        if crat == 1:
            M = Mm * (2 * rn**2 / (1+rn**2))
        else:
            M = Mm * (2 * rn**2 / (2-crat+crat*rn**2))**(1/(2.0-crat))

        V = M / (r+1e-8)  # (Add long tail to V to avoid discontinuity in vorticity at outer radius) 3/2013
        V = np.maximum(V,0)
        if se != 0:
            Mm = rm2*vm2
            rn = r/rm2
            if crat == 1:
                M = Mm*(2*rn**2/(1+rn**2))
            else:
                M = Mm*(2*rn**2/(2-crat+crat*rn**2))**(1/(2.0-crat))
            V2 = M/(r+1e-8)  # (Add long tail to V to avoid discontinuity in vorticity at outer radius) 3/2013
            V2 = np.maximum(V2,0)
            
    if se != 0 & wprofile < 4:
        u = np.maximum(r,1)/rm
        hu = np.maximum(np.sign(u-1),0)
        V = V*(1-hu) + hu*vm/u
        u = r/np.maximum(rm2,1)
        hu = np.maximum(np.sign(u-1),0)
        V2 = V2*hu + (1-hu)*vm2*u
        V = np.maximum(V,V2)

    V = V*3600/1852    # Convert wind speed to knots
    
    return V
            
            
def utrans(latstore,longstore):
    """
    This function reads in track file, calculate translation speeds (knots) and smooth
    Args:
        latstore[nxm]: Latitude of points along each track
        longstore[nxm]: Longitude of points along each track
    """
    smfac  = 0.4
    pifac  = math.acos(-1)/180
    dfac   = 60 * 1.852
    dtfac  = 1/(4*1.852)
    netfac = dtfac * dfac
    nn, m  = np.shape(latstore)
    ut     = np.nan*np.zeros((nn,m))
    vt     = np.nan*np.zeros((nn,m))
    jmax   = np.zeros((nn))
    
    for n in range(nn):
        lat = latstore[n,:]
        # jm = np.argmin(np.abs(lat))     # For Emanuel output (matlab), get the length of each event by finding the first 0 element
        lat = lat[(lat != 0) & ~np.isnan(lat)]            
        jm = len(lat)        
        jm = np.maximum(jm,1)
        jmax[n] = jm
        long = longstore[n,:]
        longdif = long[2:jm] - long[0:jm-2]    # spatial difference in longitude
        for j in range(jm-2):
            if longdif[j] < -300:
                longdif[j] = longdif[j]+360
            elif longdif[j] > 300:
                longdif[j] = longdif[j]-360

        ut[n,1:jm-1] = netfac * np.cos(pifac*lat[1:jm-1])*longdif
        vt[n,1:jm-1] = netfac * (latstore[n,2:jm]-latstore[n,0:jm-2])
        j2 = np.minimum(2,jm)
        j3 = np.maximum(jm-3,0)
        ut[n,0] = 2*ut[n,1] - ut[n,j2]
        vt[n,0] = 2*vt[n,1] - vt[n,j2]        
        if jm > 0:
            ut[n,jm-1] = 2.*ut[n,jm-2]-ut[n,j3]
            vt[n,jm-1] = 2.*vt[n,jm-2]-vt[n,j3]

        #  Smooth translation velocity
        vt[n,1:jm-1] = vt[n,1:jm-1]+smfac*(vt[n,0:jm-2]+vt[n,2:jm]-2*vt[n,1:jm-1])
        ut[n,1:jm-1] = ut[n,1:jm-1]+smfac*(ut[n,0:jm-2]+ut[n,2:jm]-2*ut[n,1:jm-1])

    return ut, vt, jmax


def pointwfield(latstore,longstore,vstore,rmstore,vsestore,
                rmsestore,ut,vt,us,vs,plat,plong,h,hx,hy):
    """
    This function takes the nn X m 'store' matrices  and, using the 
    latitude (plat) and longitude (plong) vectors of points of interest, 
    calculates the spatial distribution of vertical velocity centered at 
    plong, plat. The routine utrans should be run before this one. 
    Args:
        us, vs: vertical shears (m/s) used to estimate the baroclinic 
                components of the vertical motion
        h, hx, hy: contain topographic heights and their derivatives in x and y
        latstore[nxm], longstore[nxm]: Latitude and longitude of points along each track
        vstore [nxm]: Maximum circular wind speed along each track.
        vsestore [nxm]: The maximum circular wind speed of maximum circular wind of any 
                        secondary eyewalls that may be present, along each track.
        rmstore [nxm]: The radius (km) of maximum circular wind along each track
        rmsestore [nxm]: The radius (km) of maximum circular wind of any secondary wind 
                         maxima present (0 if absent), along each track.
    
    Returns:
        w [nxm]: spatial distribution of vertical velocity centered at plong, plat                  
    """
    # Load parameters
    deltar = params.deltar
    timeresw = params.timeresw
    Htrop = params.Htrop
    radcity = params.radcity
    wprofile = params.wprofile
    
    se = np.max(vsestore)  # Check for secondary eyewalls
    deltari = 1. / deltar
    timereswi = 1. / (3600 * timeresw)
    Hi = 1. / Htrop
    nn = 1
    m = len(ut)
    sx = len(plong)
    sy = len(plat)
    
    ngrid=0
    if sx == sy:
        ngrid = 1
    
    if ngrid == 1:
        sy = 1
    else:
        sy = np.max(np.shape(plat))
    
    w = np.zeros((nn,3,sx,sy))
    w2 = np.zeros((nn,3,sx,sy))
    
    if plong[0] < 0:
        plong=plong + 360
    
    pifac = math.acos(-1)/180         # pi number
    dfac = 60 * 1.852
    sfac = 1 / (0.25 * 60 * 1852)
    knotfac = 1852. / 3600            # 1 knots = 0.5144 m/s
    latfac = latstore[0] / (abs(latstore[0])+1e-8)
    ut = ut * knotfac                 # convert knots to m/s
    vt = vt * knotfac                 # convert knots to m/s
    Omega = math.acos(-1)/(6*3600)    # Earth angular velocity parameter
    
    # Estimate Drag coeffs:
    cdrag, cdx, cdy = EstimateDragCoefficients(plat, plong, sfac)
    cdrag = np.maximum(cdrag, 0.0015)    # To put lower bound on drag coefficeint over water
    rfull, dx, dy = CalculateDistancePOTfromTrack(plat, plong, latstore, longstore, nn, m, sx, sy, ngrid, dfac)
    rfull = np.maximum(rfull, 0.5)

    denom = 1/rfull
    costhetafull = dx*denom
    sinthetafull = dy*denom

    vfull = np.zeros((nn, 3, sx, sy))
    rmfull = np.zeros((nn, 3, sx, sy))
    vsefull = np.zeros((nn, 3, sx, sy))
    rmsefull = np.zeros((nn, 3, sx, sy))
    latfull = np.zeros((nn, 3, sx, sy))
    longfull = np.zeros((nn, 3, sx, sy))
    utfull = np.zeros((nn, 3, sx, sy))
    vtfull = np.zeros((nn, 3, sx, sy))
    usfull = np.zeros((nn, 3, sx, sy))
    vsfull = np.zeros((nn, 3, sx, sy))
    hfull = np.zeros((nn, 3, sx, sy))
    hyfull = np.zeros((nn, 3, sx, sy))
    hxfull = np.zeros((nn, 3, sx, sy))
    cdfull = np.zeros((nn, 3, sx, sy))
    cdyfull = np.zeros((nn, 3, sx, sy))
    cdxfull = np.zeros((nn, 3, sx, sy))
    
    for i in range(sx):
        for jj in range(sy):
            if ngrid == 1:
                j = i
            else:
                j = jj                        
            vfull[0,:,i,j] = vstore[:]
            rmfull[0,:,i,j] = rmstore[:]
            vsefull[0,:,i,j] = vsestore[:]
            rmsefull[0,:,i,j] = rmsestore[:]
            latfull[0,:,i,j] = latstore[:]
            longfull[0,:,i,j] = longstore[:]
            utfull[0,:,i,j] = ut[:]
            vtfull[0,:,i,j] = vt[:]
            usfull[0,:,i,j] = us[:]
            vsfull[0,:,i,j] = vs[:]

    w, cp, V, Vd, Vrp, Vrm, u1temp, u2temp, Cdp, Cdm = \
        CalculateWindPrimary(utfull, vtfull, vfull, rfull, rmfull, vsefull, rmsefull, 
                             latfull, longfull, w, nn, 3, sx, sy, h, hx, hy, cdrag, hfull,
                             hxfull, hyfull, cdfull, cdx, cdxfull, cdy, cdyfull, 
                             Hi, Htrop, Omega, pifac, deltar, deltari, knotfac, latfac, timereswi, 
                             costhetafull, sinthetafull, 2, wprofile, adj_water=True)   
    
    if se > 0:   #  If secondary eyewalls present, add in their contribution to w    
        w2, cp, V, Vd, Vrp, Vrm = \
            CalculateWindSecondary(utfull, vtfull, vfull, rfull, rmfull, vsefull, rmsefull, 
                             latfull, longfull, w, nn, 3, sx, sy, Hi, Htrop, Omega, 
                             pifac, deltar, deltari, knotfac, latfac, timereswi, 
                             costhetafull, sinthetafull,  u1temp, u2temp, Cdp, Cdm, wprofile)        
        w = np.maximum(w,w2)

    hxmod = -0.0005 * (cp + 2*V/(0.1+rfull) + deltari*(Vrp-Vrm)) * vsfull # to include shear dotted with storm entropy
    hymod =  0.0005 * (cp + 2*V/(0.1+rfull) + deltari*(Vrp-Vrm)) * usfull
    ufunc = (radcity-rfull)/50
    ufunc = np.minimum(ufunc,1)
    ufunc = np.maximum(ufunc,0)
    utfull = utfull * ufunc
    vtfull = vtfull * ufunc
    w[:,1,:,:] = w[:,1,:,:] \
        + (vtfull[:,1,:,:] + Vd[:,1,:,:]*latfac*costhetafull[:,1,:,:]) * hyfull[:,1,:,:] \
        + (utfull[:,1,:,:] - Vd[:,1,:,:]*latfac*sinthetafull[:,1,:,:]) * hxfull[:,1,:,:] \
        + Vd[:,1,:,:]*costhetafull[:,1,:,:]*hymod[:,1,:,:] \
        - Vd[:,1,:,:]*sinthetafull[:,1,:,:]*hxmod[:,1,:,:]
    w = np.minimum(w,7)
    
    return w


def vouternew(vm, fc, ro, wc, CD, q):
    """
    Numerical integration of the outer wind profile from simple ODE
    Args:
        vm: maximum wind (m/s) 
        f: Coriolis force (s^-1) 
        ro: outer radiusin (km) 
        wc: radiative subsidence rate (mm/s) 
        CD: drag coefficient 
        q: number of radial points

    Returns:  
        v: potential intensity (m/s) 
        r: radius at maximum wind (km)
        imin: minimum index for which V and r are defined

    """    
    # Asselin filter coefficient
    assl = 0.2

    # Convert units
    ro = ro * 1000.0  # Convert to meters
    wc = wc * 0.001   # Convert to m/s
    chi = CD * fc * ro / wc  # Definition of chi
    rdim = vm / fc
    rond = ro / rdim
    dr = rond / (q - 2)

    # Initialize arrays
    m = [0.0] * q
    v = [0.0] * q
    r = [0.0] * q

    # Integrate outer wind ODE
    rnd = rond - dr
    m[q - 2] = 0.25 * (rond**2 - rnd**2)
    v[q - 2] = m[q - 2] / (rond - dr)
    r[q - 2] = rond - dr

    for i in range(q - 3, 0, -1):
        r[i] = r[i+1] - dr
        m[i] = m[i+2] - 2.0 * dr * (chi * m[i+1]**2 / (rond**2 - r[i+1]**2) - r[i+1])
        m[i+1] = m[i+1] + assl * ( m[i] + m[i+2] - 2.0 * m[i+1])
        v[i] = m[i] / r[i]
        if v[i] > 1.0:
            imin = i
            break

    # Fill in values inside radius where v = potential intensity
    for i in range(imin-1):
        v[i] = 0.0
        r[i] = dr * i

    # Re-dimensionalize
    for i in range(q):
        v[i] = vm * v[i]                # v in m/s
        r[i] = 0.001 * vm * r[i] / fc   # r in km

    epa = (v[imin]-vm)/(v[imin]-v[imin+1])
    rm = r[imin]+epa*(r[imin+1]-r[imin])

    return v, rm, imin



def rainswathx(nt,latstore,longstore,rmstore,vstore,rmsestore,vsestore,ut,vt,u850store,v850store):

    magfac = params.magfac
    deltax = params.deltax
    deltay = params.deltay
    bxmin = params.bxmin
    bxmax = params.bxmax
    bymin = params.bymin
    bymax = params.bymax
    dellatlongs = params.dellatlongs
    q900 = params.q900
    eprecip = params.eprecip
    wrad = params.wrad
    timeres = params.timeres

    m_to_mm = 1000
    rowa_over_rowl = 0.00117
        
    nrm, mrm = np.shape(rmstore)
    rfac = magfac * (1+np.zeros((nrm,mrm)))
        
    if bxmin < 0:
        bxmin += 360
    if bxmax < 0:
        bxmax += 360
    
    # Initialize variables
    latdata = latstore[nt,:]
    latdata = latdata[(latdata != 0) & ~np.isnan(latdata)]
    latsize = len(latdata)
    utd = ut[nt,:latsize].reshape((1,latsize))
    vtd = vt[nt,:latsize].reshape((1,latsize))
    ush = np.zeros_like(utd)
    vsh = np.zeros_like(vtd)    
    vdrift = 1.5 * 3600 / 1852
    vdrift *= latstore[0, 0] / (np.abs(latstore[0, 0]) + 1e-8)

    if 'u850store' in locals():
        ush = 5 * 1852 / 3600 * (utd - u850store[nt, :latsize])
        vsh = 5 * 1852 / 3600 * (vtd - vdrift * np.cos(np.pi / 180 * latstore[nt, :latsize]) - v850store[nt, :latsize])

    lat = latstore[nt,:latsize].reshape((1,latsize))
    long = longstore[nt,:latsize].reshape((1,latsize))
    long[long<0] += 360         # Convert longitude to 0 to 360 degree east
    v = vstore[nt,:latsize].reshape((1,latsize))
    vse = vsestore[nt,:latsize].reshape((1,latsize))
        
    # Scale and randomize radii of maximum wind
    nrm, mrm = rmstore.shape
    jmaxd = latsize
    rfac = np.ones((nrm, mrm))
    
    temp = magfac * rfac[nt,:] * rmstore[nt,:]
    rm = temp[:latsize].reshape((1,latsize))
    nrm = np.shape(rm)[1]
    temp = 0 * magfac * rfac[nt,:] * rmsestore[nt,:]
    rmse = temp[:latsize].reshape((1,latsize))
        
    for i in range(jmaxd):
        if long[0,0] > 200 and long[0,i] < 50:
            long[0,i] += 360
    
    bxmin = np.min(long[np.nonzero(long)]) - deltax
    bxmax = np.max(long[np.nonzero(long)]) + deltax
    bymin = np.min(lat[np.nonzero(lat)]) - deltay
    bymax = np.max(lat[np.nonzero(lat)]) + deltay
        
    h,hx,hy,x,y = EstimateTopographicHeight(bxmin,bxmax,bymin,bymax,dellatlongs)
    w = pointwshortn(lat, long, v, rm, vse, rmse, utd, vtd, ush, vsh, y, x, h, hx, hy, timeres)
    wq = np.maximum(w-wrad,0) * q900
    netrain = eprecip * m_to_mm * timeres * 3600 * rowa_over_rowl * np.sum(wq,axis=(0,1))
    
    return x,y,netrain.T



def pointwshortn(latstore,longstore,vstore,rmstore,vsestore,rmsestore,ut,vt,us,vs,plat,plong,h,hx,hy,timeres):
    """
    This function take the nn X m 'store' matrices  and,
    using the latitude (plat) and longitude (plong) vectors of points of 
    interest, calculates time series of vertical velocity at the points of
    interest. The arrays 'us' and 'vs' are the vertical shears (m/s) used 
    to estimate the baroclinic components of the vertical motion. 
    """    
    
    deltar = params.deltar
    timelength = params.timelength
    timeresi = 1. / (3600 * timeres)
    Htrop = params.Htrop
    radcity = params.radcity
    wprofile = params.wprofile
        
    pifac = np.arccos(-1) / 180
    dfac = 60.0 * 1.852
    sfac = 1.0 / (0.25 * 60.0 * 1852)
    knotfac = 1852.0 / 3600.0
    Omega = math.acos(-1)/(6*3600)    # Earth angular velocity parameter
                
    se = 0  # Test for secondary eyewalls
    ntime = int(timelength / timeres + 1)
    nsteps = round(2.0 / timeres)
    nstepsi = 1.0 / nsteps
    delj = np.floor(timelength/4).astype(int)
    deltari = 1.0 / deltar
    Hi = 1.0 / Htrop
    
    nn, m = ut.shape
    sx = plong.size
    sy = plat.size
    ngrid = 0
    
    if sx == sy:
        ngrid = 1
        sy = 1    
        
    w = np.zeros((nn, ntime, sx, sy))
    w2 = np.zeros((nn, ntime, sx, sy))

    if plong[0] < 0:
        plong = plong + 360
        
    latfac = latstore[0,0] / (abs(latstore[0,0])+1e-8)
    ut = ut * knotfac
    vt = vt * knotfac

    cdrag, cdx, cdy = EstimateDragCoefficients(plat, plong, sfac)
    # cdrag = np.maximum(cdrag, 0.0015)    # To put lower bound on drag coefficeint over water
    
    # Three statements below added Dec. 2018 to reduce drag coefficient of near-coastal locations
    cdrag = np.minimum(cdrag, 1.5e-3 * (1 + np.maximum(h, 0) / 100))
    cdx = cdx * 0.01 * np.minimum(np.maximum(h, 0), 100)
    cdy = cdy * 0.01 * np.minimum(np.maximum(h, 0), 100)    
        
    radius,dx,dy = CalculateDistancePOTfromTrack(plat,plong,latstore,longstore,nn,m,sx,sy,ngrid,dfac)
    radius = np.maximum(radius,0.5)

    jmin = np.argmin(radius, axis=1)
    jmin = np.maximum(jmin, 0+delj)
    jmin = np.minimum(jmin, m-1-delj)
    jstart = jmin - delj
    jend = jmin + delj + 1
    jtot = 2 * delj + 1
    jfine = 1 + nsteps * (jtot - 1)

    # Create arrays for reduced length time series of various quantities
    vshort = np.zeros((nn,jtot,sx,sy))
    rmshort = np.zeros((nn,jtot,sx,sy))
    vseshort = np.zeros((nn,jtot,sx,sy))
    rmseshort = np.zeros((nn,jtot,sx,sy))
    rshort = np.zeros((nn,jtot,sx,sy))
    latshort = np.zeros((nn,jtot,sx,sy))
    longshort = np.zeros((nn,jtot,sx,sy))
    utshort = np.zeros((nn,jtot,sx,sy))
    vtshort = np.zeros((nn,jtot,sx,sy))
    usshort = np.zeros((nn,jtot,sx,sy))
    vsshort = np.zeros((nn,jtot,sx,sy))

    for i in range(sx):
        for j in range(sy):
            for n in range(nn):
                vshort[n,:,i,j] = vstore[n,jstart[n,i,j]:jend[n,i,j]]
                rmshort[n,:,i,j] = rmstore[n,jstart[n,i,j]:jend[n,i,j]]
                vseshort[n,:,i,j] = vsestore[n,jstart[n,i,j]:jend[n,i,j]]
                rmseshort[n,:,i,j] = rmsestore[n,jstart[n,i,j]:jend[n,i,j]]
                rshort[n,:,i,j] = radius[n,jstart[n,i,j]:jend[n,i,j],i,j]
                latshort[n,:,i,j] = latstore[n,jstart[n,i,j]:jend[n,i,j]]
                longshort[n,:,i,j] = longstore[n,jstart[n,i,j]:jend[n,i,j]]
                utshort[n,:,i,j] = ut[n,jstart[n,i,j]:jend[n,i,j]]
                vtshort[n,:,i,j] = vt[n,jstart[n,i,j]:jend[n,i,j]]
                usshort[n,:,i,j] = us[n,jstart[n,i,j]:jend[n,i,j]]
                vsshort[n,:,i,j] = vs[n,jstart[n,i,j]:jend[n,i,j]]
                
    # Create arrays filled with zeros
    vfine = np.zeros((nn,jfine,sx,sy))
    rmfine = np.zeros((nn,jfine,sx,sy))
    vsefine = np.zeros((nn,jfine,sx,sy))
    rmsefine = np.zeros((nn,jfine,sx,sy))
    latfine = np.zeros((nn,jfine,sx,sy))
    longfine = np.zeros((nn,jfine,sx,sy))
    utfine = np.zeros((nn,jfine,sx,sy))
    vtfine = np.zeros((nn,jfine,sx,sy))
    usfine = np.zeros((nn,jfine,sx,sy))
    vsfine = np.zeros((nn,jfine,sx,sy))
    hfine = np.zeros((nn,jfine,sx,sy))
    hyfine = np.zeros((nn,jfine,sx,sy))
    hxfine = np.zeros((nn,jfine,sx,sy))
    cdfine = np.zeros((nn,jfine,sx,sy))
    cdyfine = np.zeros((nn,jfine,sx,sy))
    cdxfine = np.zeros((nn,jfine,sx,sy))

    k = 0
    for j in range(jtot - 1):
        for n in range(nsteps):
            weight = (n * nstepsi)
            vfine[:, k, :, :] = (1 - weight) * vshort[:, j, :, :] + weight * vshort[:, j + 1, :, :]
            rmfine[:, k, :, :] = (1 - weight) * rmshort[:, j, :, :] + weight * rmshort[:, j + 1, :, :]
            vsefine[:, k, :, :] = (1 - weight) * vseshort[:, j, :, :] + weight * vseshort[:, j + 1, :, :]
            rmsefine[:, k, :, :] = (1 - weight) * rmseshort[:, j, :, :] + weight * rmseshort[:, j + 1, :, :]
            latfine[:, k, :, :] = (1 - weight) * latshort[:, j, :, :] + weight * latshort[:, j + 1, :, :]
            longfine[:, k, :, :] = (1 - weight) * longshort[:, j, :, :] + weight * longshort[:, j + 1, :, :]
            utfine[:, k, :, :] = (1 - weight) * utshort[:, j, :, :] + weight * utshort[:, j + 1, :, :]
            vtfine[:, k, :, :] = (1 - weight) * vtshort[:, j, :, :] + weight * vtshort[:, j + 1, :, :]
            usfine[:, k, :, :] = (1 - weight) * usshort[:, j, :, :] + weight * usshort[:, j + 1, :, :]
            vsfine[:, k, :, :] = (1 - weight) * vsshort[:, j, :, :] + weight * vsshort[:, j + 1, :, :]
            k += 1
            
    vfine[:, k, :, :] = vshort[:, jtot-1, :, :]
    rmfine[:, k, :, :] = rmshort[:, jtot-1, :, :]
    vsefine[:, k, :, :] = vseshort[:, jtot-1, :, :]
    rmsefine[:, k, :, :] = rmseshort[:, jtot-1, :, :]
    latfine[:, k, :, :] = latshort[:, jtot-1, :, :]
    longfine[:, k, :, :] = longshort[:, jtot-1, :, :]
    utfine[:, k, :, :] = utshort[:, jtot-1, :, :]
    vtfine[:, k, :, :] = vtshort[:, jtot-1, :, :]
    usfine[:, k, :, :] = usshort[:, jtot-1, :, :]
    vsfine[:, k, :, :] = vsshort[:, jtot-1, :, :]

    rmsefine = np.maximum(rmsefine,0.1)
    rfine,dx,dy = CalculateDistancePOTfromTrack(plat,plong,latfine,longfine,nn,jfine,sx,sy,ngrid,dfac)
    rfinei = 1 / np.maximum(rfine,1)

    w, cp, V, Vd, Vrp, Vrm, u1temp, u2temp, Cdp, Cdm = \
        CalculateWindPrimary(utfine, vtfine, vfine, rfine, rmfine, vsefine, rmsefine, 
                             latfine, longfine, w, nn, jfine, sx, sy, h, hx, hy, cdrag, hfine,
                             hxfine, hyfine, cdfine, cdx, cdxfine, cdy, cdyfine, 
                             Hi, Htrop, Omega, pifac, deltar, deltari, knotfac, latfac, timeresi, 
                             dx*rfinei, dy*rfinei, 10, wprofile, adj_water=False)
                        
    if se > 0:   #  If secondary eyewalls present, add in their contribution to w    
        w2, cp, V, Vd, Vrp, Vrm = \
            CalculateWindSecondary(utfine, vtfine, vfine, rfine, rmfine, vsefine, rmsefine, 
                             latfine, longfine, w, nn, jfine, sx, sy, Hi, Htrop, Omega, 
                             pifac, deltar, deltari, knotfac, latfac, timeresi, 
                             dx*rfinei, dy*rfinei,  u1temp, u2temp, Cdp, Cdm, wprofile)                       
        w = np.maximum(w, w2)       
        
    # Now add in topographic and shear components        
    hxmod = -0.0005*(cp+2*V/(0.1+rfine)+deltari*(Vrp-Vrm))*vsfine # to include shear dotted with storm entropy
    hymod =  0.0005*(cp+2*V/(0.1+rfine)+deltari*(Vrp-Vrm))*usfine
             
    # Reduce effect of translation speed outside of radcity
    ufunc = (radcity-rfine) / 50.
    ufunc = np.minimum(ufunc,1)
    ufunc = np.maximum(ufunc,0)
    utfine = utfine * ufunc
    vtfine = vtfine * ufunc
    
    # Reduce effect of orography outside of storm core
    ufunc = (150-rfine) / 30
    ufunc = np.minimum(ufunc,0.6)
    ufunc = np.maximum(ufunc,0.2)
    hxfine = hxfine * ufunc
    hyfine = hyfine * ufunc    
    
    w = w + (Vd*latfac*dx*rfinei + vtfine)*hyfine + (utfine - Vd*latfac*dy*rfinei)*hxfine + Vd*dx*rfinei*hymod - Vd*dy*rfinei*hxmod
    w = np.minimum(w,7)
    return w



def raingen(plat, plong, latstore, longstore, vstore, rmstore, vsestore, rmsestore, u850store, v850store, ut, vt):
    plat = np.array([plat])
    plong = np.array([plong])
    magfac = params.magfac
    sx = plong.size    
    sy = plat.size
    
    mat = sio.loadmat('data/bathymetry_high.mat')
    bathy = mat['bathy']
    ntopo,_ = np.shape(bathy)
    topores = 360/ntopo
    toporesi = 1/topores
    sfac = 1./(topores*60.*1852)
    pifac = math.acos(-1)/180
    knotfac = 1852./3600    
    m, n = ut.shape
    
    ush = np.zeros((n,m))
    vsh = np.zeros((n,m))
    vdrift = 1.5*3600/1852
    vdrift = vdrift*latstore[0,0]/(abs(latstore[0,0])+1e-8)
    if 'u850store' in locals():
        ush = 5 * knotfac * (ut - u850store)
        vsh = 5 * knotfac * (vt - vdrift * np.cos(pifac*latstore) - v850store)    
    
    lat = latstore.copy()
    long = longstore.copy()
    v = vstore.copy()
    vse = vsestore.copy()
    
    nrm, mrm = rmstore.shape
    rfac = magfac * np.ones((nrm, mrm))
    
    rm = rmstore * rfac 
    rmse = rmsestore * rfac 
    
    m_to_mm = 1000
    rowa_over_rowl = 0.00117
    bathy = np.maximum(bathy,-1)
    
    h, hx, hy = CalculateSpatialDerivatives(bathy,plong,plat,sx,1,sfac,pifac,ntopo,toporesi)

    return h, hx, hy
