"""
Functions for wind in PyTCR
"""


import os
import sys
import math
import numpy as np
from tcr import terrain_boundary as tcr_tb
from tcr import iodata as tcr_io
from tcr import params


def calculate_wind_primary(utf, vtf, vf, rf, rmf, vsef, rmsef, latf,
                           nn, jf, sx, sy, h, hx, hy, cdrag, hf, hxf, hyf, cdf,
                           cdx, cdxf, cdy, cdyf, Hi, Htrop, omega, pifac,
                           deltar, deltari, knotfac, latfac, timereswi,
                           costhetaf, sinthetaf, thresM, wprofile,
                           adj_water=False):
    """
    This function calculate the wind in the primary band.

    Inputs:
    -------
        - utf, vtf: west-east and storm translation velocity
        - vf    : maximum circular wind speed at each 2-hour point along each track
        - rf    : Euclidean distance from each track
        - rmf   : radius (km) of maximum circular wind of 2-hour points along each track
        - vsef  : maximum circular wind speed  of any secondary eyewalls that may be present at
                  each 2-hour point along each track
        - rmsef : radius (km) of maximum circular wind of any secondary wind maxima of 2-hour points
                  along each track
    Returns:
    --------
    """

    V = np.zeros((nn, jf, sx, sy))
    Vd = np.zeros((nn, jf, sx, sy))
    Vpp = np.zeros((nn, jf, sx, sy))
    Vmp = np.zeros((nn, jf, sx, sy))
    Vpm = np.zeros((nn, jf, sx, sy))
    Vmm = np.zeros((nn, jf, sx, sy))
    Vrp = np.zeros((nn, jf, sx, sy))
    Vrm = np.zeros((nn, jf, sx, sy))
    up = np.zeros((nn, jf, sx, sy))
    um = np.zeros((nn, jf, sx, sy))
    w = np.zeros((nn, jf, sx, sy))

    V = windprofiles(vf, rmf, rf, wprofile)
    Vd = windprofiles(vf, rmf, rf, wprofile, vsef, rmsef)
    Vpp[:, 1:jf-1, :, :] = windprofiles(vf[:, 2:jf, :, :],
                                        rmf[:, 2:jf, :, :],
                                        rf[:, 1:jf-1, :, :] + deltar,
                                        wprofile)
    Vpm[:, 1:jf-1, :, :] = windprofiles(vf[:, 2:jf, :, :],
                                        rmf[:, 2:jf, :, :],
                                        np.maximum(rf[:, 1:jf-1, :, :]-deltar, 0),
                                        wprofile)
    Vmp[:, 1:jf-1, :, :] = windprofiles(vf[:, 0:jf-2, :, :],
                                        rmf[:, 0:jf-2, :, :],
                                        rf[:, 1:jf-1, :, :] + deltar,
                                        wprofile)
    Vmm[:, 1:jf-1, :, :] = windprofiles(vf[:, 0:jf-2, :, :],
                                        rmf[:, 0:jf-2, :, :],
                                        np.maximum(rf[:, 1:jf-1, :, :]-deltar, 0),
                                        wprofile)
    Vrp[:, :, :, :] = windprofiles(vf[:, :, :, :],
                                   rmf[:, :, :, :],
                                   rf[:, :, :, :]+deltar,
                                   wprofile)
    Vrm[:, :, :, :] = windprofiles(vf[:, :, :, :],
                                   rmf[:, :, :, :],
                                   np.maximum(rf[:, :, :, :] - deltar, 0),
                                   wprofile)

    # Convert to meters per second. Done now because windprofile expects knots.
    V = knotfac * V
    Vd = knotfac * Vd
    Vpp = knotfac * Vpp
    Vpm = knotfac * Vpm
    Vmp = knotfac * Vmp
    Vmm = knotfac * Vmm
    Vrp = knotfac * Vrp
    Vrm = knotfac * Vrm

    vph = 0.5 * (V + Vrp)
    vmh = 0.5 * (V + Vrm)
    u1temp = vtf * costhetaf - utf * sinthetaf
    u2temp = vtf**2 + utf**2
    vnetp = np.sqrt(vph**2 + 2 * vph * latfac * u1temp + u2temp)
    vnetm = np.sqrt(vmh**2 + 2 * vmh * latfac * u1temp + u2temp)

    for n in range(nn):
        for j in range(jf):
            hf[n, j, :, :] = h[:, :]
            hyf[n, j, :, :] = hy[:, :]
            hxf[n, j, :, :] = hx[:, :]
            cdf[n, j, :, :] = cdrag[:, :]
            cdyf[n, j, :, :] = cdy[:, :]
            cdxf[n, j, :, :] = cdx[:, :]

    cdfac = 1e3 * 0.5 * deltar * (cdxf*costhetaf + cdyf*sinthetaf)
    Cdp = cdf+cdfac
    Cdm = cdf-cdfac
    Cdp = np.maximum(Cdp, 0)
    Cdp = np.minimum(Cdp, 0.005)
    Cdm = np.maximum(Cdm, 0)
    Cdm = np.minimum(Cdm, 0.005)

    #  These lines of code added/modified  December 2018 to account for roughness changes over water
    if adj_water:
        facp = 1 + 0.0193 * vnetp
        facm = 1 + 0.0193 * vnetm
        facp[hf > 0] = 1  # Do not let roughness change wind over land
        facm[hf > 0] = 1
        uekp = -Hi * Cdp * facp * vph * vnetp
        uekm = -Hi * Cdm * facm * vmh * vnetm
    else:
        uekp = -Hi * Cdp * vph * vnetp
        uekm = -Hi * Cdm * vmh * vnetm

    cp = 1000 * 2 * omega * np.sin(pifac*abs(latf))
    dMdrp = cp * (rf+0.5*deltar) + (rf+0.5*deltar) * deltari * (Vrp-V) + 0.5 * (Vrp+V)
    dMdrp = np.maximum(dMdrp, thresM)
    dMdrm = cp * (rf-0.5*deltar) + (rf-0.5*deltar) * deltari * (V-Vrm) + 0.5 * (Vrm+V)
    dMdrm = np.maximum(dMdrm, thresM)

    rmf_safe = np.where(rmf[:, 1:jf-1, :, :] == 0, 1e-32, rmf[:, 1:jf-1, :, :])
    efacp = np.minimum((-1 + 2*((rf[:, 1:jf-1, :, :]+deltar) / rmf_safe)**2), 1)
    efacm = np.minimum((-1 + 2*((rf[:, 1:jf-1, :, :]-deltar) / rmf_safe)**2), 1)

    up[:, 1:jf - 1, :, :] = (
        (rf[:, 1:jf - 1, :, :] + deltar)
        * (
            -0.5 * timereswi * efacp
            * (Vpp[:, 1:jf - 1, :, :] - Vmp[:, 1:jf - 1, :, :])
            + uekp[:, 1:jf - 1, :, :]
        )
        / dMdrp[:, 1:jf - 1, :, :]
    )
    um[:, 1:jf - 1, :, :] = (
        (rf[:, 1:jf - 1, :, :] - deltar)
        * (
            -0.5 * timereswi * efacm
            * (Vpm[:, 1:jf - 1, :, :] - Vmm[:, 1:jf - 1, :, :])
            + uekm[:, 1:jf - 1, :, :]
        )
        / dMdrm[:, 1:jf - 1, :, :]
    )
    w[:, 1:jf - 1, :, :] = (
        -Htrop * deltari
        * (
            (rf[:, 1:jf - 1, :, :] + deltar) * up[:, 1:jf - 1, :, :]
            - (rf[:, 1:jf - 1, :, :] - deltar) * um[:, 1:jf - 1, :, :]
        )
        / np.maximum(rf[:, 1:jf - 1, :, :], 1)
    )
    return w, cp, V, Vd, Vrp, Vrm, u1temp, u2temp, Cdp, Cdm


def calculate_wind_secondary(rf, vsef, rmsef, latf, nn, jf, sx, sy, Hi, Htrop, omega, pifac,
                             deltar, deltari, knotfac, latfac, timereswi, u1temp, u2temp, Cdp,
                             Cdm, wprofile):
    """
    This function calculate wind for secondary eyewalls

    Inputs:
        - utf:

    Returns:
        - w, cp, V, Vd, Vrp, Vrm:
    """

    V = np.zeros((nn, jf, sx, sy))
    Vd = np.zeros((nn, jf, sx, sy))
    Vpp = np.zeros((nn, jf, sx, sy))
    Vmp = np.zeros((nn, jf, sx, sy))
    Vpm = np.zeros((nn, jf, sx, sy))
    Vmm = np.zeros((nn, jf, sx, sy))
    Vrp = np.zeros((nn, jf, sx, sy))
    Vrm = np.zeros((nn, jf, sx, sy))
    up = np.zeros((nn, jf, sx, sy))
    um = np.zeros((nn, jf, sx, sy))
    w2 = np.zeros((nn, jf, sx, sy))

    V = windprofiles(vsef, rmsef, rf, wprofile)
    Vpp[:, 1:jf-1, :, :] = windprofiles(vsef[:, 2:jf, :, :],
                                        rmsef[:, 2:jf, :, :],
                                        rf[:, 1:jf-1, :, :] + deltar,
                                        wprofile)
    Vpm[:, 1:jf-1, :, :] = windprofiles(vsef[:, 2:jf, :, :], rmsef[:, 2:jf, :, :],
                                        np.maximum(rf[:, 1:jf-1, :, :]-deltar, 0),
                                        wprofile)
    Vmp[:, 1:jf-1, :, :] = windprofiles(vsef[:, 0:jf-2, :, :],
                                        rmsef[:, 0:jf-2, :, :],
                                        rf[:, 1:jf-1, :, :]+deltar,
                                        wprofile)
    Vmm[:, 1:jf-1, :, :] = windprofiles(vsef[:, 0:jf-2, :, :],
                                        rmsef[:, 0:jf-2, :, :],
                                        np.maximum(rf[:, 1:jf-1, :, :]-deltar, 0),
                                        wprofile)
    Vrp = windprofiles(vsef, rmsef, rf+deltar, wprofile)
    Vrm = windprofiles(vsef, rmsef, np.maximum(rf-deltar, 0), wprofile)

    # Convert to meters per second because windprofile expects knots.
    V = knotfac * V
    Vd = knotfac * Vd
    Vpp = knotfac * Vpp
    Vpm = knotfac * Vpm
    Vmp = knotfac * Vmp
    Vmm = knotfac * Vmm
    Vrp = knotfac * Vrp
    Vrm = knotfac * Vrm

    vph = 0.5 * (V + Vrp)
    vmh = 0.5 * (V + Vrm)
    vnetp = np.sqrt(vph**2 + 2 * vph * latfac * u1temp + u2temp)
    vnetm = np.sqrt(vmh**2 + 2 * vmh * latfac * u1temp + u2temp)
    uekp = -Hi * Cdp * vph * vnetp
    uekm = -Hi * Cdm * vmh * vnetm

    cp = 1000 * 2 * omega * np.sin(pifac * abs(latf))
    dMdrp = (
        cp * (rf + 0.5 * deltar)
        + (rf + 0.5 * deltar) * deltari * (Vrp - V)
        + 0.5 * (Vrp + V)
    )
    dMdrp = np.maximum(dMdrp, 2)
    dMdrm = (
        cp * (rf - 0.5 * deltar)
        + (rf - 0.5 * deltar) * deltari * (V - Vrm)
        + 0.5 * (Vrm + V)
    )
    dMdrm = np.maximum(dMdrm, 2)
    efacp = np.minimum(
        (-1 + 2 * ((rf[:, 1:jf-1, :, :] + deltar) / rmsef[:, 1:jf-1, :, :]) ** 2), 1
    )
    efacm = np.minimum(
        (-1 + 2 * ((rf[:, 1:jf-1, :, :] - deltar) / rmsef[:, 1:jf-1, :, :]) ** 2), 1
    )

    #  Do not consider time rate of change when either velocity is zero
    efacp = efacp*np.minimum(Vpp[:, 1:jf-1, :, :], 1) * np.minimum(Vmp[:, 1:jf-1, :, :], 1)
    efacm = efacm*np.minimum(Vpm[:, 1:jf-1, :, :], 1) * np.minimum(Vmm[:, 1:jf-1, :, :], 1)
    up[:, 1:jf-1, :, :] = (rf[:, 1:jf-1, :, :] + deltar) \
        * (-0.5 * timereswi * efacp
            * (Vpp[:, 1:jf-1, :, :] - Vmp[:, 1:jf-1, :, :])
            + uekp[:, 1:jf-1, :, :]) \
        / dMdrp[:, 1:jf-1, :, :]
    um[:, 1:jf-1, :, :] = (rf[:, 1:jf-1, :, :]-deltar) \
        * (-0.5 * timereswi * efacm
           * (Vpm[:, 1:jf-1, :, :] - Vmm[:, 1:jf-1, :, :])
           + uekm[:, 1:jf-1, :, :]) \
        / dMdrm[:, 1:jf-1, :, :]

    #  Do not calculate vertical velocities if either radial velocity is zero.
    #  This is necessary because secondary wind maximum can vanish from one time
    #  step to the next, so that time rate of change blows up
    ufac = np.minimum(np.abs(30 * up[:, 1:jf-1, :, :]), 1) \
        * np.minimum(np.abs(30 * um[:, 1:jf-1, :, :]), 1)
    w2[:, 1:jf-1, :, :] = -Htrop * ufac * deltari \
        * ((rf[:, 1:jf-1, :, :] + deltar) * up[:, 1:jf-1, :, :]
           - (rf[:, 1:jf-1, :, :] - deltar)*um[:, 1:jf-1, :, :]) \
        / np.maximum(rf[:, 1:jf-1, :, :], 1)

    return w2, cp, V, Vd, Vrp, Vrm


def windprofiles(vm, rm, r, wp, vm2=None, rm2=None, opt=True):
    """
    Calculate the radial profiles of azimuthal wind. If secondary wind parameters
    (vm2 and rm2) are provided, the function incorporates any secondary eyewalls.

    Inputs:
    -----------
    vm : float
        Maximum circular wind speed (knots).
    rm : float
        Radius of maximum wind (km).
    r : float
        Distance of each event from the point of interest (POI) (km).
    wp : float
        Wind profile shape parameter, typically used to adjust the wind decay rate.
    vm2 : float, optional
        Secondary maximum circular wind speed (knots), if a secondary eyewall is present.
        Default is 0, indicating no secondary eyewall.
    rm2 : float, optional
        Secondary radius of maximum wind (km). Default is 0, indicating no secondary eyewall.

    Returns:
    --------
    V : float
        Azimuthal wind (knots)

    Notes:
    ------
    - The function assumes vm2 and rm2 are set to 0 if no secondary eyewall is present.
    """

    wprofile = wp  # Use holland (1) or emanuel (2) or er2011 (3) wind profile
    wprofile = np.min([wprofile, 3])  # Forbid use of experimental profiles
    vm = vm * 1852 / 3600  # Convert maximum wind speed to m/s

    if vm2 is not None:
        vm2 = vm2 * 1852 / 3600     # Convert maximum wind speed to m/s

    if rm2 is not None:
        se = np.sum(rm2[rm2 != 0])  # Test if there are any secondary eyewalls

    # Holland 2010 wind profile model
    if wprofile == 1:
        bs = 1.8
        rn = 300 / 20  # Holland parameters
        xn = 0.8

        rh = r / rm
        x = 0.5 + (rh - 1) * (xn - 0.5) / (rn - 1)
        x[x < 0.5] = 0.5
        V = vm * (rh ^ -bs * np.exp(1 - rh ^ -bs)) ** x

        if rm2 is not None:
            if se != 0:
                rh = r / rm2
                x = 0.5 + (rh - 1) * (xn - 0.5) / (rn - 1)
                x[x < 0.5] = 0.5
                V2 = vm2 * (rh**-bs * np.exp(1 - rh**-bs)) ** x

    elif wprofile == 2:
        r0 = 1000  # Outer radius (km)

        # Re-scale radius of maximum winds by random number drawn from log-normal distribution
        # Shape parameters
        b = 0.25
        nb = 0.9
        mb = 1.6
        mb2 = 2 * mb
        fac1 = (1 - b) * (mb + nb)
        fac2 = b * (1 + 2 * mb)
        fac3 = 2 * (mb + nb)
        fac4 = 2 * mb + 1
        rat = r / np.maximum(rm, 1)
        V = (
            vm
            * (np.maximum((r0 - r), 0) / (r0 - rm))
            * np.sqrt(
                rat**mb2 * (fac1 / (nb + mb * rat**fac3) + fac2 / (1 + mb2 * rat**fac4))
            )
        )
        if rm2 is not None:
            if se != 0:
                rat = r / np.maximum(rm2, 1)
                V2 = (
                    vm2
                    * (np.maximum((r0 - r), 0) / (r0 - rm2))
                    * np.sqrt(rat**mb2 * (fac1 / (nb + mb * rat**fac3) + fac2 / (1 + mb2 * rat**fac4)))
                )

    elif wprofile == 3:
        crat = 1
        f = 5.0e-5 * 1000        # converts from kms to meters

        if opt:
            Mm = rm * vm
        else:
            Mm = rm * vm + 0.5 * f * rm**2
        rn = r / np.where(rm == 0, 1e-32, rm)

        if crat == 1:
            M = Mm * (2 * rn**2 / (1 + rn**2))
        else:
            M = Mm * (2 * rn**2 / (2 - crat + crat * rn**2)) ** (1 / (2 - crat))

        # Add long tail to V to avoid discontinuity in vorticity at outer radius
        V = M / (r + 1e-8)
        V[V < 0] = 0
        V = np.nan_to_num(V)

        if rm2 is not None:
            if se != 0:
                Mm = rm2 * vm2
                rn = r / rm2
                if crat == 1:
                    M = Mm * (2 * rn**2 / (1 + rn**2))
                else:
                    M = Mm * (2 * rn**2 / (2 - crat + crat * rn**2)) ** (1 / (2.0 - crat))

                # Add long tail to V to avoid discontinuity in vorticity at outer radius
                V2 = M / (r + 1e-8)
                V2 = np.maximum(V2, 0)
                V2 = np.nan_to_num(V2)

    # Merge primary and secondary wind profiles
    if rm2 is not None:
        if se != 0 & wprofile < 4:
            rm_clipped = np.where(rm == 0, 1e-10, rm)
            u = np.maximum(r, 1) / rm_clipped
            hu = np.maximum(np.sign(u - 1), 0)
            V = V * (1 - hu) + hu * vm / u
            u = r / np.maximum(rm2, 1)
            hu = np.maximum(np.sign(u - 1), 0)
            V2 = V2 * hu + (1 - hu) * vm2 * u
            V = np.maximum(V, V2)

    V = V * 3600 / 1852    # Convert wind speed to knots

    return V


def transfunction(latitude):
    """
    This function produces a translation factor (transfactor) based on latitude,
    which is added to the circular wind speeds of the synthetic track sets.
    """
    transfac = 0.8
    transcap = 1
    amplitude = 0.35
    centerlat = 35
    latscale = 10

    # Calculate the translation factor
    transfactor = transfac + amplitude * (1 + np.tanh((np.abs(latitude) - centerlat) / latscale))

    # Cap the translation factor at the transcap value
    transfactor = np.minimum(transfactor, transcap)

    return transfactor


def utrans(latstore, longstore):
    """
    This function reads in track file, calculate translation speeds (knots) and smooth

    Inputs:
        - latstore [n x m]: Latitude of points along each track
        - longstore [n x m]: Longitude of points along each track

    Returns:
        - ut, vt [nn x m]: west-east, north-south component of the storm translation velocity
    """
    smfac = 0.4                         # Smoothing factor
    pifac = math.acos(-1)/180           # pi number
    dfac = 60 * 1.852                   # 1 nautical mile = 1/60 degree = 1.852 km
    dtfac = 1/(4*1.852)
    netfac = dtfac * dfac
    nn, m = np.shape(latstore)          # nn: number of storms; m: max number of time points
    ut = np.nan*np.zeros((nn, m))
    vt = np.nan*np.zeros((nn, m))
    jmax = np.zeros((nn))               # time length of each event

    # loop for each storm
    for n in range(nn):
        lat = latstore[n, :]
        long = longstore[n, :]

        # Get the length of each event by finding the first 0 element
        lat = lat[(lat != 0) & ~np.isnan(lat)]
        jm = len(lat)
        jm = np.maximum(jm, 1)
        jmax[n] = jm

        # Calculate storm translation velocity.
        # The difference of lat/long between 2 consecutive time steps provide the distance that
        # storms move during this time period, thus the translation speed in x & y direction.

        # ...for spatial difference in longitude, need to handle when track cross the prime meridian
        longdif = long[2:jm] - long[0:jm-2]
        for j in range(jm-2):
            if longdif[j] < -300:
                longdif[j] = longdif[j]+360
            elif longdif[j] > 300:
                longdif[j] = longdif[j]-360

        # include correction for ut at high latitude
        ut[n, 1: jm - 1] = netfac * np.cos(pifac * lat[1: jm - 1]) * longdif
        vt[n, 1: jm - 1] = netfac * (latstore[n, 2:jm] - latstore[n, 0: jm - 2])
        j2 = np.minimum(2, jm)
        j3 = np.maximum(jm - 3, 0)
        ut[n, 0] = 2 * ut[n, 1] - ut[n, j2]
        vt[n, 0] = 2 * vt[n, 1] - vt[n, j2]
        if jm > 0:
            ut[n, jm - 1] = 2.0 * ut[n, jm - 2] - ut[n, j3]
            vt[n, jm - 1] = 2.0 * vt[n, jm - 2] - vt[n, j3]

        # Smooth translation velocity
        vt[n, 1 : jm - 1] = vt[n, 1 : jm - 1] + smfac * (
            vt[n, 0 : jm - 2] + vt[n, 2:jm] - 2 * vt[n, 1 : jm - 1]
        )
        ut[n, 1 : jm - 1] = ut[n, 1 : jm - 1] + smfac * (
            ut[n, 0 : jm - 2] + ut[n, 2:jm] - 2 * ut[n, 1 : jm - 1]
        )

    ut = np.nan_to_num(ut)
    vt = np.nan_to_num(vt)

    return ut, vt, jmax


def utransfull(latstore, longstore, vstore=None, u850store=None, v850store=None):
    """
    Calculates translation speed and corrections to circular wind speed.
    In addition to latstore and longstore, inputs (in order) can be vstore,
    u850store, and v850store if available.
    """
    # Calculate translation velocity
    ut, vt, _ = utrans(latstore, longstore)

    # Add latitude-dependent fraction of translation speed to circular wind
    transfunc = transfunction(latstore)
    uinc = transfunc * ut
    vinc = transfunc * vt

    # If environmental 850 hPa wind available, add baroclinic factor to wind speed
    if u850store is not None and v850store is not None:
        udrift = -0.9 * 3600 / 1852
        vdrift = 1.4 * 3600 / 1852
        vdrift = vdrift * latstore / (np.abs(latstore) + 1e-8)

        # Apply baroclinic correction (factors changed from 0.25 to 0.5 in May 2014)
        uinc = uinc + 0.5 * (ut - udrift - u850store) * (vstore / 15)
        vinc = vinc + 0.5 * (vt - vdrift - v850store) * (vstore / 15)

    # Do not allow wind speed increments to exceed circular gradient wind
    if vstore is not None:
        ufrac = vstore / (0.1 + np.abs(uinc))
        ufrac = np.minimum(ufrac, 1)
        uinc = uinc * ufrac

        vfrac = vstore / (0.1 + np.abs(vinc))
        vfrac = np.minimum(vfrac, 1)
        vinc = vinc * vfrac

    return ut, vt, uinc, vinc


def pointwfield(latstore, longstore, vstore, rmstore, vsestore,
                rmsestore, ut, vt, us, vs, plat, plong, h, hx, hy):
    """
    This function calculates the spatial distribution of vertical velocity.
    nn: number of storms; m: max number of time points (length of each storm)

    Inputs:
    ------
        - latstore[n x m], longstore[nxm]: Latitude and longitude of points along each track
        - vstore [n x m]: Maximum circular wind speed along each track.
        - vsestore [n x m]: The maximum circular wind speed of maximum circular wind of any
                secondary eyewalls that may be present, along each track.
        - rmstore [n x m]: The radius (km) of maximum circular wind along each track
        - rmsestore [n x m]: The radius (km) of maximum circular wind of any secondary wind
                maxima present (0 if absent), along each track.
        - ut, vt [nn x m]: west-east, north-south component of the storm translation velocity
        - us, vs [nn x m]: vertical shears (m/s) used to estimate the baroclinic components of the
                vertical motion
        - h, hx, hy: contain topographic heights and their derivatives in x and y
    Returns:
    ------
        - w [n x m]: spatial distribution of vertical velocity centered at plong, plat
    """

    # Load parameters
    deltar = params.deltar              # Delta radius (km)
    timeresw = params.timeresw          # Native time resolution of WRT output (km)
    Htrop = params.Htrop                # Depth of lower troposphere (m)
    radcity = params.radcity            # Distance of storm from point of interest beyond which
                                        # the influence of the storm is ignored (km)
    wprofile = params.wprofile          # Wind profile
    se = np.max(vsestore)               # Check for secondary eyewalls
    deltari = 1. / deltar
    timereswi = 1. / (3600 * timeresw)
    Hi = 1. / Htrop
    nn = 1
    m = len(ut)
    sx = len(plong)
    sy = len(plat)

    ngrid = 0
    if sx == sy:
        ngrid = 1

    if ngrid == 1:
        sy = 1
    else:
        sy = np.max(np.shape(plat))

    if plong[0] < 0:
        plong = plong + 360

    pifac = math.acos(-1)/180           # pi number
    dfac = 60 * 1.852                   # 1 nautical mile = 1/60 degree = 1.852 km
    sfac = 1 / (0.25 * 60 * 1852)       # factor converting 0.25 degree to m
    knotfac = 1852. / 3600              # 1 knots = 0.5144 m/s
    ut = ut * knotfac                   # convert knots to m/s
    vt = vt * knotfac                   # convert knots to m/s
    omega = math.acos(-1) / (6 * 3600)  # Earth angular velocity parameter
    if latstore.ndim == 1:
        latfac = latstore[0] / (abs(latstore[0]) + 1e-8)
    elif latstore.ndim == 2:
        latfac = latstore[0, 0] / (abs(latstore[0, 0]) + 1e-8)

    # Estimate Drag coeffs:
    cdrag, cdx, cdy = tcr_tb.estimate_drag_coefficients(plat, plong, sfac)
    cdrag = np.maximum(cdrag, 0.0015)   # To put lower bound on drag coefficeint over water

    rfull, dx, dy = tcr_tb.calculate_distance_to_track(
        plat, plong, latstore, longstore, nn, m, sx, sy, ngrid, dfac)
    rfull = np.maximum(rfull, 0.5)
    costhetafull = dx / rfull
    sinthetafull = dy / rfull

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
            vfull[0, :, i, j] = vstore[:]
            rmfull[0, :, i, j] = rmstore[:]
            vsefull[0, :, i, j] = vsestore[:]
            rmsefull[0, :, i, j] = rmsestore[:]
            latfull[0, :, i, j] = latstore[:]
            longfull[0, :, i, j] = longstore[:]
            utfull[0, :, i, j] = ut[:]
            vtfull[0, :, i, j] = vt[:]
            usfull[0, :, i, j] = us[:]
            vsfull[0, :, i, j] = vs[:]

    w, cp, V, Vd, Vrp, Vrm, u1temp, u2temp, Cdp, Cdm = \
        calculate_wind_primary(utfull, vtfull, vfull, rfull, rmfull, vsefull, rmsefull, latfull,
                               nn, 3, sx, sy, h, hx, hy, cdrag, hfull, hxfull, hyfull, cdfull, cdx,
                               cdxfull, cdy, cdyfull, Hi, Htrop, omega, pifac, deltar, deltari,
                               knotfac, latfac, timereswi, costhetafull, sinthetafull, 2, wprofile,
                               adj_water=True)

    # If secondary eyewalls present, add in their contribution to wind
    if se > 0:
        w2, cp, V, Vd, Vrp, Vrm = calculate_wind_secondary(rfull, vsefull, rmsefull, latfull, nn, 3,
                                                           sx, sy, Hi, Htrop, omega, pifac, deltar,
                                                           deltari, knotfac, latfac, timereswi,
                                                           u1temp, u2temp, Cdp, Cdm, wprofile)
        w = np.maximum(w, w2)

    # to include shear dotted with storm entropy
    hxmod = -0.0005 * (cp + 2 * V / (0.1+rfull) + deltari * (Vrp-Vrm)) * vsfull
    hymod = 0.0005 * (cp + 2 * V / (0.1+rfull) + deltari * (Vrp-Vrm)) * usfull
    ufunc = (radcity - rfull) / 50
    ufunc = np.minimum(ufunc, 1)
    ufunc = np.maximum(ufunc, 0)
    utfull = utfull * ufunc
    vtfull = vtfull * ufunc
    w[:, 1, :, :] = w[:, 1, :, :] \
        + (vtfull[:, 1, :, :] + Vd[:, 1, :, :] * latfac * costhetafull[:, 1, :, :]) \
        * hyfull[:, 1, :, :] \
        + (utfull[:, 1, :, :] - Vd[:, 1, :, :] * latfac * sinthetafull[:, 1, :, :]) \
        * hxfull[:, 1, :, :] \
        + Vd[:, 1, :, :] * costhetafull[:, 1, :, :] * hymod[:, 1, :, :] \
        - Vd[:, 1, :, :] * sinthetafull[:, 1, :, :] * hxmod[:, 1, :, :]
    w = np.minimum(w, 7)
    return w


def pointshortn(latstore, longstore, vstore, rmstore, vsestore, rmsestore, ut, vt,
                plat, plong, timeres):
    """
    Calculate time series of wind speed and direction at the points of interest.

    Inputs:
    ------
        - latstore, longstore: Latitudes and longitude along each track
        - vstore: maximum circuslar wind along each track
        - rmstore: the radius (km) of maximum circular wind of points along each track
        - vsestore: maximum circular wind of any secondary eyewalls that may be present
        - rmsestore: the radius (km) of maximum circular wind of any secondary eyewalls
        - ut, vt: west-east, north-south component of the storm translation velocity
        - plat: latitude of point of interest
        - plong: longitude of point of interest
        - timeres: time resolution

    Returns:
    -------
        - vs: wind speed
        - direction: wind direction
    """

    timelength = params.timelength
    wheight = params.wheight
    wprofile = params.wprofile          # Wind profile
    radcity = params.radcity

    nsteps = round(2 / timeres)
    nstepsi = 1 / nsteps
    delj = np.floor(timelength / 4).astype(int)

    nn, m = ut.shape
    sx = plong.shape[0]
    sy = plat.shape[0]
    ngrid = 0
    if sx == sy:
        ngrid = 1
        sy = 1

    if plong[0] < 0:
        plong = plong + 360

    # pifac = np.arccos(-1) / 180
    dfac = 60 * 1.852

    logfac = np.log(wheight / 500) / 0.35

    # Load bathymetry
    bathy = tcr_io.load_netcdf_2d_parameters('data', 'surface_data.nc', 'bathy')

    # Load neutral drag coefficients
    cd = tcr_io.load_netcdf_2d_parameters('data', 'surface_data.nc', 'cd')
    mincd = np.min(cd)
    cd[bathy < 0] = mincd
    rat = 1 / (1 + np.sqrt(mincd) * logfac)

    # Interpolate drag coefficient to points of interest (POI)
    cdrag = np.zeros((sx, sy))
    for i in range(sx):
        for j in range(sy):
            ib = np.floor(4 * plong[i]).astype(int)
            ibp = ib + 1
            if ibp > 1440 - 1:
                ibp = 0
            jb = np.floor(4 * (plat[j] + 90)).astype(int)
            b1 = cd[ib, jb]
            b2 = cd[ib, jb + 1]
            b3 = cd[ibp, jb]
            b4 = cd[ibp, jb + 1]
            dely = 4 * (plat[j] + 90) - jb
            delx = 4 * plong[i] - ib
            d1, d2 = (1 - delx) * (1 - dely), dely * (1 - delx)
            d3, d4 = delx * (1 - dely), delx * dely
            cdrag[i, j] = 1 / (d1 / b1 + d2 / b2 + d3 / b3 + d4 / b4)

    # Calculate distance of each POI from track
    radius, _, _ = tcr_tb.calculate_distance_to_track(
        plat, plong, latstore, longstore, nn, m, sx, sy, 0, dfac)
    radius = np.maximum(radius, 0.5)

    jmin = np.argmin(radius, axis=1)            # index where radius of the storm is smallest
    jmin = np.maximum(jmin, 0 + delj)           # cut off delj steps at the begining
    jmin = np.minimum(jmin, m - 1 - delj)       # cut off delj steps at the end
    jstart = jmin - delj                        # index start
    jend = jmin + delj + 1                      # index end
    jtot = 2 * delj + 1
    jfine = 1 + nsteps * (jtot - 1)

    # Create reduced length time series of each quantity
    vshort = np.zeros((nn, jtot, sx, sy))
    rmshort = np.zeros((nn, jtot, sx, sy))
    vseshort = np.zeros((nn, jtot, sx, sy))
    rmseshort = np.zeros((nn, jtot, sx, sy))
    latshort = np.zeros((nn, jtot, sx, sy))
    longshort = np.zeros((nn, jtot, sx, sy))
    utshort = np.zeros((nn, jtot, sx, sy))
    vtshort = np.zeros((nn, jtot, sx, sy))

    for i in range(sx):
        for j in range(sy):
            for n in range(nn):
                vshort[n, :, i, j] = vstore[n, jstart[n, i, j]: jend[n, i, j]]
                rmshort[n, :, i, j] = rmstore[n, jstart[n, i, j]: jend[n, i, j]]
                vseshort[n, :, i, j] = vsestore[n, jstart[n, i, j]: jend[n, i, j]]
                rmseshort[n, :, i, j] = rmsestore[n, jstart[n, i, j]: jend[n, i, j]]
                latshort[n, :, i, j] = latstore[n, jstart[n, i, j]: jend[n, i, j]]
                longshort[n, :, i, j] = longstore[n, jstart[n, i, j]: jend[n, i, j]]
                utshort[n, :, i, j] = ut[n, jstart[n, i, j]: jend[n, i, j]]
                vtshort[n, :, i, j] = vt[n, jstart[n, i, j]: jend[n, i, j]]

    # Create high time-resolution series
    vfine = np.zeros((nn, jfine, sx, sy))
    rmfine = np.zeros((nn, jfine, sx, sy))
    vsefine = np.zeros((nn, jfine, sx, sy))
    rmsefine = np.zeros((nn, jfine, sx, sy))
    latfine = np.zeros((nn, jfine, sx, sy))
    longfine = np.zeros((nn, jfine, sx, sy))
    utfine = np.zeros((nn, jfine, sx, sy))
    vtfine = np.zeros((nn, jfine, sx, sy))

    k = 0
    for j in range(jtot - 1):
        for n in range(nsteps):
            weight = n * nstepsi
            vfine[:, k, :, :] = (1 - weight) * vshort[:, j, :, :] + weight * vshort[:, j + 1, :, :]
            rmfine[:, k, :, :] = (1 - weight) * rmshort[:, j, :, :] + \
                weight * rmshort[:, j + 1, :, :]
            vsefine[:, k, :, :] = (1 - weight) * vseshort[:, j, :, :] + \
                weight * vseshort[:, j + 1, :, :]
            rmsefine[:, k, :, :] = (1 - weight) * rmseshort[:, j, :, :] + \
                weight * rmseshort[:, j + 1, :, :]
            latfine[:, k, :, :] = (1 - weight) * latshort[:, j, :, :] + \
                weight * latshort[:, j + 1, :, :]
            longfine[:, k, :, :] = (1 - weight) * longshort[:, j, :, :] + \
                weight * longshort[:, j + 1, :, :]
            utfine[:, k, :, :] = (1 - weight) * utshort[:, j, :, :] + \
                weight * utshort[:, j + 1, :, :]
            vtfine[:, k, :, :] = (1 - weight) * vtshort[:, j, :, :] + \
                weight * vtshort[:, j + 1, :, :]
            k += 1

    vfine[:, k, :, :] = vshort[:, jtot - 1, :, :]
    rmfine[:, k, :, :] = rmshort[:, jtot - 1, :, :]
    vsefine[:, k, :, :] = vseshort[:, jtot - 1, :, :]
    rmsefine[:, k, :, :] = rmseshort[:, jtot - 1, :, :]
    latfine[:, k, :, :] = latshort[:, jtot - 1, :, :]
    longfine[:, k, :, :] = longshort[:, jtot - 1, :, :]
    utfine[:, k, :, :] = utshort[:, jtot - 1, :, :]
    vtfine[:, k, :, :] = vtshort[:, jtot - 1, :, :]

    rfine, dx, dy = tcr_tb.calculate_distance_to_track(
        plat, plong, latfine, longfine, nn, jfine, sx, sy, ngrid, dfac)

    V = windprofiles(vfine, rmfine, rfine, wprofile, vsefine, rmsefine, opt=True)
    V = V * latfine / (np.abs(latfine)+1e-8)

    # Calculate cdfac and update V
    cdfine = np.tile(cdrag, (nn, jfine, 1, 1))
    cdfac = np.maximum(1 + np.sqrt(cdfine) * logfac, 0)
    V = V * rat * cdfac

    vn = vtfine + V * dx / np.maximum(rfine, 0.5)
    un = utfine - V * dy / np.maximum(rfine, 0.5)
    rfac = rfine / np.sqrt(rfine**2 + radcity**2)
    vs = np.sqrt(un**2 + vn**2) - rfac * np.sqrt(utfine**2 + vtfine**2)
    vs = np.maximum(vs, 0)
    tempd = 360 + np.degrees(np.arctan2(-un, -vn))
    direction = np.mod(tempd, 360)

    return vs, direction


def pointwshortn(latstore, longstore, vstore, rmstore, vsestore, rmsestore, ut,
                 vt, us, vs, plong, plat, h, hx, hy, timeres):
    """
    This function calculates time series of vertical velocity at the points of interest.

    Inputs:
    ------
        - latstore, longstore: Latitudes and longitude along each track
        - vstore: maximum circuslar wind along each track
        - rmstore: the radius (km) of maximum circular wind of points along each track
        - vsestore: maximum circular wind of any secondary eyewalls that may be present
        - rmsestore: the radius (km) of maximum circular wind of any secondary eyewalls
        - ut, vt: west-east, north-south component of the storm translation velocity
        - us, vs: vertical shears (m/s) used to estimate the baroclinic components of
                the vertical motion
        - plat: latitude of point of interest
        - plong: longitude of point of interest
        - h: topographic heights
        - hx, hy: derivatives of h in x and y
        - timeres: time resolution

    Returns:
    -------
        - wq: vertical velocity (m/s)
        - date_record: Time in date format corresponding to rainrate
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
    omega = math.acos(-1)/(6*3600)      # Earth angular velocity parameter

    se = np.max(vsestore)               # for secondary eyewalls
    ntime = int(timelength / timeres + 1)
    nsteps = round(2.0 / timeres)
    nstepsi = 1.0 / nsteps
    delj = np.floor(timelength/nsteps).astype(int)
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

    latfac = latstore[0, 0] / (abs(latstore[0, 0])+1e-8)
    ut = ut * knotfac
    vt = vt * knotfac

    cdrag, cdx, cdy = tcr_tb.estimate_drag_coefficients(plat, plong, sfac)

    # To put lower bound on drag coefficient over water
    # cdrag = np.maximum(cdrag, 0.0015)

    # Reduce drag coefficient of near-coastal locations
    cdrag = np.minimum(cdrag, 1.5e-3 * (1 + np.maximum(h, 0) / 100))
    cdx = cdx * 0.01 * np.minimum(np.maximum(h, 0), 100)
    cdy = cdy * 0.01 * np.minimum(np.maximum(h, 0), 100)

    radius, dx, dy = tcr_tb.calculate_distance_to_track(
        plat, plong, latstore, longstore, nn, m, sx, sy, ngrid, dfac)
    radius = np.maximum(radius, 0.5)

    # Get index where the radius is smallest to start
    jmin = np.argmin(radius, axis=1)
    jmin = np.maximum(jmin, 0+delj)
    jmin = np.minimum(jmin, m-1-delj)
    jstart = jmin - delj
    jend = jmin + delj + 1
    jtot = 2 * delj + 1
    jfine = 1 + nsteps * (jtot - 1)

    # Create arrays for reduced length time series of each quantity
    vshort = np.zeros((nn, jtot, sx, sy))
    rmshort = np.zeros((nn, jtot, sx, sy))
    vseshort = np.zeros((nn, jtot, sx, sy))
    rmseshort = np.zeros((nn, jtot, sx, sy))
    rshort = np.zeros((nn, jtot, sx, sy))
    latshort = np.zeros((nn, jtot, sx, sy))
    longshort = np.zeros((nn, jtot, sx, sy))
    utshort = np.zeros((nn, jtot, sx, sy))
    vtshort = np.zeros((nn, jtot, sx, sy))
    usshort = np.zeros((nn, jtot, sx, sy))
    vsshort = np.zeros((nn, jtot, sx, sy))

    for i in range(sx):
        for j in range(sy):
            for n in range(nn):
                vshort[n, :, i, j] = vstore[n, jstart[n, i, j]:jend[n, i, j]]
                rmshort[n, :, i, j] = rmstore[n, jstart[n, i, j]:jend[n, i, j]]
                vseshort[n, :, i, j] = vsestore[n, jstart[n, i, j]:jend[n, i, j]]
                rmseshort[n, :, i, j] = rmsestore[n, jstart[n, i, j]:jend[n, i, j]]
                rshort[n, :, i, j] = radius[n, jstart[n, i, j]:jend[n, i, j], i, j]
                latshort[n, :, i, j] = latstore[n, jstart[n, i, j]:jend[n, i, j]]
                longshort[n, :, i, j] = longstore[n,  jstart[n, i, j]:jend[n, i, j]]
                utshort[n, :, i, j] = ut[n, jstart[n, i, j]:jend[n, i, j]]
                vtshort[n, :, i, j] = vt[n, jstart[n, i, j]:jend[n, i, j]]
                usshort[n, :, i, j] = us[n, jstart[n, i, j]:jend[n, i, j]]
                vsshort[n, :, i, j] = vs[n, jstart[n, i, j]:jend[n, i, j]]

    # Create arrays filled with zeros
    vfine = np.zeros((nn, jfine, sx, sy))
    rmfine = np.zeros((nn, jfine, sx, sy))
    vsefine = np.zeros((nn, jfine, sx, sy))
    rmsefine = np.zeros((nn, jfine, sx, sy))
    latfine = np.zeros((nn, jfine, sx, sy))
    longfine = np.zeros((nn, jfine, sx, sy))
    utfine = np.zeros((nn, jfine, sx, sy))
    vtfine = np.zeros((nn, jfine, sx, sy))
    usfine = np.zeros((nn, jfine, sx, sy))
    vsfine = np.zeros((nn, jfine, sx, sy))
    hfine = np.zeros((nn, jfine, sx, sy))
    hyfine = np.zeros((nn, jfine, sx, sy))
    hxfine = np.zeros((nn, jfine, sx, sy))
    cdfine = np.zeros((nn, jfine, sx, sy))
    cdyfine = np.zeros((nn, jfine, sx, sy))
    cdxfine = np.zeros((nn, jfine, sx, sy))

    k = 0
    for j in range(jtot - 1):
        for n in range(nsteps):
            weight = n * nstepsi
            vfine[:, k, :, :] = (1 - weight) * vshort[:, j, :, :] \
                + weight * vshort[:, j + 1, :, :]
            rmfine[:, k, :, :] = (1 - weight) * rmshort[:, j, :, :] \
                + weight * rmshort[:, j + 1, :, :]
            vsefine[:, k, :, :] = (1 - weight) * vseshort[:, j, :, :] \
                + weight * vseshort[:, j + 1, :, :]
            rmsefine[:, k, :, :] = (1 - weight) * rmseshort[:, j, :, :] \
                + weight * rmseshort[:, j + 1, :, :]
            latfine[:, k, :, :] = (1 - weight) * latshort[:, j, :, :] \
                + weight * latshort[:, j + 1, :, :]
            longfine[:, k, :, :] = (1 - weight) * longshort[:, j, :, :] \
                + weight * longshort[:, j + 1, :, :]
            utfine[:, k, :, :] = (1 - weight) * utshort[:, j, :, :] \
                + weight * utshort[:, j + 1, :, :]
            vtfine[:, k, :, :] = (1 - weight) * vtshort[:, j, :, :] \
                + weight * vtshort[:, j + 1, :, :]
            usfine[:, k, :, :] = (1 - weight) * usshort[:, j, :, :] \
                + weight * usshort[:, j + 1, :, :]
            vsfine[:, k, :, :] = (1 - weight) * vsshort[:, j, :, :] \
                + weight * vsshort[:, j + 1, :, :]
            k += 1

    vfine[:, k, :, :] = vshort[:, jtot - 1, :, :]
    rmfine[:, k, :, :] = rmshort[:, jtot - 1, :, :]
    vsefine[:, k, :, :] = vseshort[:, jtot - 1, :, :]
    rmsefine[:, k, :, :] = rmseshort[:, jtot - 1, :, :]
    latfine[:, k, :, :] = latshort[:, jtot - 1, :, :]
    longfine[:, k, :, :] = longshort[:, jtot - 1, :, :]
    utfine[:, k, :, :] = utshort[:, jtot - 1, :, :]
    vtfine[:, k, :, :] = vtshort[:, jtot - 1, :, :]
    usfine[:, k, :, :] = usshort[:, jtot - 1, :, :]
    vsfine[:, k, :, :] = vsshort[:, jtot - 1, :, :]

    rmsefine = np.maximum(rmsefine, 0.1)
    rfine, dx, dy = tcr_tb.calculate_distance_to_track(
        plat, plong, latfine, longfine, nn, jfine, sx, sy, ngrid, dfac)
    rfinei = 1 / np.maximum(rfine, 1)

    w, cp, V, Vd, Vrp, Vrm, u1temp, u2temp, Cdp, Cdm = calculate_wind_primary(
        utfine, vtfine, vfine, rfine, rmfine, vsefine, rmsefine, latfine,
        nn, jfine, sx, sy, h, hx, hy, cdrag, hfine, hxfine,
        hyfine, cdfine, cdx, cdxfine, cdy, cdyfine, Hi, Htrop, omega, pifac,
        deltar, deltari, knotfac, latfac, timeresi, dx*rfinei, dy*rfinei, 10,
        wprofile, adj_water=False)

    if se > 0:  # If secondary eyewalls present, add in their contribution to wind
        w2, cp, V, Vd, Vrp, Vrm = calculate_wind_secondary(
            rfine, vsefine, rmsefine, latfine, nn, jfine, sx, sy, Hi, Htrop, omega,
            pifac, deltar, deltari, knotfac, latfac, timeresi, u1temp, u2temp, Cdp,
            Cdm, wprofile)
        w = np.maximum(w, w2)

    # Now add in topographic and shear components
    # to include shear dotted with storm entropy
    hxmod = -0.0005 * (cp + 2 * V / (0.1 + rfine) + deltari * (Vrp - Vrm)) * vsfine
    hymod = 0.0005 * (cp + 2 * V / (0.1 + rfine) + deltari * (Vrp - Vrm)) * usfine

    # Reduce effect of translation speed outside of radcity
    ufunc = (radcity - rfine) / 50.0
    ufunc = np.minimum(ufunc, 1)
    ufunc = np.maximum(ufunc, 0)
    utfine = utfine * ufunc
    vtfine = vtfine * ufunc

    # Reduce effect of orography outside of storm core
    ufunc = (150 - rfine) / 30
    ufunc = np.minimum(ufunc, 0.6)
    ufunc = np.maximum(ufunc, 0.2)
    hxfine = hxfine * ufunc
    hyfine = hyfine * ufunc

    w = (
        w
        + (Vd * latfac * dx * rfinei + vtfine) * hyfine
        + (utfine - Vd * latfac * dy * rfinei) * hxfine
        + Vd * dx * rfinei * hymod
        - Vd * dy * rfinei * hxmod
    )
    w = np.minimum(w, 7)
    return w


def pointwshortnqdx(latstore, longstore, datestore, dq, vstore, rmstore, vsestore, rmsestore,
                    ut, vt, us, vs, plong, plat, h, hx, hy, timeres, wrad):
    """
    This function calculates time series of vertical velocity at the points of interest.

    Inputs:
    ------
        - latstore, longstore: Latitudes and longitude along each track
        - vstore: maximum circuslar wind along each track
        - rmstore: the radius (km) of maximum circular wind of points along each track
        - vsestore: maximum circular wind of any secondary eyewalls that may be present
        - rmsestore: the radius (km) of maximum circular wind of any secondary eyewalls
        - ut, vt: west-east, north-south component of the storm translation velocity
        - us, vs: vertical shears (m/s) used to estimate the baroclinic components of
                the vertical motion
        - plat: latitude of point of interest
        - plong: longitude of point of interest
        - h: topographic heights
        - hx, hy: gradients of h in x and y
        - timeres: time resolution (hr)
        - wrad: background subsidence velocity under radiative cooling

    Returns:
    -------
        - wq: vertical velocity (m/s)
        - date_record: Time in date format corresponding to rainrate
    """

    deltar = params.deltar      # Delta radius (km) for calculating dM/dr
    deltari = 1.0 / deltar      # inverse of Delta radius
    timelength = params.timelength
    timeresi = 1.0 / (3600 * timeres)
    Htrop = params.Htrop
    wprofile = params.wprofile
    radcity = params.radcity

    se = np.nanmax(vsestore)       # for secondary eyewalls
    ntime = int(timelength / timeres + 1)
    nsteps = round(2.0 / timeres)
    nstepsi = 1.0 / nsteps
    delj = np.floor(timelength / 4).astype(int)

    Hi = 1.0 / Htrop

    nn, m = ut.shape
    sx = plong.size
    sy = 1
    w = np.zeros((nn, ntime, sx, sy))

    knotfac = 1852.0 / 3600             # convert knots to m/s (1 knots = 0.5144 m/s)
    ut = ut * knotfac                   # in m/s
    vt = vt * knotfac                   # in m/s

    dfac = 60 * 1.852                   # convert degree to km (1 nautical = 1/60 degree = 1.852 km)
    sfac = 1 / (0.25 * dfac)
    pifac = math.acos(-1) / 180         # pi number
    latfac = latstore[0, 0] / (abs(latstore[0, 0]) + 1e-8)
    omega = math.acos(-1) / (6 * 3600)  # Earth angular velocity parameter

    if plong[0] < 0:
        plong = plong + 360

    # Estimate Drag coeffs:
    cdrag, cdx, cdy = tcr_tb.estimate_drag_coefficients(plat, plong, sfac)
    cdrag = np.minimum(cdrag, 1.5e-3 * (1 + np.maximum(h, 0) / 100))
    cdx = cdx * 0.01 * np.minimum(np.maximum(h, 0), 100)
    cdy = cdy * 0.01 * np.minimum(np.maximum(h, 0), 100)

    # Calculate radius distance from storm center to POI
    radius, dx, dy = tcr_tb.calculate_distance_to_track(
        plat, plong, latstore, longstore, nn, m, sx, sy, 0, dfac)
    radius = np.maximum(radius, 0.5)

    jmin = np.argmin(radius, axis=1)            # index where radius of the storm is smallest
    jmin = np.maximum(jmin, 0 + delj)           # cut off delj steps at the begining
    jmin = np.minimum(jmin, m - 1 - delj)       # cut off delj steps at the end
    jstart = jmin - delj                        # index start
    jend = jmin + delj + 1                      # index end
    jtot = 2 * delj + 1                         # total index?
    jfine = 1 + nsteps * (jtot - 1)             # ???

    # Create reduced length time series of each quantity
    vshort = np.zeros((nn, jtot, sx, sy))
    rmshort = np.zeros((nn, jtot, sx, sy))
    vseshort = np.zeros((nn, jtot, sx, sy))
    rmseshort = np.zeros((nn, jtot, sx, sy))
    rshort = np.zeros((nn, jtot, sx, sy))
    latshort = np.zeros((nn, jtot, sx, sy))
    longshort = np.zeros((nn, jtot, sx, sy))
    dateshort = np.zeros((nn, jtot, sx, sy))
    # dateshort = np.full((nn, jtot, sx, sy), "", dtype="U19")
    utshort = np.zeros((nn, jtot, sx, sy))
    vtshort = np.zeros((nn, jtot, sx, sy))
    usshort = np.zeros((nn, jtot, sx, sy))
    vsshort = np.zeros((nn, jtot, sx, sy))
    dqshort = np.zeros((nn, jtot, sx, sy))

    for i in range(sx):
        for j in range(sy):
            for n in range(nn):
                vshort[n, :, i, j] = vstore[n, jstart[n, i, j]: jend[n, i, j]]
                rmshort[n, :, i, j] = rmstore[n, jstart[n, i, j]: jend[n, i, j]]
                vseshort[n, :, i, j] = vsestore[n, jstart[n, i, j]: jend[n, i, j]]
                rmseshort[n, :, i, j] = rmsestore[n, jstart[n, i, j]: jend[n, i, j]]
                rshort[n, :, i, j] = radius[n, jstart[n, i, j]: jend[n, i, j], i, j]
                latshort[n, :, i, j] = latstore[n, jstart[n, i, j]: jend[n, i, j]]
                longshort[n, :, i, j] = longstore[n, jstart[n, i, j]: jend[n, i, j]]
                dateshort[n, :, i, j] = datestore[n, jstart[n, i, j]: jend[n, i, j]]
                utshort[n, :, i, j] = ut[n, jstart[n, i, j]: jend[n, i, j]]
                vtshort[n, :, i, j] = vt[n, jstart[n, i, j]: jend[n, i, j]]
                usshort[n, :, i, j] = us[n, jstart[n, i, j]: jend[n, i, j]]
                vsshort[n, :, i, j] = vs[n, jstart[n, i, j]: jend[n, i, j]]
                dqshort[n, :, i, j] = dq[n, jstart[n, i, j]:jend[n, i, j]]

    # Create high time-resolution series
    vfine = np.zeros((nn, jfine, sx, sy))
    rmfine = np.zeros((nn, jfine, sx, sy))
    vsefine = np.zeros((nn, jfine, sx, sy))
    rmsefine = np.zeros((nn, jfine, sx, sy))
    latfine = np.zeros((nn, jfine, sx, sy))
    longfine = np.zeros((nn, jfine, sx, sy))
    date_record = np.zeros((nn, jfine, sx, sy))
    utfine = np.zeros((nn, jfine, sx, sy))
    vtfine = np.zeros((nn, jfine, sx, sy))
    usfine = np.zeros((nn, jfine, sx, sy))
    vsfine = np.zeros((nn, jfine, sx, sy))
    dqfine = np.zeros((nn, jfine, sx, sy))    

    hfine = np.zeros((nn, jfine, sx, sy))
    hyfine = np.zeros((nn, jfine, sx, sy))
    hxfine = np.zeros((nn, jfine, sx, sy))
    cdfine = np.zeros((nn, jfine, sx, sy))
    cdyfine = np.zeros((nn, jfine, sx, sy))
    cdxfine = np.zeros((nn, jfine, sx, sy))    

    k = 0
    for j in range(jtot - 1):
        for n in range(nsteps):
            weight = n * nstepsi
            vfine[:, k, :, :] = (1 - weight) * vshort[:, j, :, :] + weight * vshort[:, j + 1, :, :]
            rmfine[:, k, :, :] = (1 - weight) * rmshort[:, j, :, :] + \
                weight * rmshort[:, j + 1, :, :]
            vsefine[:, k, :, :] = (1 - weight) * vseshort[:, j, :, :] + \
                weight * vseshort[:, j + 1, :, :]
            rmsefine[:, k, :, :] = (1 - weight) * rmseshort[:, j, :, :] + \
                weight * rmseshort[:, j + 1, :, :]
            latfine[:, k, :, :] = (1 - weight) * latshort[:, j, :, :] + \
                weight * latshort[:, j + 1, :, :]
            longfine[:, k, :, :] = (1 - weight) * longshort[:, j, :, :] + \
                weight * longshort[:, j + 1, :, :]
            date_record[:, k, :, :] = (1 - weight) * dateshort[:, j, :, :] + \
                weight * dateshort[:, j + 1, :, :]
            utfine[:, k, :, :] = (1 - weight) * utshort[:, j, :, :] + \
                weight * utshort[:, j + 1, :, :]
            vtfine[:, k, :, :] = (1 - weight) * vtshort[:, j, :, :] + \
                weight * vtshort[:, j + 1, :, :]
            usfine[:, k, :, :] = (1 - weight) * usshort[:, j, :, :] + \
                weight * usshort[:, j + 1, :, :]
            vsfine[:, k, :, :] = (1 - weight) * vsshort[:, j, :, :] + \
                weight * vsshort[:, j + 1, :, :]
            dqfine[:, k, :, :] = (1 - weight) * dqshort[:, j, :, :] + \
                weight * dqshort[:, j + 1, :, :]
            k += 1

    vfine[:, k, :, :] = vshort[:, jtot - 1, :, :]
    rmfine[:, k, :, :] = rmshort[:, jtot - 1, :, :]
    vsefine[:, k, :, :] = vseshort[:, jtot - 1, :, :]
    rmsefine[:, k, :, :] = rmseshort[:, jtot - 1, :, :]
    latfine[:, k, :, :] = latshort[:, jtot - 1, :, :]
    longfine[:, k, :, :] = longshort[:, jtot - 1, :, :]
    date_record[:, k, :, :] = dateshort[:, jtot - 1, :, :]
    utfine[:, k, :, :] = utshort[:, jtot - 1, :, :]
    vtfine[:, k, :, :] = vtshort[:, jtot - 1, :, :]
    usfine[:, k, :, :] = usshort[:, jtot - 1, :, :]
    vsfine[:, k, :, :] = vsshort[:, jtot - 1, :, :]
    dqfine[:, k, :, :] = dqshort[:, jtot - 1, :, :]

    rmsefine = np.maximum(rmsefine, 0.1)
    rfine, dx, dy = tcr_tb.calculate_distance_to_track(
        plat, plong, latfine, longfine, nn, jfine, sx, sy, 0, dfac)
    rfinei = 1 / np.maximum(rfine, 1)

    # for primary eyewalls
    w, cp, V, Vd, Vrp, Vrm, u1temp, u2temp, Cdp, Cdm = calculate_wind_primary(
        utfine, vtfine, vfine, rfine, rmfine, vsefine, rmsefine, latfine,
        nn, jfine, sx, sy, h, hx, hy, cdrag, hfine, hxfine,
        hyfine, cdfine, cdx, cdxfine, cdy, cdyfine, Hi, Htrop, omega, pifac,
        deltar, deltari, knotfac, latfac, timeresi, dx*rfinei, dy*rfinei, 10,
        wprofile, adj_water=False)

    # if secondary eyewalls present
    if se > 0:
        w2, cp, V, Vd, Vrp, Vrm = calculate_wind_secondary(
            rfine, vsefine, rmsefine, latfine, nn, jfine, sx, 1, Hi, Htrop, omega,
            pifac, deltar, deltari, knotfac, latfac, timeresi, u1temp, u2temp,
            Cdp, Cdm, wprofile)
        w = np.maximum(w, w2)

    # Now add in topographic and shear components
    hxmod = -0.0005 * (cp + 2 * V / (0.1 + rfine) + deltari * (Vrp - Vrm)) * vsfine
    # to include shear dotted with storm entropy
    hymod = 0.0005 * (cp + 2 * V / (0.1 + rfine) + deltari * (Vrp - Vrm)) * usfine

    # Reduce effect of translation speed outside of radcity
    ufunc = (radcity - rfine) / 50
    ufunc = np.minimum(ufunc, 1)
    ufunc = np.maximum(ufunc, 0)
    utfine = utfine * ufunc
    vtfine = vtfine * ufunc

    # Reduce effect of orography outside of storm core
    ufunc = (150 - rfine) / 30
    ufunc = np.minimum(ufunc, 0.6)
    ufunc = np.maximum(ufunc, 0.2)
    hxfine = hxfine * ufunc
    hyfine = hyfine * ufunc

    w = (
        w
        + (Vd * latfac * dx * rfinei + vtfine) * hyfine
        + (utfine - Vd * latfac * dy * rfinei) * hxfine
        + Vd * dx * rfinei * hymod
        - Vd * dy * rfinei * hxmod
    )
    w = np.minimum(w, 7)

    # Add radiative cooling
    wq = dqfine * np.maximum(w - wrad, 0)       # If w-wrad is negative (downward motion): wp = 0
    return wq, date_record


def vouternew(vm, fc, ro, wc, CD, q):
    """
    Numerical integration of the outer wind profile from simple ODE.

    Inputs:
    ------
        vm: maximum wind (m/s)
        f: Coriolis force (s^-1)
        ro: outer radiusin (km)
        wc: radiative subsidence rate (mm/s)
        CD: drag coefficient
        q: number of radial points

    Returns:
    -------
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
        m[i] = m[i+2] - 2.0 * dr * \
            (chi * m[i+1]**2 / (rond**2 - r[i+1]**2) - r[i+1])
        m[i+1] = m[i+1] + assl * (m[i] + m[i+2] - 2.0 * m[i+1])
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


def estimate_radius_wind(ds, lat_tracks, vmax_tracks, id_tracks, directory, fname):
    """
    Estimate the radius of maximum circular wind from maximum circular wind speed.

    Inputs:
    ------
        ds: xarray dataset containing data
        lat_tracks: latitude (degrees) for the wind tracks
        vmax_tracks: maximum circular wind speed (m/s)
        id_tracks: id of TC tracks
        directory: folder containing the ds dataset
        fname: filename of the dataset

    Returns:
    -------
        rm_tracks: radius of maximum circular wind (km)
    """

    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_folder = os.path.join(current_dir, directory)
    fullpathname = os.path.join(data_folder, fname)

    ro = 1000.          # Outer radius
    wc = 3.0            # Radiative subsidence rate
    cdouter = 1.2e-3    # Drag coefficient
    nouter = 1000       # Number of radial point
    num_tcs = len(id_tracks)

    rm_tracks = np.zeros(vmax_tracks.shape)
    for ind in range(num_tcs):
        vmax = vmax_tracks[ind, :]
        vmax = vmax[~np.isnan(vmax)]
        jmax = len(vmax)
        for jnd in range(jmax):
            alats = lat_tracks[ind, jnd]
            vsin = vmax_tracks[ind, jnd]
            # skip at equator (coriolis force fc = 0) or vs = 0
            if alats != 0 and vsin > 0:
                fc1 = 1.45e-4 * np.abs(alats * 0.0175)
                _, rm, _ = vouternew(vsin, fc1, ro, wc, cdouter, nouter)
            else:
                rm = 0
            rm_tracks[ind, jnd] = rm
        print(f"Tracks {ind}/{num_tcs} completed...", end='\r')
        sys.stdout.flush()
    print('\nDone!')

    ds = ds.assign(rm_tracks=(['n_trk', 'time'], rm_tracks))

    # remove old ds and save updated ds to disk
    if os.path.exists(fullpathname):
        os.remove(fullpathname)
    ds.to_netcdf(fullpathname, mode='w')

    return rm_tracks


def smoothb(x, nz, jmin, jmax):
    """
    Apply a 1-2-1 smoothing filter to a 2D array.

    Inputs:
    ------
        - x (numpy.ndarray): Input 2D array to be smoothed.
        - nz (int): Number of rows in the 2D array.
        - jmin (int): Minimum column index to start smoothing.
        - jmax (int): Maximum column index to end smoothing.

    Returns:
    -------
        - xsmooth (numpy.ndarray): Smoothed 2D array.
    """
    xsmooth = np.copy(x)  # Make a copy of the input array to avoid modifying it directly

    # Apply the 1-2-1 smoothing filter
    for i in range(1, nz):
        for j in range(jmin, jmax):
            xsmooth[i, j] = 0.125 * (x[i - 1, j] + x[i + 1, j] +
                                     x[i, j - 1] + x[i, j + 1]) + 0.5 * x[i, j]

    return xsmooth


def windswathx(nt, latstore, longstore, rmstore, vstore, rmsestore, vsestore, uinc, vinc):
    """
    This script calculates the distribution of maximum point wind speed (knots) for a single storm.

    Inputs:
    ------
        - nt: Track number of the storm
        - latstore, longstore: Latitudes and longitudes along each track
        - vstore: Maximum circular wind along each storm track
        - rmstore: Radius (in km) of maximum circular wind along each track
        - vsestore: maximum circular wind of any secondary eyewalls that may be present
        - rmsestore: Radius (in km) of maximum circular wind of any secondary eyewalls
        - ut: West-east component of the storm translation velocity
        - u850store, v850store: Zonal & meridional components of the 850 hPa environmental
            wind speed (knots)
    Returns:
    --------
        - x, y: vectors containing the longitudes and latitudes of the grid
        - maxwind: storm maximum wind speed (knots) at each point on the grid
    """
    magfac = params.magfac              # overall scale factor for storm size
    deltax = params.deltax              # longitudinal distance of map boundaries from storm center
    deltay = params.deltay              # latitudinal distance of map boundaries from storm center
    bxmin = params.bxmin                # minimum longitude of map (degree)
    bxmax = params.bxmax                # maximum longitude of map (degree)
    bymin = params.bymin                # minimum latitude of map (degree)
    bymax = params.bymax                # maximum latitude of map (degree)
    dellatlongs = params.dellatlongs    # horizontal resolution of field maps
    timeres = params.timeres            # time resolution for time series at fixed points

    rmstore1 = rmstore * magfac
    rmsestore1 = rmsestore * magfac

    # Find the index q where latstore is closest to 0
    q = np.argmin(np.abs(latstore[nt, :])) - 1

    # Extract non-zero elements
    utd = np.transpose(uinc[nt, :][uinc[nt, :] != 0])[np.newaxis, :]
    vtd = np.transpose(vinc[nt, :][vinc[nt, :] != 0])[np.newaxis, :]
    lat = np.transpose(latstore[nt, :][latstore[nt, :] != 0])[np.newaxis, :]
    long = np.transpose(longstore[nt, :][longstore[nt, :] != 0])[np.newaxis, :]
    v = np.transpose(vstore[nt, :][vstore[nt, :] != 0])[np.newaxis, :]

    qv = np.max(v.shape)
    rm = np.transpose(rmstore1[nt, :][rmstore1[nt, :] != 0])[np.newaxis, :]
    vse = np.transpose(vsestore[nt, :qv])[np.newaxis, :]
    rmse = np.transpose(rmsestore1[nt, :qv])[np.newaxis, :]

    for i in range(q):
        if long[0, 0] > 200 and long[0, i] < 50:
            long[0, i] += 360

    bxmin = np.min(long[np.nonzero(long)]) - deltax
    bxmax = np.max(long[np.nonzero(long)]) + deltax
    bymin = np.min(lat[np.nonzero(lat)]) - deltay
    bymax = np.max(lat[np.nonzero(lat)]) + deltay

    # Create x and y grid
    x = np.arange(bxmin, bxmax+1e-10, dellatlongs)
    y = np.arange(bymin, bymax+1e-10, dellatlongs)

    sx = x.size
    sy = y.size

    # Adjust x and y sizes if needed
    if sx == sy:
        x = np.append(x, bxmax + dellatlongs)
        sx += 1

    # Initialize maxwind array and compute maxwind
    maxwind = np.zeros((sx, sy))

    vs1, _ = pointshortn(lat, long, v, rm, vse, rmse, utd, vtd, y, x, timeres)
    maxwind[:, :] = np.max(vs1[0, :, :, :], axis=0)

    maxwind = smoothb(maxwind, sx - 1, 1, sy - 1)
    maxwind = np.transpose(maxwind)

    return x, y, maxwind
