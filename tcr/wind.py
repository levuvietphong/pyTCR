"""
Functions for wind in PyTCR
"""


import os
import sys
import math
import glob
import numpy as np
from tcr import terrain_boundary as tcr_tb
from tcr import iodata as tcr_io
from tcr import parameters as tcr_params
from tcr.datadir import BASE_DATA_DIR, DOWNSCALED_DATA_DIR


def calculate_wind_primary(utf, vtf, vf, rf, rmf, vsef, rmsef, latf,
                           nn, jf, sx, sy, h, hx, hy, cdrag, hf, hxf, hyf, cdf,
                           cdx, cdxf, cdy, cdyf, Hi, Htrop, omega, pifac,
                           deltar, deltari, knotfac, latfac, timereswi,
                           costhetaf, sinthetaf, thresM, wprofile,
                           adj_water=False):
    """
    Calculate the wind in the primary band.

    Parameters:
    -----------
    utf, vtf : array-like
        West-east and storm translation velocity (m/s)
    vf : array-like
        Maximum circular wind speed at each 2-hour point along each track (knots)
    rf : array-like
        Euclidean distance from each track (km)
    rmf : array-like
        Radius of maximum circular wind of 2-hour points along each track (km)
    vsef : array-like
        Maximum circular wind speed of any secondary eyewalls that may be
        present at each 2-hour point along each track (knots)
    rmsef : array-like
        Radius of maximum circular wind of any secondary wind maxima of
        2-hour points along each track (km) 
    latf : array-like
        Latitude of points (degree)
    nn, jf, sx, sy : int
        Dimensions of the input arrays
    h, hx, hy : array-like
        Topographic heights (m) and their derivatives in x and y
    cdrag : array-like
        Drag coefficient
    hf, hxf, hyf : array-like
        Topographic heights (m) and their derivatives in x and y for the fine grid
    cdf, cdx, cdxf, cdy, cdyf : array-like
        Drag coefficients and their derivatives
    Hi, Htrop : float
        Depth of the lower troposphere (m) and its inverse
    omega : float
        Earth angular velocity parameter (rad/s)
    pifac : float
        Pi factor
    deltar, deltari : float
        Delta radius (km) and its inverse
    knotfac : float
        Conversion factor from knots to m/s
    latfac : float
        Latitude factor
    timereswi : float
        Inverse of the native time resolution of WRT output in kilometers
    costhetaf, sinthetaf : array-like
        Cosine and sine of the angle between the storm translation velocity and
        the wind direction
    thresM : float
        Threshold for the gradient wind
    wprofile : array-like
        Wind profile parameter
    adj_water : bool, optional
        Adjust for roughness changes over water (default is False)

    Returns:
    --------
    w : array-like
        Vertical velocity (m/s)
    cp : array-like
        Coriolis parameter (1/s)
    V, Vd, Vrp, Vrm : array-like
        Wind profiles (knots)
    u1temp, u2temp : array-like
        Temporary variables for wind calculations
    Cdp, Cdm : array-like
        Drag coefficients
    """

    # Initialize arrays
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

    # Calculate wind profiles
    V = windprofiles(vf, rmf, rf, wprofile)
    Vd = windprofiles(vf, rmf, rf, wprofile, vsef, rmsef, opt=False)
    Vpp[:, 1:jf - 1, :, :] = windprofiles(
        vf[:, 2:jf, :, :], rmf[:, 2:jf, :, :], rf[:, 1:jf - 1, :, :] + deltar, wprofile)
    Vpm[:, 1:jf - 1, :, :] = windprofiles(
        vf[:, 2:jf, :, :], rmf[:, 2:jf, :, :],
        np.maximum(rf[:, 1:jf - 1, :, :] - deltar, 0), wprofile)
    Vmp[:, 1:jf - 1, :, :] = windprofiles(
        vf[:, 0:jf - 2, :, :], rmf[:, 0:jf - 2, :, :], rf[:, 1:jf - 1, :, :] + deltar, wprofile)
    Vmm[:, 1:jf - 1, :, :] = windprofiles(
        vf[:, 0:jf - 2, :, :],
        rmf[:, 0:jf - 2, :, :],
        np.maximum(rf[:, 1:jf - 1, :, :] - deltar, 0),
        wprofile
    )
    Vrp[:, :, :, :] = windprofiles(
        vf[:, :, :, :], rmf[:, :, :, :], rf[:, :, :, :] + deltar, wprofile)
    Vrm[:, :, :, :] = windprofiles(
        vf[:, :, :, :], rmf[:, :, :, :], np.maximum(rf[:, :, :, :] - deltar, 0), wprofile)

    # Convert knot to m/s
    V = knotfac * V
    Vd = knotfac * Vd
    Vpp = knotfac * Vpp
    Vpm = knotfac * Vpm
    Vmp = knotfac * Vmp
    Vmm = knotfac * Vmm
    Vrp = knotfac * Vrp
    Vrm = knotfac * Vrm

    # Calculate wind components
    vph = 0.5 * (V + Vrp)
    vmh = 0.5 * (V + Vrm)
    u1temp = vtf * costhetaf - utf * sinthetaf
    u2temp = vtf ** 2 + utf ** 2
    vnetp = np.sqrt(vph ** 2 + 2 * vph * latfac * u1temp + u2temp)
    vnetm = np.sqrt(vmh ** 2 + 2 * vmh * latfac * u1temp + u2temp)

    # Update topographic heights and drag coefficients
    hf[:] = h[np.newaxis, np.newaxis, :, :]
    hyf[:] = hy[np.newaxis, np.newaxis, :, :]
    hxf[:] = hx[np.newaxis, np.newaxis, :, :]
    cdf[:] = cdrag[np.newaxis, np.newaxis, :, :]
    cdyf[:] = cdy[np.newaxis, np.newaxis, :, :]
    cdxf[:] = cdx[np.newaxis, np.newaxis, :, :]

    # Calculate drag coefficients
    cdfac = 1e3 * 0.5 * deltar * (cdxf * costhetaf + cdyf * sinthetaf)
    Cdp = cdf + cdfac
    Cdm = cdf - cdfac
    Cdp = np.maximum(Cdp, 0)
    Cdp = np.minimum(Cdp, 0.005)
    Cdm = np.maximum(Cdm, 0)
    Cdm = np.minimum(Cdm, 0.005)

    # Adjust for roughness changes over water
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

    # Calculate Coriolis parameter and gradient wind
    cp = 1000 * 2 * omega * np.sin(pifac * abs(latf))
    dMdrp = (cp * (rf + 0.5 * deltar) +
             (rf + 0.5 * deltar) * deltari * (Vrp - V) +
             0.5 * (Vrp + V))
    dMdrp = np.maximum(dMdrp, thresM)
    dMdrm = (cp * (rf - 0.5 * deltar) +
             (rf - 0.5 * deltar) * deltari * (V - Vrm) +
             0.5 * (Vrm + V))
    dMdrm = np.maximum(dMdrm, thresM)

    rmf_safe = np.where(rmf[:, 1:jf - 1, :, :] == 0, 1e-32,
                        rmf[:, 1:jf - 1, :, :])
    efacp = np.minimum(
        (-1 + 2 * ((rf[:, 1:jf - 1, :, :] + deltar) / rmf_safe) ** 2), 1)
    efacm = np.minimum(
        (-1 + 2 * ((rf[:, 1:jf - 1, :, :] - deltar) / rmf_safe) ** 2), 1)

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


def calculate_wind_secondary(rf, vsef, rmsef, latf, nn, jf, sx, sy, Hi, Htrop,
                             omega, pifac, deltar, deltari, knotfac, latfac,
                             timereswi, u1temp, u2temp, Cdp, Cdm, wprofile):
    """
    Calculate wind for secondary eyewalls in tropical cyclones.
    This function computes various wind-related parameters for secondary
    eyewalls, including vertical velocity, wind profiles, and Coriolis
    parameters.

    Parameters:
    -----------
    rf : array-like
        Euclidean distance from each track (km)
    vsef : array-like
        Maximum circular wind speed of any secondary eyewalls (knots)
    rmsef : array-like
        Radius of maximum circular wind of any secondary wind maxima (km)
    latf : array-like
        Latitude of points (degree)
    nn, jf, sx, sy : int
        Dimensions of the input arrays.
    Hi : float
        Inverse of the depth of the lower troposphere (1/m)
    Htrop : float
        Depth of the lower troposphere (m)
    omega : float
        Earth angular velocity parameter (rad/s)
    pifac : float
        Pi factor
    deltar : float
        Delta radius (km)
    deltari : float
        Inverse of delta radius
    knotfac : float
        Conversion factor from knots to m/s
    latfac : float
        Latitude factor
    timereswi : float
        Inverse of the native time resolution of WRT output
    u1temp, u2temp : array-like
        Temporary variables for wind calculations
    Cdp, Cdm : array-like
        Drag coefficients
    wprofile : int
        Wind profile parameter

    Returns:
    --------
    w2 : array-like
        Vertical velocity for secondary eyewalls (m/s)
    cp : array-like
        Coriolis parameter (1/s)
    V, Vrp, Vrm : array-like
        Wind profiles for secondary eyewalls (knots)
    """

    # Initialize arrays
    V = np.zeros((nn, jf, sx, sy))
    Vpp = np.zeros((nn, jf, sx, sy))
    Vmp = np.zeros((nn, jf, sx, sy))
    Vpm = np.zeros((nn, jf, sx, sy))
    Vmm = np.zeros((nn, jf, sx, sy))
    Vrp = np.zeros((nn, jf, sx, sy))
    Vrm = np.zeros((nn, jf, sx, sy))
    up = np.zeros((nn, jf, sx, sy))
    um = np.zeros((nn, jf, sx, sy))
    w2 = np.zeros((nn, jf, sx, sy))

    # Calculate wind profiles
    V = windprofiles(vsef, rmsef, rf, wprofile)
    Vpp[:, 1:jf-1, :, :] = windprofiles(
        vsef[:, 2:jf, :, :],
        rmsef[:, 2:jf, :, :],
        rf[:, 1:jf-1, :, :] + deltar,
        wprofile
    )
    Vpm[:, 1:jf-1, :, :] = windprofiles(
        vsef[:, 2:jf, :, :],
        rmsef[:, 2:jf, :, :],
        np.maximum(rf[:, 1:jf-1, :, :] - deltar, 0),
        wprofile
    )
    Vmp[:, 1:jf-1, :, :] = windprofiles(
        vsef[:, 0:jf-2, :, :],
        rmsef[:, 0:jf-2, :, :],
        rf[:, 1:jf-1, :, :] + deltar,
        wprofile
    )
    Vmm[:, 1:jf-1, :, :] = windprofiles(
        vsef[:, 0:jf-2, :, :],
        rmsef[:, 0:jf-2, :, :],
        np.maximum(rf[:, 1:jf-1, :, :] - deltar, 0),
        wprofile
    )
    Vrp = windprofiles(vsef, rmsef, rf + deltar, wprofile)
    Vrm = windprofiles(vsef, rmsef, np.maximum(rf - deltar, 0), wprofile)

    # Convert to meters per second
    V = knotfac * V
    Vpp = knotfac * Vpp
    Vpm = knotfac * Vpm
    Vmp = knotfac * Vmp
    Vmm = knotfac * Vmm
    Vrp = knotfac * Vrp
    Vrm = knotfac * Vrm

    # Calculate wind components
    vph = 0.5 * (V + Vrp)
    vmh = 0.5 * (V + Vrm)
    vnetp = np.sqrt(vph**2 + 2 * vph * latfac * u1temp + u2temp)
    vnetm = np.sqrt(vmh**2 + 2 * vmh * latfac * u1temp + u2temp)
    uekp = -Hi * Cdp * vph * vnetp
    uekm = -Hi * Cdm * vmh * vnetm

    # Calculate Coriolis parameter and gradient wind
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
        (-1 + 2 * ((rf[:, 1:jf-1, :, :] + deltar) /
                   rmsef[:, 1:jf-1, :, :]) ** 2),
        1
    )
    efacm = np.minimum(
        (-1 + 2 * ((rf[:, 1:jf-1, :, :] - deltar) /
                   rmsef[:, 1:jf-1, :, :]) ** 2),
        1
    )

    # Do not consider time rate of change when either velocity is zero
    efacp = efacp * np.minimum(Vpp[:, 1:jf-1, :, :], 1) * np.minimum(
        Vmp[:, 1:jf-1, :, :], 1)
    efacm = efacm * np.minimum(Vpm[:, 1:jf-1, :, :], 1) * np.minimum(
        Vmm[:, 1:jf-1, :, :], 1)
    up[:, 1:jf-1, :, :] = (rf[:, 1:jf-1, :, :] + deltar) * (
        -0.5 * timereswi * efacp * (Vpp[:, 1:jf-1, :, :] - Vmp[:, 1:jf-1, :, :])
        + uekp[:, 1:jf-1, :, :]) / dMdrp[:, 1:jf-1, :, :]
    um[:, 1:jf-1, :, :] = (rf[:, 1:jf-1, :, :] - deltar) * (
        -0.5 * timereswi * efacm * (Vpm[:, 1:jf-1, :, :] - Vmm[:, 1:jf-1, :, :])
        + uekm[:, 1:jf-1, :, :]) / dMdrm[:, 1:jf-1, :, :]

    # Do not calculate vertical velocities if either radial velocity is zero.
    # This is necessary because secondary wind maximum can vanish from one time
    # step to the next, so that time rate of change blows up
    ufac = np.minimum(np.abs(30 * up[:, 1:jf-1, :, :]), 1) * np.minimum(
        np.abs(30 * um[:, 1:jf-1, :, :]), 1)
    w2[:, 1:jf-1, :, :] = (
        -Htrop * ufac * deltari
        * ((rf[:, 1:jf-1, :, :] + deltar) * up[:, 1:jf-1, :, :]
           - (rf[:, 1:jf-1, :, :] - deltar) * um[:, 1:jf-1, :, :])
        / np.maximum(rf[:, 1:jf-1, :, :], 1)
    )

    return w2, cp, V, Vrp, Vrm


def windprofiles(vm, rm, r, wp, vm2=None, rm2=None, opt=True):
    """
    Calculate the radial profiles of azimuthal wind. If secondary wind
    parameters (vm2 and rm2) are provided, the function incorporates
    any secondary eyewalls.

    Parameters:
    -----------
    vm : float
        Maximum circular wind speed (knots)
    rm : float
        Radius of maximum wind (km)
    r : float
        Distance of each event from the point of interest (POI) (km)
    wp : int
        Wind profile shape parameter, used to adjust the wind decay rate
    vm2 : float, optional
        Secondary maximum circular wind speed (knots), if a secondary eyewall is present
        Default is 0, indicating no secondary eyewall
    rm2 : float, optional
        Secondary radius of maximum wind (km)
        Default is 0, indicating no secondary eyewall
    opt : bool, optional
        Option to include additional scaling for the Emanuel wind profile
        Default is True

    Returns:
    --------
    V : array-like
        Azimuthal wind (knots).

    Notes:
    ------
    - The function assumes vm2 and rm2 are set to 0 if no secondary eyewall.
    """

    wprofile = wp  # Use holland (1) or emanuel (2) or er2011 (3) wind profile
    wprofile = np.minimum(wprofile, 3)  # Forbid use of experimental profiles
    vm = vm * 1852 / 3600  # Convert maximum wind speed to m/s

    if vm2 is not None:
        vm2 = vm2 * 1852 / 3600  # Convert maximum wind speed to m/s

    se = 0  # Initialize se to 0
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
        V = vm * (rh ** -bs * np.exp(1 - rh ** -bs)) ** x

        if rm2 is not None and se != 0:
            rh = r / rm2
            x = 0.5 + (rh - 1) * (xn - 0.5) / (rn - 1)
            x[x < 0.5] = 0.5
            V2 = vm2 * (rh ** -bs * np.exp(1 - rh ** -bs)) ** x

    elif wprofile == 2:
        r0 = 1000  # Outer radius (km)

        # Re-scale radius of maximum winds by random number drawn from
        # log-normal distribution
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
                rat ** mb2 * (fac1 / (nb + mb * rat ** fac3) + fac2 / (1 + mb2 * rat ** fac4))
            )
        )
        if rm2 is not None and se != 0:
            rat = r / np.maximum(rm2, 1)
            V2 = (
                vm2
                * (np.maximum((r0 - r), 0) / (r0 - rm2))
                * np.sqrt(
                    rat ** mb2
                    * (fac1 / (nb + mb * rat ** fac3) + fac2 / (1 + mb2 * rat ** fac4))
                )
            )

    elif wprofile == 3:
        crat = 1
        f = 5.0e-5 * 1000  # converts from kms to meters

        if opt:
            Mm = rm * vm + 0.5 * f * rm ** 2
        else:
            Mm = rm * vm
        rn = r / np.where(rm == 0, 1e-32, rm)

        if crat == 1:
            M = Mm * (2 * rn ** 2 / (1 + rn ** 2))
        else:
            M = Mm * (2 * rn ** 2 / (2 - crat + crat * rn ** 2)) ** (1 / (2 - crat))

        # Add long tail to V to avoid discontinuity in vorticity at outer radius
        V = M / (r + 1e-8)
        V[V < 0] = 0
        V = np.nan_to_num(V)

        if rm2 is not None and se != 0:
            Mm = rm2 * vm2
            rn = r / rm2
            if crat == 1:
                M = Mm * (2 * rn ** 2 / (1 + rn ** 2))
            else:
                M = Mm * (2 * rn ** 2 / (2 - crat + crat * rn ** 2)) ** (1 / (2.0 - crat))

            # Add long tail to V to avoid discontinuity in vorticity at outer radius
            V2 = M / (r + 1e-8)
            V2 = np.maximum(V2, 0)
            V2 = np.nan_to_num(V2)

    # Merge primary and secondary wind profiles
    if rm2 is not None and se != 0 and wprofile < 4:
        rm_clipped = np.where(rm == 0, 1e-10, rm)
        u = np.maximum(r, 1) / rm_clipped
        hu = np.maximum(np.sign(u - 1), 0)
        V = V * (1 - hu) + hu * vm / u
        u = r / np.maximum(rm2, 1)
        hu = np.maximum(np.sign(u - 1), 0)
        V2 = V2 * hu + (1 - hu) * vm2 * u
        V = np.maximum(V, V2)

    V = V * 3600 / 1852  # Convert wind speed to knots

    return V


def transfunction(latitude):
    """
    Calculate a translation factor based on latitude.
    This function produces a translation factor (transfactor) based on the
    given latitude, which is added to the circular wind speeds of the synthetic
    track sets.

    Parameters:
    -----------
    latitude : float or array-like
        Latitude(s) for which the translation factor is to be calculated.

    Returns:
    --------
    transfactor : float or array-like
        The calculated translation factor(s).
    """
    transfac = 0.8
    transcap = 1.0
    amplitude = 0.35
    centerlat = 35.0
    latscale = 10.0

    # Calculate the translation factor
    transfactor = transfac + amplitude * (
        1 + np.tanh((np.abs(latitude) - centerlat) / latscale)
    )

    # Cap the translation factor at the transcap value
    transfactor = np.minimum(transfactor, transcap)

    return transfactor


def get_translation_speeds(latitude, longitude):
    """
    Calculate translation speeds (knots) and smooth the results.

    This function reads in track data, calculates the translation speeds in the
    west-east (ut) and north-south (vt) directions, and applies smoothing to
    the results.

    Parameters:
    -----------
    latitude : array-like, shape (n, m)
        Latitude of points along each track.
    longitude : array-like, shape (n, m)
        Longitude of points along each track.

    Returns:
    --------
    ut : array-like, shape (n, m)
        West-east component of the storm translation velocity (knots).
    vt : array-like, shape (n, m)
        North-south component of the storm translation velocity (knots).
    jmax : array-like, shape (n,)
        Time length of each event.
    """
    smfac = 0.4                         # Smoothing factor
    pifac = math.acos(-1) / 180         # Pi number
    dfac = 60 * 1.852                   # 1 nautical mile = 1/60 degree = 1.852 km
    dtfac = 1 / (4 * 1.852)
    netfac = dtfac * dfac
    nn, m = np.shape(latitude)          # nn: number of storms; m: max number of time points
    ut = np.nan * np.zeros((nn, m))
    vt = np.nan * np.zeros((nn, m))
    jmax = np.zeros((nn))               # Time length of each event

    # Loop for each storm
    for n in range(nn):
        lat = latitude[n, :]
        long = longitude[n, :]

        # Get the length of each event by finding the first 0 element
        lat = lat[(lat != 0) & ~np.isnan(lat)]
        jm = len(lat)
        jm = np.maximum(jm, 1)
        jmax[n] = jm

        # Calculate storm translation velocity.
        # The difference of lat/long between 2 consecutive time steps provides
        # the distance that storms move during this time period, thus the
        # translation speed in x & y direction.

        # Handle longitude difference when track crosses the prime meridian
        longdif = long[2:jm] - long[0:jm-2]
        for j in range(jm-2):
            if longdif[j] < -300:
                longdif[j] = longdif[j] + 360
            elif longdif[j] > 300:
                longdif[j] = longdif[j] - 360

        # Include correction for ut at high latitude
        ut[n, 1: jm - 1] = netfac * np.cos(pifac * lat[1: jm - 1]) * longdif
        vt[n, 1: jm - 1] = netfac * (latitude[n, 2:jm] - latitude[n, 0: jm - 2])
        j2 = np.minimum(2, jm)
        j3 = np.maximum(jm - 3, 0)
        ut[n, 0] = 2 * ut[n, 1] - ut[n, j2]
        vt[n, 0] = 2 * vt[n, 1] - vt[n, j2]
        if jm > 0:
            ut[n, jm - 1] = 2.0 * ut[n, jm - 2] - ut[n, j3]
            vt[n, jm - 1] = 2.0 * vt[n, jm - 2] - vt[n, j3]

        # Smooth translation velocity
        vt[n, 1: jm - 1] = vt[n, 1: jm - 1] + smfac * (
            vt[n, 0: jm - 2] + vt[n, 2:jm] - 2 * vt[n, 1: jm - 1]
        )
        ut[n, 1: jm - 1] = ut[n, 1: jm - 1] + smfac * (
            ut[n, 0: jm - 2] + ut[n, 2:jm] - 2 * ut[n, 1: jm - 1]
        )

    ut = np.nan_to_num(ut)
    vt = np.nan_to_num(vt)

    return ut, vt, jmax


def get_translation_speeds_full(lats, longs, vs=None, u850=None, v850=None):
    """
    Calculate translation speed and apply corrections to circular wind speed.

    Parameters:
    -----------
    lats : array-like
        Latitude values along the storm track (degree)
    longs : array-like
        Longitude values along the storm track (degree)
    vs : array-like, optional
        Circular wind speed values along the storm track
    u850 : array-like, optional
        Zonal component of the environmental 850 hPa wind (knots)
    v850 : array-like, optional
        Meridional component of the environmental 850 hPa wind (knots)

    Returns:
    --------
    ut : array-like
        Zonal translation velocity (knots)
    vt : array-like
        Meridional translation velocity (knots)
    uinc : array-like
        Corrected zonal wind speed (knots)
    vinc : array-like
        Corrected meridional wind speed (knots)
    """
    # Calculate translation velocity (knots)
    ut, vt, _ = get_translation_speeds(lats, longs)

    # Add latitude-dependent fraction of translation speed to circular wind
    transfunc = transfunction(lats)
    uinc = transfunc * ut
    vinc = transfunc * vt

    # If env 850 hPa wind available, add baroclinic factor to wind speed
    if u850 is not None and v850 is not None:
        # convert from m/s to knots
        udrift = -0.9 * 3600 / 1852
        vdrift = 1.4 * 3600 / 1852
        vdrift = vdrift * lats / (np.abs(lats) + 1e-8)

        # Apply baroclinic correction
        uinc += 0.5 * (ut - udrift - u850) * (vs / 15)
        vinc += 0.5 * (vt - vdrift - v850) * (vs / 15)

    # Wind speed increments is not allowed to exceed circular gradient wind
    if vs is not None:
        ufrac = vs / (0.1 + np.abs(uinc))
        ufrac = np.minimum(ufrac, 1)
        uinc *= ufrac

        vfrac = vs / (0.1 + np.abs(vinc))
        vfrac = np.minimum(vfrac, 1)
        vinc *= vfrac

    return ut, vt, uinc, vinc


def calculate_upward_velocity_field(
    latitude, longitude, velocity, radius_storm, velocity_secondary,
    radius_storm_secondary, ut, vt, us, vs, plat, plong, h, hx, hy,
    deltar=tcr_params.deltar, timeresw=tcr_params.timeresw, Htrop=tcr_params.Htrop,
    radcity=tcr_params.radcity, wprofile=tcr_params.wprofile
):
    """
    Calculate the spatial distribution of vertical velocity.

    Parameters:
    -----------
    latitude, longitude : array-like
        Latitude and longitude of points along each track (degree)
    velocity : array-like
        Maximum circular wind speed along each track (knots)
    velocity_secondary : array-like
        Maximum circular wind speed of any secondary eyewalls along each track (knots)
    radius_storm : array-like
        Radius of maximum circular wind along each track (km)
    radius_storm_secondary : array-like
        Radius of maximum circular wind of any secondary wind maxima along each track (km)
    ut, vt : array-like
        West-east and north-south components of the storm translation velocity (knots)
    us, vs : array-like
        Vertical shears used to estimate the baroclinic components of the vertical motion (m/s)
    h, hx, hy : array-like
        Topographic heights (m) and their derivatives in x and y.
    deltar : float, optional
        Delta radius (km) for calculating dM/dr, default is 2
    timeresw : float, optional
        Native time resolution (hr) of WRT output, default is 2
    Htrop : float, optional
        Depth of the lower troposphere (m), default is 4000
    radcity : float, optional
        Distance from point of interest (km) beyond which storm influence is
        ignored, default is 300
    wprofile : int, optional
        Wind profile parameter, default is 3

    Returns:
    --------
    w : array-like
        Spatial distribution of vertical velocity (m/s) centered at plong, plat
    """

    se = np.max(velocity_secondary)
    deltari = 1. / deltar
    timereswi = 1. / (3600 * timeresw)
    Hi = 1. / Htrop
    nn = 1
    m = len(ut)
    sx = len(plong)
    sy = len(plat)

    ngrid = 1 if sx == sy else 0
    sy = 1 if ngrid == 1 else np.max(np.shape(plat))

    plong = plong + 360 if plong[0] < 0 else plong
    res = 0.25                          # spatial resolution in degrees
    pifac = tcr_params.RAD2DEG           # convert radians to degrees (pi/180)
    dfac = tcr_params.DEG2KM             # factor converting degree to km
    sfac = 1.0 / (dfac * res * 1000)    # convert resolution from degree to meter
    knotfac = tcr_params.KNOTS2MPS       # convert knots to m/s (1 knots = 0.5144 m/s)
    ut = ut * knotfac   # convert to m/s
    vt = vt * knotfac   # convert to m/s
    omega = tcr_params.OMEGA  # Earth angular velocity parameter (rad/s)
    latfac = (latitude[0] / (abs(latitude[0]) + 1e-8) if latitude.ndim == 1
              else latitude[0, 0] / (abs(latitude[0, 0]) + 1e-8))

    # Estimate Drag coefficients
    cdrag, cdx, cdy = tcr_tb.estimate_drag_coefficients(plat, plong, sfac)
    cdrag = np.maximum(cdrag, 0.0015)

    rfull, dx, dy = tcr_tb.calculate_distance_to_track(
        plat, plong, latitude, longitude, nn, m, sx, sy, ngrid, dfac)
    rfull = np.maximum(rfull, 0.5)
    costhetafull = dx / rfull
    sinthetafull = dy / rfull

    # Initialize arrays
    shape = (nn, 3, sx, sy)
    vfull, rmfull, vsefull, rmsefull = (np.zeros(shape) for _ in range(4))
    latfull, longfull, utfull, vtfull = (np.zeros(shape) for _ in range(4))
    usfull, vsfull, hfull, hyfull, hxfull = (np.zeros(shape) for _ in range(5))
    cdfull, cdyfull, cdxfull = (np.zeros(shape) for _ in range(3))

    # Reshape input arrays to match the output shape
    velocity_reshaped = velocity[:, np.newaxis, np.newaxis]
    radius_storm_reshaped = radius_storm[:, np.newaxis, np.newaxis]
    velocity_secondary_reshaped = velocity_secondary[:, np.newaxis, np.newaxis]
    radius_storm_secondary_reshaped = radius_storm_secondary[:, np.newaxis, np.newaxis]
    latitude_reshaped = latitude[:, np.newaxis, np.newaxis]
    longitude_reshaped = longitude[:, np.newaxis, np.newaxis]
    ut_reshaped = ut[:, np.newaxis, np.newaxis]
    vt_reshaped = vt[:, np.newaxis, np.newaxis]
    us_reshaped = us[:, np.newaxis, np.newaxis]
    vs_reshaped = vs[:, np.newaxis, np.newaxis]

    # Broadcast the reshaped arrays to fill the output arrays
    vfull[0, :, :, :] = np.broadcast_to(velocity_reshaped, (3, sx, sy))
    rmfull[0, :, :, :] = np.broadcast_to(radius_storm_reshaped, (3, sx, sy))
    vsefull[0, :, :, :] = np.broadcast_to(velocity_secondary_reshaped, (3, sx, sy))
    rmsefull[0, :, :, :] = np.broadcast_to(radius_storm_secondary_reshaped, (3, sx, sy))
    latfull[0, :, :, :] = np.broadcast_to(latitude_reshaped, (3, sx, sy))
    longfull[0, :, :, :] = np.broadcast_to(longitude_reshaped, (3, sx, sy))
    utfull[0, :, :, :] = np.broadcast_to(ut_reshaped, (3, sx, sy))
    vtfull[0, :, :, :] = np.broadcast_to(vt_reshaped, (3, sx, sy))
    usfull[0, :, :, :] = np.broadcast_to(us_reshaped, (3, sx, sy))
    vsfull[0, :, :, :] = np.broadcast_to(vs_reshaped, (3, sx, sy))

    w, cp, V, Vd, Vrp, Vrm, u1temp, u2temp, Cdp, Cdm = calculate_wind_primary(
        utfull, vtfull, vfull, rfull, rmfull, vsefull, rmsefull, latfull,
        nn, 3, sx, sy, h, hx, hy, cdrag, hfull, hxfull, hyfull, cdfull, cdx,
        cdxfull, cdy, cdyfull, Hi, Htrop, omega, pifac, deltar, deltari,
        knotfac, latfac, timereswi, costhetafull, sinthetafull, 2, wprofile,
        adj_water=True)

    # If secondary eyewalls present, add in their contribution to wind
    if se > 0:
        w2, cp, V, Vrp, Vrm = calculate_wind_secondary(
            rfull, vsefull, rmsefull, latfull, nn, 3, sx, sy, Hi, Htrop,
            omega, pifac, deltar, deltari, knotfac, latfac, timereswi,
            u1temp, u2temp, Cdp, Cdm, wprofile)
        w = np.maximum(w, w2)

    # Include shear dotted with storm entropy
    hxmod = -0.0005 * (cp + 2 * V / (0.1 + rfull) + deltari * (Vrp - Vrm)) * vsfull
    hymod = 0.0005 * (cp + 2 * V / (0.1 + rfull) + deltari * (Vrp - Vrm)) * usfull
    ufunc = np.clip((radcity - rfull) / 50, 0, 1)
    utfull *= ufunc
    vtfull *= ufunc
    w[:, 1, :, :] = (
        w[:, 1, :, :]
        + (vtfull[:, 1, :, :] + Vd[:, 1, :, :] * latfac *
           costhetafull[:, 1, :, :]) * hyfull[:, 1, :, :]
        + (utfull[:, 1, :, :] - Vd[:, 1, :, :] * latfac *
           sinthetafull[:, 1, :, :]) * hxfull[:, 1, :, :]
        + Vd[:, 1, :, :] * costhetafull[:, 1, :, :] * hymod[:, 1, :, :]
        - Vd[:, 1, :, :] * sinthetafull[:, 1, :, :] * hxmod[:, 1, :, :]
    )
    w = np.minimum(w, 7)
    return w


def calculate_wind_time_series(latitude, longitude, velocity, radius_storm,
                               velocity_secondary, radius_storm_secondary, ut, vt,
                               plat, plong, timeres, timelength=tcr_params.timelength,
                               wheight=tcr_params.wheight, wprofile=tcr_params.wprofile,
                               radcity=tcr_params.radcity):
    """
    Calculate the time series of wind speed and direction at specified points
    of interest.

    Parameters:
    -----------
    latitude : array-like
        Latitudes along each track (degree)
    longitude : array-like
        Longitudes along each track (degree)
    velocity : array-like
        Maximum circular wind speed along each track (knots)
    radius_storm : array-like
        Radius (km) of maximum circular wind along each track (km)
    velocity_secondary : array-like
        Maximum circular wind speed of any secondary eyewalls present (knots)
    radius_storm_secondary : array-like
        Radius (km) of maximum circular wind of any secondary eyewalls (km)
    ut : array-like
        West-east component of the storm translation velocity (knots)
    vt : array-like
        North-south component of the storm translation velocity (knots)
    plat : array-like
        Latitudes of points of interest (degree)
    plong : array-like
        Longitudes of points of interest (degree)
    timeres : float
        Time resolution (hr)
    timelength : int, optional
        Length of the time series at fixed points (default is 96)
    wheight : float, optional
        Altitude above local terrain to estimate surface winds (default is 30 m)
    wprofile : int, optional
        Wind profile parameter (default is 3)
    radcity : int, optional
        Distance from point of interest beyond which storm influence is ignored (km)
        Default is 300

    Returns:
    --------
    vs : array-like
        Wind speed at points of interest.
    direction : array-like
        Wind direction at points of interest.
    """

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
    bathy = tcr_io.load_netcdf_2d_parameters(
        BASE_DATA_DIR, "surface_data.nc", "bathymetry"
    )

    # Load neutral drag coefficients
    cd = tcr_io.load_netcdf_2d_parameters(
        BASE_DATA_DIR, "surface_data.nc", "cdrag"
    )
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
        plat, plong, latitude, longitude, nn, m, sx, sy, 0, dfac)
    radius = np.maximum(radius, 0.5)

    jmin = np.argmin(radius, axis=1)  # where radius of the storm is smallest
    jmin = np.maximum(jmin, 0 + delj)  # cut off delj steps at the begining
    jmin = np.minimum(jmin, m - 1 - delj)  # cut off delj steps at the end
    jstart = np.maximum(jmin - delj, 0)  # index start
    jend = jmin + delj + 1  # index end
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
                if jend[n, i, j] - jstart[n, i, j] < jtot:
                    tind = jend[n, i, j] - jstart[n, i, j]
                else:
                    tind = jtot
                vshort[n, :tind, i, j] = velocity[n, jstart[n, i, j]: jend[n, i, j]]
                rmshort[n, :tind, i, j] = radius_storm[n, jstart[n, i, j]: jend[n, i, j]]
                vseshort[n, :tind, i, j] = velocity_secondary[n, jstart[n, i, j]: jend[n, i, j]]
                rmseshort[n, :tind, i, j] = radius_storm_secondary[n, jstart[n, i, j]: jend[n, i, j]]
                latshort[n, :tind, i, j] = latitude[n, jstart[n, i, j]: jend[n, i, j]]
                longshort[n, :tind, i, j] = longitude[n, jstart[n, i, j]: jend[n, i, j]]
                utshort[n, :tind, i, j] = ut[n, jstart[n, i, j]: jend[n, i, j]]
                vtshort[n, :tind, i, j] = vt[n, jstart[n, i, j]: jend[n, i, j]]

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
            vfine[:, k, :, :] = (1 - weight) * vshort[:, j, :, :] + \
                weight * vshort[:, j + 1, :, :]
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


def calculate_upward_velocity_time_series(
    latitude, longitude, velocity, radius_storm,
    velocity_secondary, radius_storm_secondary, ut, vt, us, vs, plong, plat,
    h, hx, hy, timeres, deltar=tcr_params.deltar, timelength=tcr_params.timelength,
    Htrop=tcr_params.Htrop, wprofile=tcr_params.wprofile, radcity=tcr_params.radcity,
    date_records=None, dq=None, wrad=None,
):
    """
    Calculate the time series of vertical velocity multiplied by saturation specific humidity
    at specified points of interest (POI)

    Parameters:
    -----------
    latitude : array-like
        Latitudes along each track
    longitude : array-like
        Longitudes along each track
    velocity : array-like
        Maximum circular wind along each track (m/s)
    radius_storm : array-like
        Radius of maximum circular wind along each track (km)
    velocity_secondary : array-like
        Maximum circular wind of any secondary eyewalls that may be present (m/s)
    radius_storm_secondary : array-like
        Radius of maximum circular wind of any secondary eyewalls (km)
    ut : array-like
        West-east component of the storm translation velocity (knots)
    vt : array-like
        North-south component of the storm translation velocity (knots)
    us : array-like
        Vertical shear in the west-east direction (m/s)
    vs : array-like
        Vertical shear in the north-south direction (m/s)
    plat : array-like
        Latitudes of POI (degree)
    plong : array-like
        Longitudes of POI (degree)
    h : array-like
        Topographic heights (m)
    hx : array-like
        Gradients of topographic heights in the x direction (-)
    hy : array-like
        Gradients of topographic heights in the y direction (-)
    timeres : float
        Temporal resolution (hour)
    date_records : array-like
        Dates corresponding to each track point
    dq : array-like
        Saturation specific humidity at 900 hPa (g/kg)
    wrad : float
        Background subsidence velocity under radiative cooling
    deltar : float, optional
        Delta radius for calculating dM/dr (km). Default is 2
    timelength : int, optional
        Length of time series at fixed points (hours)
        Default is 96
    Htrop : float, optional
        Depth of lower troposphere (m). Default is 4000
    wprofile : int, optional
        Wind profile (1=Holland, 2=Emanuel, 3=Emanuel & Rotunno 2011)
        Default is 3
    radcity : int, optional
        Distance from point of interest beyond which storm influence is ignored (km)
        Default is 300

    Returns:
    --------
    wq : array-like
        Vertical velocity (m/s).
    date_record : array-like
        Time in date format corresponding to the rain rate.
    """

    timeresi = 1.0 / (3600 * timeres)

    # convert degree to km (1 nautical mile = 1/60 degree = 1.852 km)
    pifac = tcr_params.RAD2DEG
    dfac = 60 * 1.852
    res = 0.25                          # spatial resolution in degrees
    sfac = 1.0 / (dfac * res * 1000)    # convert resolution from degree to meter
    knotfac = tcr_params.KNOTS2MPS       # convert knots to m/s (1 knots = 0.5144 m/s)
    omega = tcr_params.OMEGA             # Earth angular velocity parameter

    se = np.nanmax(velocity_secondary)  # for secondary eyewalls
    ntime = int(timelength / timeres + 1)
    nsteps = round(2.0 / timeres)
    nstepsi = 1.0 / nsteps
    delj = np.floor(timelength / nsteps).astype(int)
    deltari = 1.0 / deltar  # inverse of Delta radius
    Hi = 1.0 / Htrop

    nn, m = ut.shape
    sx = plong.size
    ngrid = 0

    if date_records is not None:
        sy = 1
    else:
        sy = plat.size
        if sx == sy:
            ngrid = 1
            sy = 1

    w = np.zeros((nn, ntime, sx, sy))

    plong = np.where(plong < 0, plong + 360, plong)

    latfac = latitude[0, 0] / (abs(latitude[0, 0]) + 1e-8)
    ut = ut * knotfac  # in m/s
    vt = vt * knotfac  # in m/s

    # Estimate Drag coefficients:
    cdrag, cdx, cdy = tcr_tb.estimate_drag_coefficients(plat, plong, sfac)

    # Reduce drag coefficient of near-coastal locations
    cdrag = np.minimum(cdrag, 1.5e-3 * (1 + np.maximum(h, 0) / 100))
    cdx = cdx * 0.01 * np.minimum(np.maximum(h, 0), 100)
    cdy = cdy * 0.01 * np.minimum(np.maximum(h, 0), 100)

    # Calculate radius distance from storm center to POI
    radius, dx, dy = tcr_tb.calculate_distance_to_track(
        plat, plong, latitude, longitude, nn, m, sx, sy, ngrid, dfac)
    radius = np.maximum(radius, 0.5)

    # Get index where the radius is smallest to start
    jmin = np.argmin(radius, axis=1)  # where the radius of the storm smallest (peak of the storm)
    jmin = np.maximum(jmin, 0 + delj)  # cut off delj steps at the begining
    jmin = np.minimum(jmin, m - 1 - delj)  # cut off delj steps at the end
    jstart = np.maximum(jmin - delj, 0)  # index start
    jend = jmin + delj + 1  # index end
    jtot = 2 * delj + 1  # total index?
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
    if date_records is not None:
        dateshort = np.zeros((nn, jtot, sx, sy))
        dqshort = np.zeros((nn, jtot, sx, sy))

    for i in range(sx):
        for j in range(sy):
            for n in range(nn):
                if jend[n, i, j] - jstart[n, i, j] < jtot:
                    tind = jend[n, i, j] - jstart[n, i, j]
                else:
                    tind = jtot
                vshort[n, :tind, i, j] = velocity[n, jstart[n, i, j]:jend[n, i, j]]
                rmshort[n, :tind, i, j] = radius_storm[n, jstart[n, i, j]:jend[n, i, j]]
                vseshort[n, :tind, i, j] = velocity_secondary[n, jstart[n, i, j]:jend[n, i, j]]
                rmseshort[n, :tind, i, j] = radius_storm_secondary[n, jstart[n, i, j]:jend[n, i, j]]
                rshort[n, :tind, i, j] = radius[n, jstart[n, i, j]:jend[n, i, j], i, j]
                latshort[n, :tind, i, j] = latitude[n, jstart[n, i, j]:jend[n, i, j]]
                longshort[n, :tind, i, j] = longitude[n, jstart[n, i, j]:jend[n, i, j]]
                utshort[n, :tind, i, j] = ut[n, jstart[n, i, j]:jend[n, i, j]]
                vtshort[n, :tind, i, j] = vt[n, jstart[n, i, j]:jend[n, i, j]]
                usshort[n, :tind, i, j] = us[n, jstart[n, i, j]:jend[n, i, j]]
                vsshort[n, :tind, i, j] = vs[n, jstart[n, i, j]:jend[n, i, j]]
                if date_records is not None:
                    dateshort[n, :tind, i, j] = date_records[n, jstart[n, i, j]:jend[n, i, j]]
                    dqshort[n, :tind, i, j] = dq[n, jstart[n, i, j]:jend[n, i, j]]

    # Create high time-resolution series
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

    if date_records is not None:
        date_record = np.zeros((nn, jfine, sx, sy))
        dqfine = np.zeros((nn, jfine, sx, sy))

    k = 0
    for j in range(jtot - 1):
        for n in range(nsteps):
            weight = n * nstepsi
            vfine[:, k, :, :] = (1 - weight) * vshort[:, j, :, :] + \
                weight * vshort[:, j + 1, :, :]
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
            usfine[:, k, :, :] = (1 - weight) * usshort[:, j, :, :] + \
                weight * usshort[:, j + 1, :, :]
            vsfine[:, k, :, :] = (1 - weight) * vsshort[:, j, :, :] + \
                weight * vsshort[:, j + 1, :, :]
            if date_records is not None:
                date_record[:, k, :, :] = (1 - weight) * dateshort[:, j, :, :] + \
                    weight * dateshort[:, j + 1, :, :]
                dqfine[:, k, :, :] = (1 - weight) * dqshort[:, j, :, :] + \
                    weight * dqshort[:, j + 1, :, :]
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
    if date_records is not None:
        date_record[:, k, :, :] = dateshort[:, jtot - 1, :, :]
        dqfine[:, k, :, :] = dqshort[:, jtot - 1, :, :]

    rmsefine = np.maximum(rmsefine, 0.1)
    rfine, dx, dy = tcr_tb.calculate_distance_to_track(
        plat, plong, latfine, longfine, nn, jfine, sx, sy, ngrid, dfac)
    rfinei = 1 / np.maximum(rfine, 1)

    # for primary eyewalls
    w, cp, V, Vd, Vrp, Vrm, u1temp, u2temp, Cdp, Cdm = calculate_wind_primary(
        utfine, vtfine, vfine, rfine, rmfine, vsefine, rmsefine, latfine,
        nn, jfine, sx, sy, h, hx, hy, cdrag, hfine, hxfine, hyfine,
        cdfine, cdx, cdxfine, cdy, cdyfine, Hi, Htrop, omega, pifac,
        deltar, deltari, knotfac, latfac, timeresi, dx * rfinei, dy * rfinei, 10,
        wprofile, adj_water=False)

    # if secondary eyewalls present
    if se > 0:
        w2, cp, V, Vrp, Vrm = calculate_wind_secondary(
            rfine, vsefine, rmsefine, latfine, nn, jfine, sx, sy, Hi, Htrop, omega,
            pifac, deltar, deltari, knotfac, latfac, timeresi, u1temp, u2temp, Cdp,
            Cdm, wprofile)
        w = np.maximum(w, w2)

    # Now add in topographic and shear components to include shear dotted with storm entropy
    hxmod = -0.0005 * (cp + 2 * V / (0.1 + rfine) + deltari * (Vrp - Vrm)) * vsfine
    hymod = 0.0005 * (cp + 2 * V / (0.1 + rfine) + deltari * (Vrp - Vrm)) * usfine

    # Reduce effect of translation speed outside of radcity
    ufunc = np.clip((radcity - rfine) / 50.0, 0, 1)
    utfine = utfine * ufunc
    vtfine = vtfine * ufunc

    # Reduce effect of orography outside of storm core
    ufunc = np.clip((150 - rfine) / 30, 0.2, 0.6)
    hxfine *= ufunc
    hyfine *= ufunc

    w = (
        w
        + (Vd * latfac * dx * rfinei + vtfine) * hyfine
        + (utfine - Vd * latfac * dy * rfinei) * hxfine
        + Vd * dx * rfinei * hymod
        - Vd * dy * rfinei * hxmod
    )
    w = np.minimum(w, 7)

    # Add radiative cooling. If w+wrad is negative (downward motion): wp = 0
    if date_records is not None:
        wq = dqfine * np.maximum(w + wrad, 0)
        return wq, date_record
    else:
        return w


def integrate_outer_wind_profile(vm, fc, ro, wc, Cd, q):
    """
    Numerically integrates the outer wind profile from a simple ordinary
    differential equation (ODE).

    Parameters:
    -----------
    vm : float
        Maximum wind speed (m/s)
    fc : float
        Coriolis parameter (s^-1)
    ro : float
        Outer radius (km)
    wc : float
        Radiative subsidence rate (mm/s)
    Cd : float
        Drag coefficient.
    q : int
        Number of radial points.

    Returns:
    --------
    v : array of float
        Potential intensity at each radial point (m/s)
    r : array of float
        Radius at each radial point (km)
    imin : int
        Minimum index for which v and r are defined
    """

    # Asselin filter coefficient
    assl = 0.2

    # Convert units
    ro *= 1000.0  # Convert to meters
    wc *= 0.001   # Convert to m/s
    chi = Cd * fc * ro / wc  # Definition of chi
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
        r[i] = r[i + 1] - dr
        m[i] = (m[i + 2] - 2.0 * dr * (chi * m[i + 1]**2 /
                (rond**2 - r[i + 1]**2) - r[i + 1]))
        m[i + 1] += assl * (m[i] + m[i + 2] - 2.0 * m[i + 1])
        v[i] = m[i] / r[i]
        if v[i] > 1.0:
            imin = i
            break

    # Fill in values inside radius where v = potential intensity
    for i in range(imin - 1):
        v[i] = 0.0
        r[i] = dr * i

    # Re-dimensionalize
    for i in range(q):
        v[i] = vm * v[i]                # v in m/s
        r[i] = 0.001 * vm * r[i] / fc   # r in km

    epa = (v[imin] - vm) / (v[imin] - v[imin + 1])
    rm = r[imin] + epa * (r[imin + 1] - r[imin])

    return v, rm, imin


def estimate_radius_wind(ds, lat_tracks, vmax_tracks, id_tracks,
                         data_directory=os.path.join(DOWNSCALED_DATA_DIR, 'downscaled'),
                         model=None, basin=None, expmnt=None, force_recompute=False):
    """
    Estimate the radius of maximum circular wind from maximum circular wind speed.

    Parameters:
    -----------
    ds : xarray.Dataset
        Dataset containing the data
    lat_tracks : np.ndarray
        Latitude for the wind tracks (degree)
    vmax_tracks : np.ndarray
        Maximum wind or intensity (m/s)
    id_tracks : np.ndarray
        IDs of the TC tracks
    data_directory : str
        Folder containing the dataset
    model : str, optional
        Model name.
    basin : str, optional
        Basin name.
    expmnt : str, optional
        Experiment name
    force_recompute : bool, optional
        If True, force recomputation of the radius of maximum circular wind

    Returns:
    --------
    np.ndarray
        Radius of maximum circular wind (km).
    """

    ro = 1000.0       # Outer radius
    wc = 3.0          # Radiative subsidence rate
    cdouter = 1.2e-3  # Drag coefficient
    nouter = 1000     # Number of radial point
    num_tcs = len(id_tracks)

    if 'rm_trks' not in ds.keys() or force_recompute:
        print('Estimating the radius of maximum circular wind... ')
        rm_trks = np.zeros(vmax_tracks.shape)
        for ind in range(num_tcs):
            vmax = vmax_tracks[ind, :]
            vmax = vmax[~np.isnan(vmax)]
            jmax = len(vmax)
            for jnd in range(jmax):
                alats = lat_tracks[ind, jnd]
                vsin = vmax_tracks[ind, jnd]
                # Skip at equator (coriolis force fc = 0) or vs = 0
                if alats != 0 and vsin > 0:
                    fc1 = 2 * tcr_params.OMEGA * np.sin(np.radians(abs(alats)))
                    _, rm, _ = integrate_outer_wind_profile(vsin, fc1, ro, wc, cdouter, nouter)
                else:
                    rm = 0
                rm_trks[ind, jnd] = rm
            print(f"Tracks {ind + 1}/{num_tcs} completed...", end='\r')
            sys.stdout.flush()
        print('\nDone!')

        ds = ds.assign(rm_trks=(['n_trk', 'time'], rm_trks))

        # Remove old dataset and save updated dataset to disk
        ncfile = glob.glob(os.path.join(data_directory, expmnt, f"tracks_{basin}_{model}_*.nc"))[0]
        if os.path.exists(ncfile):
            os.remove(ncfile)
        ds.to_netcdf(ncfile, mode='w')
    else:
        print('The radius of maximum circular wind was found in file!')
        print('Loading rm_trks...')
        rm_trks = ds['rm_trks'].values

    return rm_trks


def smooth_2d_array(x, nz, jmin, jmax):
    """
    Apply a 1-2-1 smoothing filter to a 2D array

    Parameters:
    -----------
    x : numpy.ndarray
        Input 2D array to be smoothed
    nz : int
        Number of rows in the 2D array
    jmin : int
        Minimum column index to start smoothing
    jmax : int
        Maximum column index to end smoothing

    Returns:
    --------
    numpy.ndarray
        Smoothed 2D array.
    """
    xsmooth = np.copy(x)  # Make a copy of the input array to avoid modiying it

    # Apply the 1-2-1 smoothing filter
    for i in range(1, nz):
        for j in range(jmin, jmax):
            xsmooth[i, j] = (
                0.125 * (x[i - 1, j] + x[i + 1, j] + x[i, j - 1] + x[i, j + 1])
                + 0.5 * x[i, j]
            )

    return xsmooth


def calculate_wind_swath(
    nt, latitude, longitude, radius_storm, velocity, radius_storm_secondary,
    velocity_secondary, uinc, vinc, extent=None, shapefile=None, 
    magfac=tcr_params.magfac, deltax=tcr_params.deltax, deltay=tcr_params.deltay,
    dellatlongs=tcr_params.dellatlongs, timeres=tcr_params.timeres
):
    """
    Calculate the distribution of maximum point wind speed (in knots) for a
    single storm, including the effects of primary and secondary eyewalls.

    Parameters:
    -----------
    nt : int
        Track number of the storm
    latitude : numpy.ndarray
        Latitudes along each track (degree)
    longitude : numpy.ndarray
        Longitudes along each track (degree)
    radius_storm : numpy.ndarray
        Radius of maximum circular wind along each track (km)
    velocity : numpy.ndarray
        Maximum circular wind along each storm track (knots)
    radius_storm_secondary : numpy.ndarray
        Radius of maximum circular wind of any secondary eyewalls (km)
    velocity_secondary : numpy.ndarray
        Maximum circular wind of any secondary eyewalls that may be present (knots)
    uinc : numpy.ndarray
        West-east component of the storm translation velocity (knots)
    vinc : numpy.ndarray
        Zonal & meridional components of the 850 hPa env wind speed (knots)
    extent : tuple, optional
        Bounding box and spacing in data coordinates (left, right, bottom, top)
        Defines the spatial extent of the map, default is None
    shapefile : str or shapefile-like object, optional
        A shapefile to overlay on the map
    magfac : float, optional
        Overall scale factor for storm size, default is 1.0
    deltax : float, optional
        Longitudinal distance of map boundaries from storm center (degree), default is 5
    deltay : float, optional
        Latitudinal distance of map boundaries from storm center (degree), default is 4
    dellatlongs : float, optional
        Horizontal resolution of swath maps (degree), default is 0.15
    timeres : float, optional
        Time resolution for time series at fixed points (hr), default is 2

    Returns:
    --------
    x : numpy.ndarray
        Longitudes of the grid (degree)
    y : numpy.ndarray
        Latitudes of the grid (degree)
    maxwind : numpy.ndarray
        Storm maximum wind speed at each point on the grid (knots)
    """

    if extent is None and shapefile is None:
        bxmin, bxmax, bymin, bymax = 20, 380, -60, 60
    else:
        if shapefile is None:
            bxmin, bxmax, bymin, bymax = extent
        else:
            bxmin, bxmax, bymin, bymax = tcr_io.get_bbox_from_shapefile(
                shapefile)

    # Apply the overall scale factor to the radius of maximum circular wind
    scaled_radius_storm = radius_storm * magfac
    scaled_radius_storm_secondary = radius_storm_secondary * magfac

    # Find the index q where latitude is closest to 0
    q = np.argmin(np.abs(latitude[nt, :])) - 1

    # Extract non-zero elements and transpose them
    utd = uinc[nt, :][uinc[nt, :] != 0][np.newaxis, :]
    vtd = vinc[nt, :][vinc[nt, :] != 0][np.newaxis, :]
    lat = latitude[nt, :][latitude[nt, :] != 0][np.newaxis, :]
    long = longitude[nt, :][longitude[nt, :] != 0][np.newaxis, :]
    v = velocity[nt, :][velocity[nt, :] != 0][np.newaxis, :]

    qv = np.max(v.shape)
    rm = scaled_radius_storm[nt, :][scaled_radius_storm[nt, :] != 0][
        np.newaxis, :]
    vse = velocity_secondary[nt, :qv][np.newaxis, :]
    rmse = scaled_radius_storm_secondary[nt, :qv][np.newaxis, :]

    # Adjust longitudes crossing the 0/360 boundary
    for i in range(q):
        if long[0, 0] > 200 and long[0, i] < 50 and long[0, i] != 0:
            long[0, i] += 360

    # Calculate map boundaries
    bxmin = np.min(long[np.nonzero(long)]) - deltax
    bxmax = np.max(long[np.nonzero(long)]) + deltax
    bymin = np.min(lat[np.nonzero(lat)]) - deltay
    bymax = np.max(lat[np.nonzero(lat)]) + deltay

    # Create x and y grid
    x = np.arange(bxmin, bxmax + 1e-10, dellatlongs)
    y = np.arange(bymin, bymax + 1e-10, dellatlongs)

    sx, sy = x.size, y.size

    # Adjust x size if needed
    if sx == sy:
        x = np.append(x, bxmax + dellatlongs)
        sx += 1

    # Initialize maxwind array and compute maxwind
    maxwind = np.zeros((sx, sy))

    wind_speed, _ = calculate_wind_time_series(lat, long, v, rm, vse, rmse, utd, vtd, y, x, timeres)
    maxwind[:, :] = np.max(wind_speed[0, :, :, :], axis=0)

    # Smooth the maxwind array
    maxwind = smooth_2d_array(maxwind, sx - 1, 1, sy - 1)
    maxwind = np.transpose(maxwind)

    return x, y, maxwind
