# -*- coding: utf-8 -*-
# Copyright (C) Katelyn Breivik (2017 - 2019)
#
#
# GW_calcs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GW_calcs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#


'''`GW_calcs`
'''

import numpy as np
import pandas as pd
import scipy.special as ss
from scipy import integrate
from scipy.interpolate import interp1d

# Constants
Rsun_au = 1/215.032
day_yr = 1/365.25
yr_sec = 3.155e7
G = 6.67e-11
c = 2.998e+8
Msun = 1.989e30
parsec = 3.086e+19

def chirp_mass(m1, m2):
    """Computes the chirp mass
    
    Parameters
    ----------
    m1 : float/array
        primary mass
        
    m2 : float/array
        secondary mass
        
    Returns
    -------
    m_chirp : float/array
        chirp mass in units supplied
    """
    m_chirp = (m1*m2)**(3./5.)/(m1+m2)**(1./5.)
    return m_chirp

def inspiral_time_peters(a0,e0,m1,m2,af=0):
    """
    Computes the inspiral time, in yr, for a binary
    a0 in Au, and masses in solar masses
    
    if different af is given, computes the time from a0,e0
    to that af
    
    for af=0, just returns inspiral time
    for af!=0, returns (t_insp,af,ef)
    """
    coef = 6.086768e-11 #G^3 / c^5 in au, gigayear, solar mass units
    beta = (64./5.) * coef * m1 * m2 * (m1+m2)
    
    if e0 == 0:
        #print e0,a0
        if not af == 0:
            print("ERROR: doesn't work for circular binaries")
            return 0
        return a0**4 / (4*beta) * 1e9
    
    c0 = a0 * (1.-e0**2.) * e0**(-12./19.) * (1.+(121./304.)*e0**2.)**(-870./2299.)
    
    if af == 0:
        eFinal = 0.
    else:
        r = ode(deda_peters)
        r.set_integrator('lsoda')
        r.set_initial_value(e0,a0)
        r.integrate(af)
        if not r.successful():
            print("ERROR, Integrator failed!")
        else:
            eFinal = r.y[0]      
    
    time_integrand = lambda e: e**(29./19.)*(1.+(121./304.)*e**2.)**(1181./2299.) / (1.-e**2.)**1.5
    integral,abserr = integrate.quad(time_integrand,eFinal,e0)
    
    if af==0:
        return integral * (12./19.) * c0**4. / beta * 1e9
    else:
        return (integral * (12./19.) * c0**4. / beta), af, eFinal * 1e9


def peters_t_merge(m1, m2, a0, ecc):
    """Computes the merger time of a very eccentric binary
    following Peters 63
    
    Params
    ------
    a_0 : float/array
        initial separation in au
    m1 : float/array
        primary mass in solar masses
    m2 : float/array
        secondary mass in solar masses   
    ecc : float/array
        initial eccentricity
        
    Returns
    -------
    t_ecc : float/array
        merger time in years
    """
    G_au_msun_yr = 39.42
    c_au_yr = 6.3198e4
    beta = 64./5.*G_au_msun_yr**3*m1*m2*(m1+m2)/c_au_yr**5
    
    t_circ = a0**4/(4*beta)
    t_ecc = 768/425*t_circ*(1-ecc**2)**(7./2)
    
    return t_ecc

def peters_t_merge_circ(m1, m2, a0):
    """Computes the merger time of a circular binary
    following Peters 63
    
    Params
    ------
    a_0 : float/array
        initial separation in au
    m1 : float/array
        primary mass in solar masses
    m2 : float/array
        secondary mass in solar masses   
    ecc : float/array
        initial eccentricity
        
    Returns
    -------
    t_circ : float/array
        merger time in years
    """
    G_au_msun_yr = 39.42
    c_au_yr = 6.3198e4
    beta = 64./5.*G_au_msun_yr**3*m1*m2*(m1+m2)/c_au_yr**5
    
    t_circ = a0**4/(4*beta)    
    return t_circ
    
def F_e(e):
    """Computes the eccentricity factor in GW power 
    from Peters & Matthews 1964
    
    Params
    ------ 
    e : float/array
        initial eccentricity
        
    Returns
    -------
    F_e : float/array
        ecc factor
    """
    
    nominator = 1 + 73./24.*e**2 + 37./96.*e**4
    denominator = (1 - e**2)**(7./2.)
    F_e = nominator/denominator
    return F_e

def peters_g(e, n):
    """Computes the factor g_n from Peters & Mathews 1963

    Parameters
    ----------
    ecc : float or array
        eccentricity
    n : int
        nth frequency harmonic

    Returns
    -------
    g_fac : array
        array of g_n
    """
    g_fac = (n**4 / 32.0)*((ss.jv((n-2), (n*e)) -\
                            2*e*ss.jv((n-1), (n*e)) +\
                            2.0/n*ss.jv((n), (n*e)) +\
                            2*e*ss.jv((n+1), (n*e)) -\
                            ss.jv((n+2), (n*e)))**2 +\
                           (1-e**2)*(ss.jv((n-2), (n*e)) -\
                                     2*ss.jv((n), (n*e)) +\
                                     ss.jv((n+2), (n*e)))**2 +\
                           (4/(3.0*n**2))*(ss.jv((n), (n*e)))**2)

    return g_fac

def peak_gw_freq(m1, m2, ecc, porb):
    """Computes the peak gravitational-wave frequency for an
    eccentric binary system. Units are SI

    Parameters
    ----------
    m1 : float or array
        primary mass [kg]
    m2 : float or array    
        secondary mass [kg]
    ecc : float or array
        eccentricity
    porb : float or array
        orbital period [s]

    Returns
    -------
    f_gw_peak : float or array
        peak gravitational-wave frequency [Hz]
    """

    # convert the orbital period into a separation using Kepler III
    sep_m = (G/(4*np.pi**2)*porb**2*(m1+m2))**(1./3.)

    f_gw_peak = ((G*(m1+m2))**0.5/np.pi) * (1+ecc)**1.1954/(sep_m*(1-ecc)**2)**1.5
    return f_gw_peak


def de_dt(m1, m2, a_0, e_0):
    """Computes the semimajor axis after t_evol years

    Parameters
    ----------
    m1 : float/array
        primary mass in solar masses
    m2 : float/array
        secondary mass in solar masses   
    a_0 : float/array
        initial separation in au
    e_0 : float or array
        initial eccentricity
    t_evol : float/array
        evolution time in years
    Parameters
    ----------
    de_dt : float or array
        rate of change of eccentricity
    """
    
    # assumes all units are SI
    c_0 = a_0*(1-e_0**2)*e_0**(-12./19)*(1+(121./304)*e_0**2)**(-870./2299)
    G_au_msun_yr = 39.42
    c_au_yr = 6.3198e4
    beta = 64./5.*G_au_msun_yr**3*m1*m2*(m1+m2)/c_au_yr**5
    
    front_fac = -19./12.*beta*c_0**(-4)
    nominator = e_0**(-29./19)*(1-e_0**2)**(3./2)
    denominator = (1+(121./304)*e_0**2)**(1181./2299)
    
    return front_fac*nominator/denominator

def peters_sep_ecc(m1, m2, a_0, e_0, t_evol, evol=False):
    """Computes the semimajor axis after t_evol years

    Parameters
    ----------
    m1 : float/array
        primary mass in solar masses
    m2 : float/array
        secondary mass in solar masses   
    a_0 : float/array
        initial separation in au
    e_0 : float
        initial eccentricity
    t_evol : float/array
        evolution time in years
    evol : bool
        if True: return lists of evolution
        
    Returns
    -------
    sep : float
        final semimajor axis
    ecc : float
        final eccentricity
    """
    t_merge = peters_t_merge_circ(m1, m2, a_0)
    if t_merge <= t_evol:
        return 0.0, 0.0
    else:
        c_0 = a_0*(1-e_0**2)*e_0**(-12./19)*(1+(121./304)*e_0**2)**(-870./2299)
        times = np.linspace(0, t_evol, 100)
        delta = times[1] - times[0]
        if evol:
            e_evol = []
            a_evol = []
            for t in times[1:]:
                e = e_0 + de_dt(m1, m2, a_0, e_0)*delta
                a = c_0*e**(12./19)/(1-e**2)*(1+(121./304)*e**2)**(870./2299)
                e_evol.append(e)
                a_evol.append(a)
                e_0 = e
                a_0 = a
                c_0 = a_0*(1-e_0**2)*e_0**(-12./19)*(1+(121./304)*e_0**2)**(-870./2299)
                if a < 0.001:
                    a_evol.append(0.0)
                    e_evol.append(0.0)
                    break
        else:    
            for t in times[1:]:
                e_evol = e_0 + de_dt(m1, m2, a_0, e_0)*delta
                a_evol = c_0*e_evol**(12./19)/(1-e_evol**2)*(1+(121./304)*e_evol**2)**(870./2299)
                a_0 = a_evol
                e_0 = e_evol
                c_0 = a_0*(1-e_0**2)*e_0**(-12./19)*(1+(121./304)*e_0**2)**(-870./2299)
                if a_evol < 0.001:
                    a_evol=0.0
                    e_evol=0.0
                    break
        return a_evol, e_evol

def peters_a_circ(a_0, m1, m2, t_evol):
    """Computes the evolution of a circular binary due to GWs
    
    Parameters
    ----------
    a_0 : float/array
        initial separation in au
    m1 : float/array
        primary mass in solar masses
    m2 : float/array
        secondary mass in solar masses   
    t_evol : float/array
        evolution time in years
        
    Returns
    -------
    a_final : float/array
        final separation in au
    """
    G_au_msun_yr = 39.42
    c_au_yr = 6.3198e4
    
    beta = 64./5.*G_au_msun_yr**3*m1*m2*(m1+m2)/c_au_yr**5
    a_final = (a_0**4 - 4*beta*t_evol)**(1./4.)
    
    return a_final
def f_dot_n(m1, m2, f_orb, e, n):
    """Computes the frequency evolution of the nth harmonic
    of the orbital frequency
    
    Parameters
    ----------
    m1 : float/array
        primary mass in solar masses
    m2 : float/array
        secondary mass in solar masses   
    f_orb : float/array
        orbital frequency [Hz]
    ecc : float/array
        eccentricity
    n : float array
        harmonic number (n=2 for cirular binaries)
        
    Returns
    -------
    f_dot_n : float/array
        frequency evolution in s^(-2)
    """
    m_c = chirp_mass(m1, m2)*Msun
    front_fac = 96./(10*np.pi)*G**(5./3.)/c**5
    nominator = (m_c)**(5./3.)*(2*np.pi*f_orb)**(11./3.)*n*F_e(e)
    f_dot_n = front_fac*nominator
    
    return f_dot_n

def h_circ_stationary(m1, m2, porb, d):
    """Computes the dimneionsless strain of a stationary source
    This is also written as h_n in the COSMIC paper
    
    Parameters
    ----------
    m1 : float/array
        primary mass in solar masses
    m2 : float/array
        secondary mass in solar masses   
    porb : float/array
        orbital period [sec]
    d : float/array
        luminosity distance [kpc]
        
    Returns
    -------
    h_c_circ : float/array
        characteristic strain of stationary sources
    """
    m_c = chirp_mass(m1, m2)
    porb = porb/3600
    h_c_circ = 0.5e-21*(m_c)**(5./3.)*(porb)**(-2./3)/(d)
    
    return h_c_circ

def h_c_circ_chirping(m1, m2, sep, d, Tobs):
    """Computes the characteristic strain of a chirping source
    This is also written as h_c,2
    
    Parameters
    ----------
    m1 : float
        primary mass in solar masses
    m2 : float
        secondary mass in solar masses   
    f_orb : float
        orbital frequency [Hz]
    sep : float
        orbital separation [rsun]
    d : float
        luminosity distance [kpc]
    Tobs : float
         Lisa observation time [yrs]
        
    Returns
    -------
    h_c_circ : float/array
        characteristic strain of stationary sources
    """
    a_final = peters_a_circ(a_0=sep*Rsun_au, m1=m1, m2=m2, t_evol=Tobs)
    p_final = porb_from_sep(m1=m1, m2=m2, sep=a_final)
    f_final = 2/(p_final*yr_sec)

    p_initial = porb_from_sep(m1=m1, m2=m2, sep=sep*Rsun_au)
    
    f_initial = 2/(p_initial*yr_sec)
    
    f_evol = np.linspace(f_initial, f_final, 5)
    m_c = chirp_mass(m1, m2)*Msun
    d = d*parsec*1000
    front_fac = 2/(3*np.pi**(4./3.))*(G**(5./3))/c**3
    h_c_circ_2 = front_fac*(m_c)**(5./3.)*(2*f_evol)**(-1./3)/(d**2)
    h_c_circ = h_c_circ_2
    
    return h_c_circ, f_evol

def h_ecc_stationary(m1, m2, porb, d, ecc, n):
    """Computes the dimensionelss strain of a stationary source
    This is also written as h_n
    
    Parameters
    ----------
    m1 : float/array
        primary mass in solar masses
    m2 : float/array
        secondary mass in solar masses   
    porb : float/array
        orbital period [sec]
    d : float/array
        luminosity distance [kpc]
    ecc : float/array
        eccentricity
    n : int/array
        nth orbital harmonic
    Tobs : float
        Lisa observation time [yrs]
    
    Returns
    -------
    h_c_ecc : float/array
        characteristic strain of stationary sources
    """
    
    m_c = chirp_mass(m1, m2)
    porb = porb/3600
    h_ecc = 1e-21*(m_c)**(5./3.)*(porb)**(-2./3)/(d*n)*(peters_g(n=n, e=ecc))**0.5
    
    return h_ecc

def h_c_ecc_chirping(m1, m2, f_orb, sep, ecc, d, Tobs):
    """Computes the characteristic strain of a chirping source
    This is also written as h_c,2
    
    Parameters
    ----------
    m1 : float
        primary mass in solar masses
    m2 : float
        secondary mass in solar masses   
    f_orb : float
        orbital frequency [Hz]
    sep : float
        orbital separation [rsun]
    ecc : float
        orbital eccentricity
    d : float
        luminosity distance [kpc]
        
    Returns
    -------
    h_c_circ : float/array
        characteristic strain of stationary sources
    """
    
    t_evol = np.linspace(1, Tobs, 10)
    a_evol = []
    e_evol = []
    a_init = sep*Rsun_au
    e_init = ecc
    a_evol, e_evol = peters_sep_ecc(a_0=a_init, e_0=e_init, m1=m1, m2=m2, t_evol=Tobs, evol=True)
    a_evol = np.array(a_evol)
    e_evol = np.array(e_evol)
    a_evol = a_evol/Rsun_au
    
    p_evol = porb_from_sep(m1=m1, m2=m2, sep=np.array(a_evol)*Rsun_au)
    p_evol = p_evol*day_yr
    f_evol = 2/(p_evol*86400)
    m_c = chirp_mass(m1, m2)*Msun
    d = d*parsec*1000
    h_c_ecc = []
    f_peak = []
    n_array = np.arange(1,100,1)
    for f, e in zip(f_evol, e_evol):
        front_fac = 2/(3*np.pi**(4./3.))*(G**(5./3))/c**3*peters_g(n=n_array, e=e)/F_e(e=e)
        h_c_ecc_power = front_fac*(m_c)**(5./3.)*(n_array*f)**(-1./3)/(d**2)*(2/n_array)**(2./3)
        
        f_peak_gw = peak_gw_freq(m1=m1*Msun, m2=m2*Msun, ecc=e, porb=1/f)
        h_c_ecc.append(np.sum(h_c_ecc_power)**0.5)
        f_peak.append(f_peak_gw)
    
    return h_c_ecc, f_evol

def lisa_characteristic_noise_power():
    '''Computes LISA characteristic strain sensitivity curve according to `Cornish and Robson 2018 <https://arxiv.org/pdf/1803.01944.pdf>`_ without the Galactic foreground
    Parameters
    ----------
    none
    Returns
    -------
    LISA_hc : interpolation of LISA sensitivity curve
    '''
    freq = np.logspace(-9,1,10000)
    # note: freq [Hz], L_arm [m], S_n [Hz^-0.5]
    L_arm = 2.5e9
    f_star = 19.09*1e-3

    P_oms = (1.5e-11)**2*(1. + (2.0e-3/freq)**4)
    P_acc = (3.0e-15)**2*(1. + (0.4e-3/freq)**2)*(1. + (freq/(8.0e-3))**4)

    P_n = (P_oms + 2.*(1. + np.cos(freq/f_star)**2)*P_acc/(2.*np.pi*freq)**4)/L_arm**2
    R = 3./10./(1. + 6./10.*(freq/f_star)**2)
    S_n = (P_n/R*freq)

    LISA_hc = interp1d(freq, S_n)
    return LISA_hc

def lisa_PSD():
    '''Computes LISA characteristic strain sensitivity curve according to `Cornish and Robson 2018 <https://arxiv.org/pdf/1803.01944.pdf>`_ without the Galactic foreground
    Parameters
    ----------
    none
    Returns
    -------
    LISA_hc : interpolation of LISA sensitivity curve
    '''
    freq = np.logspace(-9,1,10000)
    # note: freq [Hz], L_arm [m], S_n [Hz^-0.5]
    L_arm = 2.5e9
    f_star = 19.09*1e-3

    P_oms = (1.5e-11)**2*(1. + (2.0e-3/freq)**4)
    P_acc = (3.0e-15)**2*(1. + (0.4e-3/freq)**2)*(1. + (freq/(8.0e-3))**4)

    P_n = (P_oms + 2.*(1. + np.cos(freq/f_star)**2)*P_acc/(2.*np.pi*freq)**4)/L_arm**2
    R = 3./10./(1. + 6./10.*(freq/f_star)**2)
    S_n = (P_n/R)

    LISA_hc = interp1d(freq, S_n)
    return LISA_hc

def LISA_power(dat, Tobs):
    LISA_power_tot = []
    LISA_freq = np.arange(1e-7, 1e-1, 1/(Tobs*yr_sec))
    dat_circ = dat.loc[dat.ecc < 0.01]
    dat_circ['strain_2'] = h_circ_stationary(m1=np.array(dat_circ.mass_1), 
                                             m2=np.array(dat_circ.mass_2), 
                                             porb=np.array((dat_circ.porb_final*86400)), 
                                             d=np.array(dat_circ.dist))**2
    dat_ecc = dat.loc[dat.ecc >= 0.01]
    dat_circ['digits'] = np.digitize(dat_circ.f_gw_peak, LISA_freq)
    power = dat_circ.groupby('digits').sum()['strain_2']
    power_tot = np.zeros(len(LISA_freq))
    power_tot[np.array(power.index.astype(int))] = power
    power_dat = pd.DataFrame(np.vstack([LISA_freq, power_tot]).T, 
                             columns=['f_gw', 'strain_2'])
    if len(LISA_power_tot) == 0:
        LISA_power_tot = power_dat
    else:
        LISA_power_tot['strain_2'] = LISA_power_tot['strain_2'] + power_dat['strain_2']
    
    if len(dat_ecc) > 0:
        power_tot = []
        freq_tot = []
        for n in range(1, 100):
            h_ecc_n = h_ecc_stationary(m1=np.array(dat_ecc.mass_1), 
                                       m2=np.array(dat_ecc.mass_2), 
                                       porb=np.array((dat_ecc.porb_final*86400)), 
                                       d=np.array(dat_ecc.dist), 
                                       ecc=np.array(dat_ecc.ecc_final),
                                       n=n)**2
            f_gw = np.array(n/(dat_ecc.porb_final*86400))
                            
            power_tot.extend(h_ecc_n)
            freq_tot.extend(f_gw)
    
        power_ecc = pd.DataFrame(np.vstack([freq_tot, power_tot]).T,
                                 columns=['f_gw', 'strain_2'])
        power_ecc['digits'] = np.digitize(power_ecc.f_gw, LISA_freq)
        power_ecc_dat = power_ecc.groupby('digits').sum()['strain_2']
        power_tot_ecc = np.zeros(len(LISA_freq))
        power_tot_ecc[np.array(power_ecc_dat.index.astype(int))-1] = power_ecc_dat
        power_ecc_df = pd.DataFrame(np.array([LISA_freq, power_tot_ecc]).T, 
                             columns=['f_gw', 'strain_2'])
        LISA_power_tot['strain_2'] = LISA_power_tot['strain_2'] + power_ecc_df['strain_2']
    return LISA_power_tot