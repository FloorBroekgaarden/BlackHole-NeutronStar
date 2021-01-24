#!/usr/bin/env python3
import numpy as np
import scipy.special
from  scipy.optimize    import newton
from  astropy.cosmology import FlatLambdaCDM
from   astropy.cosmology import WMAP9 #as cosmo
import astropy.units as u


import pandas as pd
import ClassGWcalculations

# Constants
Rsun_au = 1/215.032
day_yr = 1/365.25
yr_sec = 3.155e7
G = 6.67e-11
c = 2.998e+8
Msun = 1.989e30
parsec = 3.086e+19







# # class MSSFR(object):
#     """
#     This class is to calculate the metallicity specific star formation
#     rate in a specific metallicity bin.

#     It combines a
#     -  MZ -relation     : galaxy stellar mass - metallicity relatiion
#                           This translates a metallicity to a 
#                           a galaxy stellar mass which has the same 
#                           average metallicity.
#     -  SFR-prescription : Star formation rate prescription.
#                           Amount of solar mass that goes into forming
#                           stars.
#     -  GSMF function    : Galaxy stellar mass function
#                           A density function of the distribution of galaxy stellar masses.


#     The entire GSMF is the mass in all stars in all galaxies. Hence the fraction of the GSMF
#     is the fraction of the stellar mass in all galaxies at that time. 
#     We assume that the average cosmic SFR spreads evenly amongs all the galaxy stellar mass.
#     Meaning a galaxy with a stellar mass twice that of another galaxy, has twice
#     the SFR. Hence this method does not include local star burts (SMC/LMC for example)
    

#     The overall idea. Given a metallicity bin with upper and lower bound.
#     Translate this to a galaxy stellar mass bin with upper and lower bound using
#     a MZ-relation. Then find the fraction of the GSMF density function that is taken by 
#     this glaxy stellar mass bin. We assumed that the fraction of the stellar mass is 
#     proportional to the fraction of the SFR. Hence we multiply this fraction by the SFR.
#     Note that because we take "cosmological averages" this method does not take into 
#     account galaxy specific metallicity distributions.


#     Alternatively we can also directly use a metallicity density/probability function
#     and calculate the fraction occipied by the metallicity bin in this distribution

#     """

#     # def __init__(self, verbose=False, MilkyWayType="",\
#     #              ):

#     #     #if True we print statements
#     #     self.verbose               = verbose



#     #     if self.metallicityGrid is not None:
#     #         self.calculateMetallicityBinEdges()
#     #     else:
#     #         print("warning MSSFR instance without metallicityGrid")
#     #     if self.totalMassEvolvedPerZ is None:
#     #         print("warning no normalisation of mass evolved assuming 1")
#     #         try:
#     #             self.totalMassEvolvedPerZ = np.ones(len(metallicityGrid))
#     #         except:
#     #             print("cant do it because we have no Metallicity grid")


        
# # -*- coding: utf-8 -*-
# # Copyright (C) Katelyn Breivik (2017 - 2019)
# #
# #
# # MW Maker is free software: you can redistribute it and/or modify
# # it under the terms of the GNU General Public License as published by
# # the Free Software Foundation, either version 3 of the License, or
# # (at your option) any later version.
# #
# # MW Maker is distributed in the hope that it will be useful,
# # but WITHOUT ANY WARRANTY; without even the implied warranty of
# # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# # GNU General Public License for more details.
# #


# '''`MW_Maker`
# '''




def porb_from_sep(m1, m2, sep):
    """Gives orbital period from m1, m2, sep
    following Kepler III
    
    Parameters
    ----------
    m1 : float/array
        primary mass in msun
        
    m2 : float/array
        secondary mass in msun
        
    sep : float/array
        binary separation in au
        

    Returns
    -------
    porb : float/array
        orbital period in years
    """
    
    porb_2 = sep**3/(m1+m2)
    return porb_2**0.5
    
def SFH(n_pop, model):
    if model == 'ThinDisk':
        # times in Myr following COSMIC
        t_birth = np.random.uniform(0, 10000, n_pop)
    elif model == 'Bulge':
        # times in Myr following COSMIC
        t_birth = np.random.uniform(0, 1000, n_pop)
    elif model == 'ThickDisk':
        # times in Myr following COSMIC
        t_birth = np.random.uniform(0, 1000, n_pop)
    return t_birth

def t_evol_GW(birth_time, formation_time, model):
    if model == 'ThinDisk':
        # times in Myr following COSMIC
        t_evol_GW = 10000 - (birth_time + formation_time)
    elif model == 'Bulge':
        # times in Myr following COSMIC
        t_evol_GW = 10000 - (birth_time + formation_time)
    elif model == 'ThickDisk':
        # times in Myr following COSMIC
        t_evol_GW = 11000 - (birth_time + formation_time)
    return t_evol_GW

def R_RL(q, a):
    """ Computes the Eggleton 1983 Roche Lobe radius
    in units of supplied semimajor axis
    
    Parameters
    ----------
    q : float or array
        m_RL/m_companion
    
    a : float or array
        semimajor axis
        
    Returns
    -------
    RL : float or array
        Roche radius
    """
    numerator = 0.49*q**(2./3.)
    denominator = (0.6*q**(2./3.)+ np.log(1+q**(1./3.)))
    RL = a*numerator/denominator
    return RL


def read_dat(dat_path, kstar_types):
    if len(kstar_types) == 2:
        kstar_1 = kstar_types[0]
        kstar_2 = kstar_types[1]
        kstar_plot = '{0}_{1}'.format(kstar_1, kstar_2)
        kstar_1 = [kstar_1]
        kstar_2 = [kstar_2]
    else:
        kstar_1 = kstar_types[0]
        kstar_2_lo = kstar_types[1]
        kstar_2_hi = kstar_types[2]
        kstar_plot = '{0}_{1}_{2}'.format(kstar_1, kstar_2_lo, kstar_2_hi)
        kstar_1 = [kstar_1]
        kstar_2 = np.arange(kstar_2_lo, kstar_2_hi)

    
    # NOTE: this is the dataframe so to get the total mass we need to take the final mass
    m_sim = fdata['doubleCompactObjects']['M1'][...].squeeze()
    m_sim_tot = m_sim.max()[0]
    conv = pd.read_hdf(dat_path+'dat_DeltaBurst_short_'+kstar_plot+'.h5', key='conv')
    return conv, m_sim, m_sim_tot

def GW_evol(pop):
    """Computes the final state of binaries in pop
    based on GW evolution over t_evol_GW years
    
    Parameters
    ----------
    pop : pandas.DataFrame
        contains all the binary info for a binary pop that
        has been scaled to the Galaxy
        
    Returns
    -------
    sep_final : array
        final separation after GW evolution
    ecc_final : array
        final eccentricity after GW evolution
    
    """
    sep_final = np.zeros(len(pop))
    ecc_final = np.zeros(len(pop))

    # units: msun, au, years
    ecc = np.array(pop.ecc)
    m1 = np.array(pop.mass_1)
    m2 = np.array(pop.mass_2)
    sep = np.array(pop.sep*Rsun_au)
    times = np.array(pop.t_evol_GW*1e6)
    circ_ind, = np.where(ecc <= 0.01)
    ecc_ind, = np.where(ecc > 0.01)
    
    
    for ind in ecc_ind:
        sep_final[ind], ecc_final[ind] = GW_calcs.peters_sep_ecc(m1=m1[ind], m2=m2[ind], 
                                                                 a_0=sep[ind], 
                                                                 e_0=ecc[ind],
                                                                 t_evol=times[ind])
    t_merge_circ = GW_calcs.peters_t_merge_circ(m1=m1[circ_ind], 
                                       m2=m2[circ_ind], 
                                       a0=sep[circ_ind])
    ind_merge, = np.where(t_merge_circ < times[circ_ind])
    ind_alive, = np.where(t_merge_circ > times[circ_ind])
    sep_final[circ_ind[ind_merge]] = 0.0
    sep_final[circ_ind[ind_alive]] = GW_calcs.peters_a_circ(a_0=sep[circ_ind[ind_alive]], 
                                                            m1=m1[circ_ind[ind_alive]], 
                                                            m2=m2[circ_ind[ind_alive]], 
                                                            t_evol=times[circ_ind[ind_alive]])
    ecc_final[circ_ind] = 0.0

    #convert the separation back to Rsun
    sep_final = sep_final/Rsun_au
    return sep_final, ecc_final

def gx_sample(conv, model, n_pop):
    # sample the Gx population and assign birth and GW evol times
    pop = conv.sample(n_pop, replace=True)
    pop['t_birth'] = SFH(model=model, n_pop=n_pop)
    pop['t_evol_GW'] = t_evol_GW(birth_time=pop.t_birth, 
                                 formation_time=pop.tphys, 
                                 model=model)
    
    # select out systems which haven't evolved yet
    pop = pop.loc[pop.t_evol_GW >= 0]
    
    # evolve the systems according to Peters 64 evolution
    sep_final, ecc_final = GW_evol(pop)
    porb_final = porb_from_sep(m1=pop.mass_1, m2=pop.mass_2, sep=sep_final*Rsun_au)/day_yr
    pop['sep_final'] = sep_final
    pop['ecc_final'] = ecc_final
    pop['porb_final'] = porb_final

    # select out those which have merged
    pop = pop.loc[(pop.sep_final > 0) & (pop.porb_final*86400 < 1e7)]
    # select out those which fill roche lobes
    if pop.kstar_1.isin([13,14]).all():
        pop = pop
    else:
        pop['RL_1_final'] = R_RL(pop.mass_1/pop.mass_2, pop.sep_final)
        pop = pop.loc[pop.rad_1 <= pop.RL_1_final]
        
    if pop.kstar_2.isin([13,14]).all():
        pop = pop
    else:
        pop['RL_2_final'] = R_RL(pop.mass_2/pop.mass_1, pop.sep_final)
        pop = pop.loc[pop.rad_2 <= pop.RL_2_final]
    
    pop['f_gw_peak'] = GW_calcs.peak_gw_freq(m1=pop.mass_1*2e30, 
                                                     m2=pop.mass_2*2e30, 
                                                     ecc=pop.ecc_final, 
                                                     porb=pop.porb_final*86400)
    return pop

def LISA_Galaxy(Data, model, m_sim_tot, kstars):
    n_pop = MC_samp.mass_weighted_number(component_mass=MC_samp.select_component_mass(model),
                                         dat=conv,
                                         total_sampled_mass=m_sim_tot)
    print('The number of {0} binaries in the {1} is: {2}'.format(kstars, model, n_pop))
    if np.all(conv.kstar_1 < 14):
        conv = conv.loc[conv.sep < 100]
        n_pop = MC_samp.mass_weighted_number(component_mass=MC_samp.select_component_mass(model),
                                             dat=conv,
                                             total_sampled_mass=m_sim_tot)
        print('filtered numbers {}'.format(n_pop))

    if n_pop >= 5e5:
        pop_LISA = []
        for ii in range(0,10):
            # do this in chunks for memory
            n_samp = int(n_pop/10)
            if len(pop_LISA) == 0:
                gx_samp = gx_sample(conv, model, n_samp)
                gx_samp = gx_samp.loc[gx_samp.f_gw_peak > 1e-5]
                pop_LISA = gx_sample(conv, model, n_samp)
            else:
                gx_samp = gx_sample(conv, model, n_samp)
                gx_samp = gx_samp.loc[gx_samp.f_gw_peak > 1e-5]
                pop_LISA = pop_LISA.append(gx_samp)
            print(ii)
    else:
        pop_LISA = gx_sample(conv, model, n_pop)
        if pop_LISA.ecc.all() == 0.0:
            LISA_fgw = pop_LISA.porb_final
    xGx, yGx, zGx, inc, omega, OMEGA = MC_samp.galactic_positions(gx_component=model, model='McMillan', size=len(pop_LISA))
    pop_LISA['xGx'] = xGx
    pop_LISA['yGx'] = yGx
    pop_LISA['zGx'] = zGx
    pop_LISA['dist'] = ((xGx - 8.0)**2 + (yGx - 0)**2 + (zGx - 0.2)**2)**0.5
    
    return pop_LISA







