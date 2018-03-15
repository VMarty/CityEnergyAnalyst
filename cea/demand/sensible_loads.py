# -*- coding: utf-8 -*-
"""
Sensible space heating and space cooling loads
EN-13970
"""
from __future__ import division
import numpy as np
from cea.utilities.physics import BOLTZMANN
from cea.demand import control_heating_cooling_systems

__author__ = "Jimeno A. Fonseca"
__copyright__ = "Copyright 2016, Architecture and Building Systems - ETH Zurich"
__credits__ = ["Jimeno A. Fonseca", "Shanshan Hsieh", "Daren Thomas", "Martin Mosteiro"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Daren Thomas"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


# capacity of emission/control system

def calc_Qhs_Qcs_sys_max(Af, prop_HVAC):
    # TODO: Documentation
    # Refactored from CalcThermalLoads

    IC_max = -prop_HVAC['Qcsmax_Wm2'] * Af
    IH_max = prop_HVAC['Qhsmax_Wm2'] * Af
    return IC_max, IH_max


# solar and heat gains

def calc_Qgain_sen(t, tsd, bpr, gv):
    # TODO

    # internal loads
    tsd['I_sol_and_I_rad'][t], tsd['I_rad'][t], tsd['I_sol'][t] = calc_I_sol(t, bpr, tsd, gv)

    return tsd


def calc_Qgain_lat(schedules, bpr):
    # TODO: Documentation
    # Refactored from CalcThermalLoads
    """
    :param list_uses: The list of uses used in the project
    :type list_uses: list
    :param schedules: The list of schedules defined for the project - in the same order as `list_uses`
    :type schedules: list[ndarray[float]]
    :param X_ghp: humidity gain from people in g/h/p for each occupancy type
    :type X_ghp: list[float]
    :param occupancy: for each use in `list_uses`, the percentage of that use for this building. Sum of values is 1.0
    :type occupancy: dict[str, float]
    :param Af: total conditioned floor area
    :type Af: float

    :param sys_e_heating: cooling system code as defined in the systems database (e.g. 'T0' if no cooling)
    :param sys_e_heating: string
    :param sys_e_cooling: cooling system code as defined in the systems database (e.g. 'T0' if no cooling)
    :param sys_e_cooling: string

    :return w_int: yearly schedule

    """
    # calc yearly humidity gains based on occupancy schedule and specific humidity gains for each occupancy type in the
    # building
    humidity_schedule = schedules['X'] * bpr.internal_loads['X_ghp']  # in g/h/m2
    if control_heating_cooling_systems.heating_system_is_ac(
            bpr) or control_heating_cooling_systems.cooling_system_is_ac(bpr):
        w_int = humidity_schedule * bpr.rc_model['Af'] / (1000 * 3600)  # kg/s
    else:
        w_int = np.zeros(8760)
        # FIXME: should humidity gains also be considered for cooling = 'T2' ? I think so. Changed.

    return w_int


def calc_I_sol(t, bpr, tsd, gv):
    """
    This function calculates the net solar radiation (incident -reflected - re-irradiated) according to ISO 13790

    :param t: hour of the year
    :param bpr: building properties object
    :param tsd: time series dataframe
    :param gv: global variables class
    :return:
        I_sol_net: vector of net solar radiation to the building
        I_rad: vector solar radiation re-irradiated to the sky.
        I_sol_gross : vector of incident radiation to the building.
    """

    # calc irradiation to the sky
    I_rad = calc_I_rad(t, tsd, bpr, gv.Rse)

    # get incident radiation
    I_sol_gross = bpr.solar.I_sol[t]

    I_sol_net = I_sol_gross + I_rad

    return I_sol_net, I_rad, I_sol_gross  # vector in W


def calc_I_rad(t, tsd, bpr, Rse):
    """
    This function calculates the solar radiation re-irradiated from a building to the sky according to ISO 13790

    :param t: hour of the year
    :param tsd: time series dataframe
    :param bpr:  building properties object
    :param gv: global variables class
    :return:
        I_rad: vector solar radiation re-irradiated to the sky.
    """

    temp_s_prev = tsd['theta_c'][t - 1]
    if np.isnan(tsd['theta_c'][t - 1]):
        temp_s_prev = tsd['T_ext'][t - 1]

    theta_ss = tsd['T_sky'][t] - temp_s_prev
    Fform_wall, Fform_win, Fform_roof = 0.5, 0.5, 1  # 50% reiradiated by vertical surfaces and 100% by horizontal.
    I_rad_win = Rse * bpr.rc_model['U_win'] * calc_hr(bpr.architecture.e_win, theta_ss) * bpr.rc_model[
        'Aw'] * theta_ss
    I_rad_roof = Rse * bpr.rc_model['U_roof'] * calc_hr(bpr.architecture.e_roof, theta_ss) * bpr.rc_model[
        'Aroof'] * theta_ss
    I_rad_wall = Rse * bpr.rc_model['U_wall'] * calc_hr(bpr.architecture.e_wall, theta_ss) * bpr.rc_model[
        'Aop_sup'] * theta_ss
    I_rad = Fform_wall * I_rad_wall + Fform_win * I_rad_win + Fform_roof * I_rad_roof

    return I_rad


def calc_hr(emissivity, theta_ss):
    """
    This function calculates the external radiative heat transfer coefficient according to ISO 13790

    :param emissivity: emissivity of the considered surface
    :param theta_ss: delta of temperature between building surface and the sky.
    :return:
        hr:

    """
    return 4.0 * emissivity * BOLTZMANN * (theta_ss + 273.0) ** 3.0


# temperature of emission/control system

def calc_temperatures_emission_systems(tsd, bpr, Qcsf_0, Qhsf_0, gv):
    from cea.technologies import radiators, heating_coils, tabs
<<<<<<< HEAD
    # local variables
    Ta_0 = np.nanmax(tsd['ta_hs_set'])
    if bpr.hvac['type_hs'] == 'T0':
        Ths_sup = np.zeros(8760)  # in C
        Ths_re = np.zeros(8760)  # in C
        mcphs = np.zeros(8760)  # in KW/C

    if bpr.hvac['type_cs'] == 'T0':
        Tcs_re = np.zeros(8760)  # in C
        Tcs_sup = np.zeros(8760)  # in C
        mcpcs = np.zeros(8760)  # in KW/C

    if bpr.hvac['type_hs'] == 'T1' or bpr.hvac['type_hs'] == 'T2':  # radiators

        Ths_sup, Ths_re, mcphs = np.vectorize(radiators.calc_radiator)(tsd['Qhsf'], tsd['T_int'], Qhsf_0, Ta_0,
                                                                       bpr.building_systems['Ths_sup_0'],
                                                                       bpr.building_systems['Ths_re_0'])

    if bpr.hvac['type_hs'] == 'T3':  # air conditioning
        index = np.where(tsd['Qhsf'] == Qhsf_0)
        ma_sup_0 = tsd['ma_sup_hs'][index[0][0]]
        Ta_sup_0 = tsd['Ta_sup_hs'][index[0][0]] + 273
        Ta_re_0 = tsd['Ta_re_hs'][index[0][0]] + 273
        Ths_sup, Ths_re, mcphs = np.vectorize(heating_coils.calc_heating_coil)(tsd['Qhsf'], Qhsf_0, tsd['Ta_sup_hs'],
                                                                               tsd['Ta_re_hs'],
                                                                               bpr.building_systems['Ths_sup_0'],
                                                                               bpr.building_systems['Ths_re_0'],
                                                                               tsd['ma_sup_hs'], ma_sup_0,
                                                                               Ta_sup_0, Ta_re_0, gv.Cpa, gv)

    if bpr.hvac['type_cs'] == 'T2':  # mini-split units

        index = np.where(tsd['Qcsf'] == Qcsf_0)
        ma_sup_0 = tsd['ma_sup_cs'][index[0][0]]
        Ta_sup_0 = tsd['Ta_sup_cs'][index[0][0]] + 273
        Ta_re_0 = tsd['Ta_re_cs'][index[0][0]] + 273
        Tcs_sup, Tcs_re, mcpcs = np.vectorize(heating_coils.calc_cooling_coil)(tsd['Qcsf'], Qcsf_0, tsd['Ta_sup_cs'],
                                                                               tsd['Ta_re_cs'],
                                                                               bpr.building_systems['Tcs_sup_0'],
                                                                               bpr.building_systems['Tcs_re_0'],
                                                                               tsd['ma_sup_cs'], ma_sup_0,
                                                                               Ta_sup_0, Ta_re_0, gv.Cpa, gv)

    if bpr.hvac['type_cs'] == 'T3':  # air conditioning

        index = np.where(tsd['Qcsf'] == Qcsf_0)
        ma_sup_0 = tsd['ma_sup_cs'][index[0][0]]
        Ta_sup_0 = tsd['Ta_sup_cs'][index[0][0]] + 273
        Ta_re_0 = tsd['Ta_re_cs'][index[0][0]] + 273
        Tcs_sup, Tcs_re, mcpcs = np.vectorize(heating_coils.calc_cooling_coil)(tsd['Qcsf'], Qcsf_0, tsd['Ta_sup_cs'],
                                                                               tsd['Ta_re_cs'],
                                                                               bpr.building_systems['Tcs_sup_0'],
                                                                               bpr.building_systems['Tcs_re_0'],
                                                                               tsd['ma_sup_cs'], ma_sup_0,
                                                                               Ta_sup_0, Ta_re_0, gv.Cpa, gv)

    if bpr.hvac['type_hs'] == 'T4':  # floor heating

        Ths_sup, Ths_re, mcphs = np.vectorize(tabs.calc_floorheating)(tsd['Qhsf'], tsd['theta_m'], Qhsf_0,
                                                                      bpr.building_systems['Ths_sup_0'],
                                                                      bpr.building_systems['Ths_re_0'],
                                                                      bpr.rc_model['Af'])
=======

    #
    # TEMPERATURES HEATING SYSTEMS
    #
    if not control_heating_cooling_systems.has_heating_system(bpr):
        # if no heating system

        tsd['Thsf_sup_ahu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['Thsf_re_ahu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['mcphsf_ahu'] = np.zeros(8760)
        tsd['Thsf_sup_aru'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['Thsf_re_aru'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['mcphsf_aru'] = np.zeros(8760)
        tsd['Thsf_sup_shu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['Thsf_re_shu'] = np.zeros(8760) * np.nan # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['mcphsf_shu'] = np.zeros(8760)

    elif control_heating_cooling_systems.has_radiator_heating_system(bpr):
        # if radiator heating system
        Ta_heating_0 = np.nanmax(tsd['ta_hs_set'])
        Qhsf_0 = np.nanmax(tsd['Qhsf'])  # in W

        tsd['Thsf_sup_ahu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['Thsf_re_ahu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['mcphsf_ahu'] = np.zeros(8760)
        tsd['Thsf_sup_aru'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['Thsf_re_aru'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['mcphsf_aru'] = np.zeros(8760)

        Ths_sup, Ths_re, mcphs = np.vectorize(radiators.calc_radiator)(tsd['Qhsf'], tsd['T_int'], Qhsf_0, Ta_heating_0,
                                                                       bpr.building_systems['Ths_sup_shu_0'],
                                                                       bpr.building_systems['Ths_re_shu_0'])

        tsd['Thsf_sup_shu'] = Ths_sup
        tsd['Thsf_re_shu'] = Ths_re
        tsd['mcphsf_shu'] = mcphs

    elif control_heating_cooling_systems.has_central_ac_heating_system(bpr):

        # ahu
        # consider losses according to loads of systems
        frac_ahu = [ahu / sys if sys > 0 else 0 for ahu, sys in zip(tsd['Qhs_sen_ahu'], tsd['Qhs_sen_sys'])]
        qhsf_ahu = tsd['Qhs_sen_ahu'] + (tsd['Qhs_em_ls'] + tsd['Qhs_dis_ls']) * frac_ahu

        Qhsf_ahu_0 = np.nanmax(qhsf_ahu)  # in W

        index = np.where(qhsf_ahu == Qhsf_ahu_0)
        ma_sup_0 = tsd['ma_sup_hs_ahu'][index[0][0]]
        Ta_sup_0 = tsd['ta_sup_hs_ahu'][index[0][0]] + 273
        Ta_re_0 = tsd['ta_re_hs_ahu'][index[0][0]] + 273
        Ths_sup, Ths_re, mcphs = np.vectorize(heating_coils.calc_heating_coil)(qhsf_ahu, Qhsf_ahu_0, tsd['ta_sup_hs_ahu'],
                                                                               tsd['ta_re_hs_ahu'],
                                                                               bpr.building_systems['Ths_sup_ahu_0'],
                                                                               bpr.building_systems['Ths_re_ahu_0'],
                                                                               tsd['ma_sup_hs_ahu'], ma_sup_0,
                                                                               Ta_sup_0, Ta_re_0, C_P_A)
        tsd['Thsf_sup_ahu'] = Ths_sup  # in C
        tsd['Thsf_re_ahu'] = Ths_re  # in C
        tsd['mcphsf_ahu'] = mcphs

        # ARU
        # consider losses according to loads of systems
        frac_aru = [aru / sys if sys > 0 else 0 for aru, sys in zip(tsd['Qhs_sen_aru'], tsd['Qhs_sen_sys'])]
        qhsf_aru = tsd['Qhs_sen_aru'] + (tsd['Qhs_em_ls'] + tsd['Qhs_dis_ls']) * frac_aru

        Qhsf_aru_0 = np.nanmax(qhsf_aru)  # in W

        index = np.where(qhsf_aru == Qhsf_aru_0)
        ma_sup_0 = tsd['ma_sup_hs_aru'][index[0][0]]
        Ta_sup_0 = tsd['ta_sup_hs_aru'][index[0][0]] + 273
        Ta_re_0 = tsd['ta_re_hs_aru'][index[0][0]] + 273
        Ths_sup, Ths_re, mcphs = np.vectorize(heating_coils.calc_heating_coil)(qhsf_aru, Qhsf_aru_0,
                                                                               tsd['ta_sup_hs_aru'],
                                                                               tsd['ta_re_hs_aru'],
                                                                               bpr.building_systems['Ths_sup_aru_0'],
                                                                               bpr.building_systems['Ths_re_aru_0'],
                                                                               tsd['ma_sup_hs_aru'], ma_sup_0,
                                                                               Ta_sup_0, Ta_re_0, C_P_A)
        tsd['Thsf_sup_aru'] = Ths_sup  # in C
        tsd['Thsf_re_aru'] = Ths_re  # in C
        tsd['mcphsf_aru'] = mcphs

        # SHU
        tsd['Thsf_sup_shu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['Thsf_re_shu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['mcphsf_shu'] = np.zeros(8760)

    elif control_heating_cooling_systems.has_floor_heating_system(bpr):

        Qhsf_0 = np.nanmax(tsd['Qhsf'])  # in W

        tsd['Thsf_sup_ahu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['Thsf_re_ahu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['mcphsf_ahu'] = np.zeros(8760)
        tsd['Thsf_sup_aru'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['Thsf_re_aru'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['mcphsf_aru'] = np.zeros(8760)

        Ths_sup, Ths_re, mcphs = np.vectorize(tabs.calc_floorheating)(tsd['Qhsf'], tsd['theta_m'], Qhsf_0,
                                                                      bpr.building_systems['Ths_sup_shu_0'],
                                                                      bpr.building_systems['Ths_re_shu_0'],
                                                                      bpr.rc_model['Af'])
        tsd['Thsf_sup_shu'] = Ths_sup
        tsd['Thsf_re_shu'] = Ths_re
        tsd['mcphsf_shu'] = mcphs

    else:
        raise Exception('Heating system not defined in function: "calc_temperatures_emission_systems"')


    #
    # TEMPERATURES COOLING SYSTEMS
    #
    if not control_heating_cooling_systems.has_cooling_system(bpr):
        # if no heating system

        tsd['Tcsf_sup_ahu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['Tcsf_re_ahu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['mcpcsf_ahu'] = np.zeros(8760)
        tsd['Tcsf_sup_aru'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['Tcsf_re_aru'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['mcpcsf_aru'] = np.zeros(8760)
        tsd['Tcsf_sup_scu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['Tcsf_re_scu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['mcpcsf_scu'] = np.zeros(8760)

    elif control_heating_cooling_systems.has_central_ac_cooling_system(bpr):

        # AHU
        # consider losses according to loads of systems
        frac_ahu = [ahu / sys if sys < 0 else 0 for ahu, sys in zip(tsd['Qcs_sen_ahu'], tsd['Qcs_sen_sys'])]
        qcsf_ahu = tsd['Qcs_sen_ahu'] + tsd['Qcs_lat_ahu'] + (tsd['Qcs_em_ls'] + tsd['Qcs_dis_ls']) * frac_ahu

        Qcsf_ahu_0 = np.nanmin(qcsf_ahu)  # in W

        index = np.where(qcsf_ahu == Qcsf_ahu_0)
        ma_sup_0 = tsd['ma_sup_cs_ahu'][index[0][0]]
        Ta_sup_0 = tsd['ta_sup_cs_ahu'][index[0][0]] + 273
        Ta_re_0 = tsd['ta_re_cs_ahu'][index[0][0]] + 273
        Tcs_sup, Tcs_re, mcpcs = np.vectorize(heating_coils.calc_cooling_coil)(qcsf_ahu, Qcsf_ahu_0, tsd['ta_sup_cs_ahu'],
                                                                               tsd['ta_re_cs_ahu'],
                                                                               bpr.building_systems['Tcs_sup_ahu_0'],
                                                                               bpr.building_systems['Tcs_re_ahu_0'],
                                                                               tsd['ma_sup_cs_ahu'], ma_sup_0,
                                                                               Ta_sup_0, Ta_re_0, C_P_A)
        tsd['Tcsf_sup_ahu'] = Tcs_sup  # in C
        tsd['Tcsf_re_ahu'] = Tcs_re  # in C
        tsd['mcpcsf_ahu'] = mcpcs

        # ARU
        # consider losses according to loads of systems
        frac_aru = [aru / sys if sys < 0 else 0 for aru, sys in zip(tsd['Qcs_sen_aru'], tsd['Qcs_sen_sys'])]
        qcsf_aru = tsd['Qcs_sen_aru'] + tsd['Qcs_lat_aru'] + (tsd['Qcs_em_ls'] + tsd['Qcs_dis_ls']) * frac_aru

        Qcsf_aru_0 = np.nanmin(qcsf_aru)  # in W

        index = np.where(qcsf_aru == Qcsf_aru_0)
        ma_sup_0 = tsd['ma_sup_cs_aru'][index[0][0]]
        Ta_sup_0 = tsd['ta_sup_cs_aru'][index[0][0]] + 273
        Ta_re_0 = tsd['ta_re_cs_aru'][index[0][0]] + 273
        Tcs_sup, Tcs_re, mcpcs = np.vectorize(heating_coils.calc_cooling_coil)(qcsf_aru, Qcsf_aru_0, tsd['ta_sup_cs_aru'],
                                                                               tsd['ta_re_cs_aru'],
                                                                               bpr.building_systems['Tcs_sup_aru_0'],
                                                                               bpr.building_systems['Tcs_re_aru_0'],
                                                                               tsd['ma_sup_cs_aru'], ma_sup_0,
                                                                               Ta_sup_0, Ta_re_0, C_P_A)
        tsd['Tcsf_sup_aru'] = Tcs_sup  # in C
        tsd['Tcsf_re_aru'] = Tcs_re  # in C
        tsd['mcpcsf_aru'] = mcpcs

        # SCU
        tsd['Tcsf_sup_scu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['Tcsf_re_scu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['mcpcsf_scu'] = np.zeros(8760)

    elif control_heating_cooling_systems.has_local_ac_cooling_system(bpr):

        # AHU
        tsd['Tcsf_sup_ahu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['Tcsf_re_ahu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['mcpcsf_ahu'] = np.zeros(8760)

        # ARU
        # consider losses according to loads of systems
        qcsf_aru = tsd['Qcs_sen_aru'] + tsd['Qcs_lat_aru'] + (tsd['Qcs_em_ls'] + tsd['Qcs_dis_ls'])
        qcsf_aru = np.nan_to_num(qcsf_aru)

        # Calc nominal temperatures of systems
        Qcsf_aru_0 = np.nanmin(qcsf_aru)  # in W

        index = np.where(qcsf_aru == Qcsf_aru_0)
        ma_sup_0 = tsd['ma_sup_cs_aru'][index[0][0]]
        Ta_sup_0 = tsd['ta_sup_cs_aru'][index[0][0]] + 273
        Ta_re_0 = tsd['ta_re_cs_aru'][index[0][0]] + 273
        Tcs_sup, Tcs_re, mcpcs = np.vectorize(heating_coils.calc_cooling_coil)(qcsf_aru, Qcsf_aru_0,
                                                                               tsd['ta_sup_cs_aru'],
                                                                               tsd['ta_re_cs_aru'],
                                                                               bpr.building_systems['Tcs_sup_aru_0'],
                                                                               bpr.building_systems['Tcs_re_aru_0'],
                                                                               tsd['ma_sup_cs_aru'], ma_sup_0,
                                                                               Ta_sup_0, Ta_re_0, C_P_A)
        tsd['Tcsf_sup_aru'] = Tcs_sup  # in C
        tsd['Tcsf_re_aru'] = Tcs_re  # in C
        tsd['mcpcsf_aru'] = mcpcs

        # SCU
        tsd['Tcsf_sup_scu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['Tcsf_re_scu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['mcpcsf_scu'] = np.zeros(8760)

    elif control_heating_cooling_systems.has_3for2_cooling_system(bpr):

        # AHU
        # consider losses according to loads of systems
        frac_ahu = [ahu/sys if sys < 0 else 0 for ahu, sys in zip(tsd['Qcs_sen_ahu'],tsd['Qcs_sen_sys'])]
        qcsf_ahu = tsd['Qcs_sen_ahu'] + tsd['Qcs_lat_ahu'] + (tsd['Qcs_em_ls'] + tsd['Qcs_dis_ls']) * frac_ahu
        qcsf_ahu = np.nan_to_num(qcsf_ahu)

        Qcsf_ahu_0 = np.nanmin(qcsf_ahu)  # in W

        index = np.where(qcsf_ahu == Qcsf_ahu_0)
        ma_sup_0 = tsd['ma_sup_cs_ahu'][index[0][0]]
        Ta_sup_0 = tsd['ta_sup_cs_ahu'][index[0][0]] + 273
        Ta_re_0 = tsd['ta_re_cs_ahu'][index[0][0]] + 273
        Tcs_sup, Tcs_re, mcpcs = np.vectorize(heating_coils.calc_cooling_coil)(qcsf_ahu, Qcsf_ahu_0,
                                                                               tsd['ta_sup_cs_ahu'],
                                                                               tsd['ta_re_cs_ahu'],
                                                                               bpr.building_systems['Tcs_sup_ahu_0'],
                                                                               bpr.building_systems['Tcs_re_ahu_0'],
                                                                               tsd['ma_sup_cs_ahu'], ma_sup_0,
                                                                               Ta_sup_0, Ta_re_0, C_P_A)
        tsd['Tcsf_sup_ahu'] = Tcs_sup  # in C
        tsd['Tcsf_re_ahu'] = Tcs_re  # in C
        tsd['mcpcsf_ahu'] = mcpcs

        # ARU
        # consider losses according to loads of systems
        frac_aru = [aru / sys if sys < 0 else 0 for aru, sys in zip(tsd['Qcs_sen_aru'], tsd['Qcs_sen_sys'])]
        qcsf_aru = tsd['Qcs_sen_aru'] + tsd['Qcs_lat_aru'] + (tsd['Qcs_em_ls'] + tsd['Qcs_dis_ls']) * frac_aru
        qcsf_aru = np.nan_to_num(qcsf_aru)

        # Calc nominal temperatures of systems
        Qcsf_aru_0 = np.nanmin(qcsf_aru)  # in W

        index = np.where(qcsf_aru == Qcsf_aru_0)
        ma_sup_0 = tsd['ma_sup_cs_aru'][index[0][0]]
        Ta_sup_0 = tsd['ta_sup_cs_aru'][index[0][0]] + 273
        Ta_re_0 = tsd['ta_re_cs_aru'][index[0][0]] + 273
        Tcs_sup, Tcs_re, mcpcs = np.vectorize(heating_coils.calc_cooling_coil)(qcsf_aru, Qcsf_aru_0,
                                                                               tsd['ta_sup_cs_aru'],
                                                                               tsd['ta_re_cs_aru'],
                                                                               bpr.building_systems['Tcs_sup_aru_0'],
                                                                               bpr.building_systems['Tcs_re_aru_0'],
                                                                               tsd['ma_sup_cs_aru'], ma_sup_0,
                                                                               Ta_sup_0, Ta_re_0, C_P_A)
        tsd['Tcsf_sup_aru'] = Tcs_sup  # in C
        tsd['Tcsf_re_aru'] = Tcs_re  # in C
        tsd['mcpcsf_aru'] = mcpcs

        # SCU
        # consider losses according to loads of systems
        frac_scu = [scu / sys if sys < 0 else 0 for scu, sys in zip(tsd['Qcs_sen_scu'], tsd['Qcs_sen_sys'])]
        qcsf_scu = tsd['Qcs_sen_scu'] + (tsd['Qcs_em_ls'] + tsd['Qcs_dis_ls']) * frac_scu
        qcsf_scu = np.nan_to_num(qcsf_scu)

        Qcsf_scu_0 = np.nanmin(qcsf_scu)  # in W
        Ta_cooling_0 = np.nanmin(tsd['ta_cs_set'])

        Tcs_sup, Tcs_re, mcpcs = np.vectorize(radiators.calc_radiator)(qcsf_scu, tsd['T_int'], Qcsf_scu_0, Ta_cooling_0,
                                                                       bpr.building_systems['Tcs_sup_scu_0'],
                                                                       bpr.building_systems['Tcs_re_scu_0'])
        tsd['Tcsf_sup_scu'] = Tcs_sup  # in C
        tsd['Tcsf_re_scu'] = Tcs_re  # in C
        tsd['mcpcsf_scu'] = mcpcs

    elif control_heating_cooling_systems.has_ceiling_cooling_system(bpr):

        # SCU
        # consider losses according to loads of systems
        qcsf_scu = tsd['Qcs_sen_scu'] + (tsd['Qcs_em_ls'] + tsd['Qcs_dis_ls'])
        qcsf_scu = np.nan_to_num(qcsf_scu)

        Qcsf_scu_0 = np.nanmin(qcsf_scu)  # in W
        Ta_cooling_0 = np.nanmin(tsd['ta_cs_set'])

        # use radiator for ceiling cooling calculation
        Tcs_sup, Tcs_re, mcpcs = np.vectorize(radiators.calc_radiator)(qcsf_scu, tsd['T_int'], Qcsf_scu_0, Ta_cooling_0,
                                                                       bpr.building_systems['Tcs_sup_scu_0'],
                                                                       bpr.building_systems['Tcs_re_scu_0'])

        tsd['Tcsf_sup_scu'] = Tcs_sup  # in C
        tsd['Tcsf_re_scu'] = Tcs_re  # in C
        tsd['mcpcsf_scu'] = mcpcs

        # AHU
        tsd['Tcsf_sup_ahu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['Tcsf_re_ahu'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['mcpcsf_ahu'] = np.zeros(8760)

        # ARU
        tsd['Tcsf_sup_aru'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['Tcsf_re_aru'] = np.zeros(8760) * np.nan  # in C  #FIXME: I don't like that non-existing temperatures are 0
        tsd['mcpcsf_aru'] = np.zeros(8760)
>>>>>>> parent of d11b194c... Merge pull request #1120 from architecture-building-systems/i1070-update-building-energy-balance-dashboard

    return Tcs_re, Tcs_sup, Ths_re, Ths_sup, mcpcs, mcphs  # C,C, C,C, W/C, W/C


# space heating/cooling losses

def calc_Qhs_Qcs_dis_ls(tair, text, Qhs, Qcs, tsh, trh, tsc, trc, Qhs_max, Qcs_max, D, Y, SystemH, SystemC, Bf, Lv):
    """calculates distribution losses based on ISO 15316"""
    # Calculate tamb in basement according to EN
    tamb = tair - Bf * (tair - text)
    if SystemH != 'T0' and Qhs > 0:
        Qhs_d_ls = ((tsh + trh) / 2 - tamb) * (Qhs / Qhs_max) * (Lv * Y)
    else:
        Qhs_d_ls = 0
    if SystemC != 'T0' and Qcs < 0:
        Qcs_d_ls = ((tsc + trc) / 2 - tamb) * (Qcs / Qcs_max) * (Lv * Y)
    else:
        Qcs_d_ls = 0

    return Qhs_d_ls, Qcs_d_ls
