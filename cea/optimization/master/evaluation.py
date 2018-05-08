"""
====================================
Evaluation function of an individual
====================================
"""
from __future__ import division

import os
import pandas as pd
import numpy as np
import cea.optimization.master.generation as generation
import cea.optimization.master.summarize_network as nM
from cea.optimization.constants import *
import cea.optimization.master.cost_model as eM
import cea.optimization.slave.cooling_main as coolMain
import cea.optimization.slave.slave_main as sM
import cea.optimization.supportFn as sFn
import cea.technologies.substation as sMain
import check as cCheck
from cea.optimization import slave_data


# +++++++++++++++++++++++++++++++++++++
# Main objective function evaluation
# ++++++++++++++++++++++++++++++++++++++

def evaluation_main(individual, building_names, locator, extraCosts, extraCO2, extraPrim, solar_features,
                    network_features, gv, config, prices, ind_num, gen):
    """
    This function evaluates an individual
    :param individual: list with values of the individual
    :param building_names: list with names of buildings
    :param locator: locator class
    :param extraCosts: costs calculated before optimization of specific energy services
     (process heat and electricity)
    :param extraCO2: green house gas emissions calculated before optimization of specific energy services
     (process heat and electricity)
    :param extraPrim: primary energy calculated before optimization ofr specific energy services
     (process heat and electricity)
    :param solar_features: solar features call to class
    :param network_features: network features call to class
    :param gv: global variables class
    :param optimization_constants: class containing constants used in optimization
    :param config: configuration file
    :param prices: class of prices used in optimization
    :type individual: list
    :type building_names: list
    :type locator: string
    :type extraCosts: float
    :type extraCO2: float
    :type extraPrim: float
    :type solar_features: class
    :type network_features: class
    :type gv: class
    :type optimization_constants: class
    :type config: class
    :type prices: class
    :return: Resulting values of the objective function. costs, CO2, prim
    :rtype: tuple
    """
    # Check the consistency of the individual or create a new one
    individual = check_invalid(individual, len(building_names), config)


    # Initialize objective functions costs, CO2 and primary energy
    costs = extraCosts
    CO2 = extraCO2
    prim = extraPrim
    QUncoveredDesign = 0
    QUncoveredAnnual = 0

    # Create the string representation of the individual
    DHN_barcode, DCN_barcode, DHN_configuration, DCN_configuration = sFn.individual_to_barcode(individual, building_names)

    if DHN_barcode.count("1") == gv.num_tot_buildings:
        network_file_name_heating = "Network_summary_result_all.csv"
        Q_DHNf_W = pd.read_csv(locator.get_optimization_network_all_results_summary('all'), usecols=["Q_DHNf_W"]).values
        Q_heating_max_W = Q_DHNf_W.max()
    elif DHN_barcode.count("1") == 0:
        network_file_name_heating = "Network_summary_result_none.csv"
        Q_heating_max_W = 0
    else:
        network_file_name_heating = "Network_summary_result_" + hex(int(str(DHN_barcode), 2)) + ".csv"
        Q_DHNf_W = pd.read_csv(locator.get_optimization_network_results_summary(DHN_barcode), usecols=["Q_DHNf_W"]).values
        Q_heating_max_W = Q_DHNf_W.max()

    if DCN_barcode.count("1") == gv.num_tot_buildings:
        network_file_name_cooling = "Network_summary_result_all.csv"
        if individual[N_HEAT * 2] == 1: # if heat recovery is ON, then only need to satisfy cooling load of space cooling and refrigeration
            Q_DCNf_W = pd.read_csv(locator.get_optimization_network_all_results_summary('all'), usecols=["Q_DCNf_space_cooling_and_refrigeration_W"]).values
        else:
            Q_DCNf_W = pd.read_csv(locator.get_optimization_network_all_results_summary('all'), usecols=["Q_DCNf_space_cooling_data_center_and_refrigeration_W"]).values
        Q_cooling_max_W = Q_DCNf_W.max()
    elif DCN_barcode.count("1") == 0:
        network_file_name_cooling = "Network_summary_result_none.csv"
        Q_cooling_max_W = 0
    else:
        network_file_name_cooling = "Network_summary_result_" + hex(int(str(DCN_barcode), 2)) + ".csv"

        if individual[N_HEAT * 2] == 1: # if heat recovery is ON, then only need to satisfy cooling load of space cooling and refrigeration
            Q_DCNf_W = pd.read_csv(locator.get_optimization_network_results_summary(DCN_barcode), usecols=["Q_DCNf_space_cooling_and_refrigeration_W"]).values
        else:
            Q_DCNf_W = pd.read_csv(locator.get_optimization_network_results_summary(DCN_barcode), usecols=["Q_DCNf_space_cooling_data_center_and_refrigeration_W"]).values
        Q_cooling_max_W = Q_DCNf_W.max()


    Q_heating_nom_W = Q_heating_max_W * (1 + Q_MARGIN_FOR_NETWORK)
    Q_cooling_nom_W = Q_cooling_max_W * (1 + Q_MARGIN_FOR_NETWORK)

    # Modify the individual with the extra GHP constraint
    try:
        cCheck.GHPCheck(individual, locator, Q_heating_nom_W, gv)
    except:
        print "No GHP constraint check possible \n"

    # Export to context
    master_to_slave_vars = calc_master_to_slave_variables(individual, Q_heating_max_W, Q_cooling_max_W, building_names, ind_num, gen)
    master_to_slave_vars.network_data_file_heating = network_file_name_heating
    master_to_slave_vars.network_data_file_cooling = network_file_name_cooling
    master_to_slave_vars.total_buildings = len(building_names)

    if master_to_slave_vars.number_of_buildings_connected_heating > 1:
        if DHN_barcode.count("0") == gv.num_tot_buildings:
            master_to_slave_vars.fNameTotalCSV = locator.get_total_demand()
        else:
            master_to_slave_vars.fNameTotalCSV = os.path.join(locator.get_optimization_network_totals_folder(),
                                                              "Total_%(DHN_barcode)s.csv" % locals())
    else:
        master_to_slave_vars.fNameTotalCSV = locator.get_optimization_substations_total_file(DHN_barcode)

    if master_to_slave_vars.number_of_buildings_connected_cooling > 1:
        if DCN_barcode.count("0") == gv.num_tot_buildings:
            master_to_slave_vars.fNameTotalCSV = locator.get_total_demand()
        else:
            master_to_slave_vars.fNameTotalCSV = os.path.join(locator.get_optimization_network_totals_folder(),
                                                              "Total_%(DCN_barcode)s.csv" % locals())
    else:
        master_to_slave_vars.fNameTotalCSV = locator.get_optimization_substations_total_file(DCN_barcode)

    if config.optimization.isheating:

        if DHN_barcode.count("1") > 0:

            (slavePrim, slaveCO2, slaveCosts, QUncoveredDesign, QUncoveredAnnual) = sM.slave_main(locator,
                                                                                                  master_to_slave_vars,
                                                                                                  solar_features, gv, config, prices)
        else:

            slaveCO2 = 0
            slaveCosts = 0
            slavePrim = 0
    else:
        slaveCO2 = 0
        slaveCosts = 0
        slavePrim = 0

    costs += slaveCosts
    CO2 += slaveCO2
    prim += slavePrim



    print "Add extra costs"
    (addCosts, addCO2, addPrim) = eM.addCosts(DHN_barcode, DCN_barcode, building_names, locator, master_to_slave_vars, QUncoveredDesign,
                                              QUncoveredAnnual, solar_features, network_features, gv, config, prices)

    if gv.ZernezFlag == 1:
        coolCosts, coolCO2, coolPrim = 0, 0, 0
    elif config.optimization.iscooling and DCN_barcode.count("1") > 0:
        (coolCosts, coolCO2, coolPrim) = coolMain.coolingMain(locator, master_to_slave_vars, network_features, gv, prices, config)
    else:
        coolCosts, coolCO2, coolPrim = 0, 0, 0



    costs += addCosts + coolCosts
    CO2 += addCO2 + coolCO2
    prim += addPrim + coolPrim
    # Converting costs into float64 to avoid longer values
    costs = np.float64(costs)
    CO2 = np.float64(CO2)
    prim = np.float64(prim)

    print ('Additional costs = ' + str(addCosts))
    print ('Additional CO2 = ' + str(addCO2))
    print ('Additional prim = ' + str(addPrim))

    print ('Total costs = ' + str(costs))
    print ('Total CO2 = ' + str(CO2))
    print ('Total prim = ' + str(prim))

    return costs, CO2, prim, master_to_slave_vars, individual

#+++++++++++++++++++++++++++++++++++
# Boundary conditions
#+++++++++++++++++++++++++++++


def check_invalid(individual, nBuildings, config):
    """
    This function rejects individuals out of the bounds of the problem
    It can also generate a new individual, to replace the rejected individual
    :param individual: individual sent for checking
    :param nBuildings: number of buildings
    :param gv: global variables class
    :type individual: list
    :type nBuildings: int
    :type gv: class
    :return: new individual if necessary
    :rtype: list
    """
    valid = True

    for i in range(N_HEAT):
        if individual[2 * i] > 0 and individual[2 * i + 1] < 0.01:
            oldValue = individual[2 * i + 1]
            shareGain = oldValue - 0.01
            individual[2 * i + 1] = 0.01

            for rank in range(N_HEAT):
                if individual[2 * rank] > 0 and i != rank:
                    individual[2 * rank + 1] += individual[2 * rank + 1] / (1 - oldValue) * shareGain

        elif individual[2*i] == 0:
            individual[2*i + 1] = 0

    frank = N_HEAT * 2 + N_HR
    for i in range(N_SOLAR):
        if individual[frank + 2 * i + 1] < 0:
            individual[frank + 2 * i + 1] = 0

    sharePlants = 0
    for i in range(N_HEAT):
        sharePlants += individual[2 * i + 1]
    if abs(sharePlants - 1) > 1E-3:
        valid = False

    shareSolar = 0
    nSol = 0
    for i in range(N_SOLAR):
        nSol += individual[frank + 2 * i]
        shareSolar += individual[frank + 2 * i + 1]
    if nSol > 0 and abs(shareSolar - 1) > 1E-3:
        valid = False

    heating_part = 2 * N_HEAT + N_HR + 2 * N_SOLAR + INDICES_CORRESPONDING_TO_DHN
    for i in range(N_COOL):
        if individual[heating_part + 2 * i] > 0 and individual[heating_part + 2 * i + 1] < 0.01:
            oldValue = individual[heating_part + 2 * i + 1]
            shareGain = oldValue - 0.01
            individual[heating_part + 2 * i + 1] = 0.01

            for rank in range(N_COOL):
                if individual[heating_part + 2 * rank] > 0 and i != rank:
                    individual[heating_part + 2 * rank + 1] += individual[heating_part + 2 * rank + 1] / (1 - oldValue) * shareGain
        elif individual[heating_part + 2*i] == 0:
            individual[heating_part + 2 * i + 1] = 0

    sharePlants = 0
    for i in range(N_COOL):
        sharePlants += individual[heating_part + 2 * i + 1]
    if abs(sharePlants - 1) > 1E-3:
        valid = False

    if not valid:
        newInd = generation.generate_main(nBuildings, config)

        L = (N_HEAT + N_SOLAR) * 2 + N_HR
        for i in range(L):
            individual[i] = newInd[i]

    return individual


def calc_master_to_slave_variables(individual, Q_heating_max_W, Q_cooling_max_W, building_names, ind_num, gen):
    """
    This function reads the list encoding a configuration and implements the corresponding
    for the slave routine's to use
    :param individual: list with inidividual
    :param Q_heating_max_W:  peak heating demand
    :param locator: locator class
    :param gv: global variables class
    :type individual: list
    :type Q_heating_max_W: float
    :type locator: string
    :type gv: class
    :return: master_to_slave_vars : class MasterSlaveVariables
    :rtype: class
    """
    # initialise class storing dynamic variables transfered from master to slave optimization
    master_to_slave_vars = slave_data.SlaveData()
    configkey = "".join(str(e)[0:4] for e in individual)

    DHN_barcode, DCN_barcode, DHN_configuration, DCN_configuration = sFn.individual_to_barcode(individual, building_names)
    configkey = configkey[:-2*len(DHN_barcode)] + hex(int(str(DHN_barcode),2)) + hex(int(str(DCN_barcode),2))
    master_to_slave_vars.configKey = configkey
    master_to_slave_vars.number_of_buildings_connected_heating = DHN_barcode.count("1") # counting the number of buildings connected in DHN
    master_to_slave_vars.number_of_buildings_connected_cooling = DCN_barcode.count("1") # counting the number of buildings connectedin DCN
    master_to_slave_vars.individual_number = ind_num
    master_to_slave_vars.generation_number = gen

    Q_heating_nom_W = Q_heating_max_W * (1 + Q_MARGIN_FOR_NETWORK)
    Q_cooling_nom_W = Q_cooling_max_W * (1 + Q_MARGIN_FOR_NETWORK)
    
    # Heating systems

    #CHP units with NG & furnace with biomass wet
    if individual[0] == 1 or individual[0] == 3:
        if FURNACE_ALLOWED == True:
            master_to_slave_vars.Furnace_on = 1
            master_to_slave_vars.Furnace_Q_max = max(individual[1] * Q_heating_nom_W, Q_MIN_SHARE * Q_heating_nom_W)
            master_to_slave_vars.Furn_Moist_type = "wet"
        elif CC_ALLOWED == True:
            master_to_slave_vars.CC_on = 1
            master_to_slave_vars.CC_GT_SIZE = max(individual[1] * Q_heating_nom_W * 1.3, Q_MIN_SHARE * Q_heating_nom_W * 1.3)
            #1.3 is the conversion factor between the GT_Elec_size NG and Q_DHN
            master_to_slave_vars.gt_fuel = "NG"

    #CHP units with BG& furnace with biomass dry
    if individual[0] == 2 or individual[0] == 4:
        if FURNACE_ALLOWED == True:
            master_to_slave_vars.Furnace_on = 1
            master_to_slave_vars.Furnace_Q_max = max(individual[1] * Q_heating_nom_W, Q_MIN_SHARE * Q_heating_nom_W)
            master_to_slave_vars.Furn_Moist_type = "dry"
        elif CC_ALLOWED == True:
            master_to_slave_vars.CC_on = 1
            master_to_slave_vars.CC_GT_SIZE = max(individual[1] * Q_heating_nom_W * 1.5, Q_MIN_SHARE * Q_heating_nom_W * 1.5)
            #1.5 is the conversion factor between the GT_Elec_size BG and Q_DHN
            master_to_slave_vars.gt_fuel = "BG"

    # Base boiler NG
    if individual[2] == 1:
        master_to_slave_vars.Boiler_on = 1
        master_to_slave_vars.Boiler_Q_max = max(individual[3] * Q_heating_nom_W, Q_MIN_SHARE * Q_heating_nom_W)
        master_to_slave_vars.BoilerType = "NG"

    # Base boiler BG
    if individual[2] == 2:
        master_to_slave_vars.Boiler_on = 1
        master_to_slave_vars.Boiler_Q_max = max(individual[3] * Q_heating_nom_W, Q_MIN_SHARE * Q_heating_nom_W)
        master_to_slave_vars.BoilerType = "BG"

    # peak boiler NG
    if individual[4] == 1:
        master_to_slave_vars.BoilerPeak_on = 1
        master_to_slave_vars.BoilerPeak_Q_max = max(individual[5] * Q_heating_nom_W, Q_MIN_SHARE * Q_heating_nom_W)
        master_to_slave_vars.BoilerPeakType = "NG"

    # peak boiler BG
    if individual[4] == 2:
        master_to_slave_vars.BoilerPeak_on = 1
        master_to_slave_vars.BoilerPeak_Q_max = max(individual[5] * Q_heating_nom_W, Q_MIN_SHARE * Q_heating_nom_W)
        master_to_slave_vars.BoilerPeakType = "BG"

    # lake - heat pump
    if individual[6] == 1  and HP_LAKE_ALLOWED == True:
        master_to_slave_vars.HP_Lake_on = 1
        master_to_slave_vars.HPLake_maxSize = max(individual[7] * Q_heating_nom_W, Q_MIN_SHARE * Q_heating_nom_W)

    # sewage - heatpump
    if individual[8] == 1 and HP_SEW_ALLOWED == True:
        master_to_slave_vars.HP_Sew_on = 1
        master_to_slave_vars.HPSew_maxSize = max(individual[9] * Q_heating_nom_W, Q_MIN_SHARE * Q_heating_nom_W)

    # Gwound source- heatpump
    if individual[10] == 1 and GHP_ALLOWED == True:
        master_to_slave_vars.GHP_on = 1
        GHP_Qmax = max(individual[11] * Q_heating_nom_W, Q_MIN_SHARE * Q_heating_nom_W)
        master_to_slave_vars.GHP_number = GHP_Qmax / GHP_HMAX_SIZE

    # heat recovery servers and compresor
    irank = N_HEAT * 2
    master_to_slave_vars.WasteServersHeatRecovery = individual[irank]
    master_to_slave_vars.WasteCompressorHeatRecovery = individual[irank + 1]

    # Solar systems
    shareAvail = 1  # all buildings in the neighborhood are connected to the solar potential

    irank = N_HEAT * 2 + N_HR

    heating_block = N_HEAT * 2 + N_HR + N_SOLAR * 2 + INDICES_CORRESPONDING_TO_DHN
    master_to_slave_vars.DHN_supplyunits = DHN_configuration
    # cooling systems

    # Lake Cooling
    if individual[heating_block] == 1 and LAKE_COOLING_ALLOWED is True:
        master_to_slave_vars.Lake_cooling_on = 1
        master_to_slave_vars.Lake_cooling_size = max(individual[heating_block + 1] * Q_cooling_nom_W, Q_MIN_SHARE * Q_cooling_nom_W)

    # VCC Cooling
    if individual[heating_block + 2] == 1 and VCC_ALLOWED is True:
        master_to_slave_vars.VCC_on = 1
        master_to_slave_vars.VCC_cooling_size = max(individual[heating_block + 3] * Q_cooling_nom_W, Q_MIN_SHARE * Q_cooling_nom_W)

    # Absorption Chiller Cooling
    if individual[heating_block + 4] == 1 and ABSORPTION_CHILLER_ALLOWED is True:
        master_to_slave_vars.Absorption_Chiller_on = 1
        master_to_slave_vars.Absorption_chiller_size = max(individual[heating_block + 5] * Q_cooling_nom_W, Q_MIN_SHARE * Q_cooling_nom_W)

    # Storage Cooling
    if individual[heating_block + 6] == 1 and STORAGE_COOLING_ALLOWED is True:
        if (individual[heating_block + 2] == 1 and VCC_ALLOWED is True) or (individual[heating_block + 4] == 1 and ABSORPTION_CHILLER_ALLOWED is True):
            master_to_slave_vars.storage_cooling_on = 1
            master_to_slave_vars.Storage_cooling_size = max(individual[heating_block + 7] * Q_cooling_nom_W, Q_MIN_SHARE * Q_cooling_nom_W)
            if master_to_slave_vars.Storage_cooling_size > STORAGE_COOLING_SHARE_RESTRICTION * Q_cooling_nom_W:
                master_to_slave_vars.Storage_cooling_size = STORAGE_COOLING_SHARE_RESTRICTION * Q_cooling_nom_W

    master_to_slave_vars.DCN_supplyunits = DCN_configuration
    master_to_slave_vars.SOLAR_PART_PV = max(individual[irank] * individual[irank + 1] * individual[irank + 8] * shareAvail,0)
    master_to_slave_vars.SOLAR_PART_PVT = max(individual[irank + 2] * individual[irank + 3] * individual[irank + 8] * shareAvail,0)
    master_to_slave_vars.SOLAR_PART_SC_ET = max(individual[irank + 4] * individual[irank + 5] * individual[irank + 8] * shareAvail,0)
    master_to_slave_vars.SOLAR_PART_SC_FP = max(individual[irank + 6] * individual[irank + 7] * individual[irank + 8] * shareAvail,0)

    return master_to_slave_vars


def checkNtw(individual, DHN_network_list, DCN_network_list, locator, gv, config, building_names):
    """
    This function calls the distribution routine if necessary
    :param individual: network configuration considered
    :param DHN_network_list: list of DHN configurations previously encounterd in the master
    :param locator: path to the folder
    :type individual: list
    :type DHN_network_list: list
    :type locator: string
    :return: None
    :rtype: Nonetype
    """
    DHN_barcode, DCN_barcode, DHN_configuration, DCN_configuration = sFn.individual_to_barcode(individual, building_names)

    if not (DHN_barcode in DHN_network_list) and DHN_barcode.count("1") > 0:
        DHN_network_list.append(DHN_barcode)

        total_demand = sFn.createTotalNtwCsv(DHN_barcode, locator)
        building_names = total_demand.Name.values

        # Run the substation and distribution routines
        sMain.substation_main(locator, total_demand, building_names, DHN_configuration, DCN_configuration, Flag=True)

        nM.network_main(locator, total_demand, building_names, config, gv, DHN_barcode)


    if not (DCN_barcode in DHN_network_list) and DCN_barcode.count("1") > 0:
        DCN_network_list.append(DCN_barcode)

        total_demand = sFn.createTotalNtwCsv(DCN_barcode, locator)
        building_names = total_demand.Name.values

        # Run the substation and distribution routines
        sMain.substation_main(locator, total_demand, building_names, DHN_configuration, DCN_configuration, Flag=True)

        nM.network_main(locator, total_demand, building_names, config, gv, DCN_barcode)

def epsIndicator(frontOld, frontNew):
    """
    This function computes the epsilon indicator
    :param frontOld: Old Pareto front
    :param frontNew: New Pareto front
    :type frontOld: list
    :type frontNew:list
    :return: epsilon indicator between the old and new Pareto fronts
    :rtype: float
    """
    epsInd = 0
    firstValueAll = True

    for indNew in frontNew:
        tempEpsInd = 0
        firstValue = True

        for indOld in frontOld:
            (aOld, bOld, cOld) = indOld.fitness.values
            (aNew, bNew, cNew) = indNew.fitness.values
            compare = max(aOld-aNew, bOld-bNew, cOld-cNew)

            if firstValue:
                tempEpsInd = compare
                firstValue = False

            if compare < tempEpsInd:
                tempEpsInd = compare

        if firstValueAll:
            epsInd = tempEpsInd
            firstValueAll = False

        if tempEpsInd > epsInd:
            epsInd = tempEpsInd

    return epsInd