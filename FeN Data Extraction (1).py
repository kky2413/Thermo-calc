# -*- coding: utf-8 -*-

import math
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from tc_python import *
from tc_python.quantity import SiteFractionOfComponentInPhase
from tc_python.quantity_factory import ThermodynamicQuantity
#from tc_python.utils import ALL_PHASES



with TCPython() as session:
    spB = str("Fe")     # Balanced chemical species in the alloy
    sp1 = str("N")     # 1st chemical species in the alloy

    spV = str("Va")     # Vaccancy used in sublattice logic


    phases = ["BCC_A2", "HCP_A3", "FE4N_LP1"] # Phase for which gibbs free energy need to be extracted    ["BCC_A2", "HCP_A3", "FE4N_LP1"]
    for ph in phases:
        if ph == "BCC_A2":
            maxminDiff = float(0.75) # 각 상별 molefraction 조건
        #elif ph == "BCC_B2":
        #    maxminDiff = float(0.33)
        elif ph == "HCP_A3":
            maxminDiff = float(0.33)
        elif ph == "FE4N_LP1":
            maxminDiff = float(0.2)

        system = session.select_database_and_elements("TCFE12.0", [spB, sp1]).get_system()
        calc = system.with_single_equilibrium_calculation().disable_global_minimization()
        
        dc = float(0.001)
        
        c1_start = float(0.00)
        c1_end   = float(c1_start + maxminDiff + dc)

        cA_start = float(0.00)
        cA_end   = float(cA_start + 1 - maxminDiff + dc)

        

        tempStart = int(813)
        tempStop = int(813)
        tempIntv = int(10)
        temp = np.arange(tempStart,tempStop+tempIntv,tempIntv)
        

        #print("[TempStart : TempIntv : TempStop] : ", tempStart,tempIntv,tempStop)
        print("For : ",sp1)
        print(c1_start, c1_end)

        N_c1 = int((c1_end - c1_start)/dc)

        print("N_c1: ",N_c1)
        c1_loops = np.arange(0,N_c1,1)

        for tempI in temp:
            (calc
                .set_condition("T", tempI + 0.0)
                .set_phase_to_suspended(ALL_PHASES)
                .set_phase_to_entered(ph, 1.0))

            file_name = str(f"gibbsFree_{ph}_") + str(tempI) + str(".dat")
            print("Temperature: ",tempI)

            with open(file_name,"w") as fp:
                #cB_content = float( 1 - maxminDiff) 
                #calc.set_condition(ThermodynamicQuantity.mole_fraction_of_a_component(spB), cB_content)

                #N_c2 = int(N_c1 - (c1_content/dc) ) # C1_content 고정으로 인한 loop횟수 차감
                #N_c1 = int(N_c2 - (c2_content/dc) ) # C2_content 고정으로 인한 loop횟수 차감

                N_c1 = int((c1_end - c1_start)/dc)
                c1_loops = np.arange(0,N_c1,1)
                
                for c1_index in c1_loops:

                    c1_content = float(c1_start + c1_index * dc)
                    calc.set_condition(ThermodynamicQuantity.mole_fraction_of_a_component(sp1), c1_content)
                    result = calc.calculate()

                    gibbs_energy = float(result.get_value_of(ThermodynamicQuantity.gibbs_energy_of_a_phase(ph)))


                    # Site fraction of blanced species in phase p at sublattice 1
                    siteFrac_pcB_1 = result.get_value_of(SiteFractionOfComponentInPhase(ph,spB,1))  # Fe

                    # Site fraction of c1 in phase p at sublattice 2
                    siteFrac_pc1_2 = result.get_value_of(SiteFractionOfComponentInPhase(ph,sp1,2))   # C

                    # Site fraction of Va in phase p at sublattice 2
                    siteFrac_pcVa_2 = result.get_value_of(SiteFractionOfComponentInPhase(ph,spV,2)) 

                            #print("GM: ",c1_content, c2_content, gibbs_energy)
                            #print("-------------")

                    out_string  = ""
                    #out_string += " "  + str(cA_content)
                    out_string += " "  + str(c1_content)
                    out_string += " "  + str(siteFrac_pcB_1)
                    out_string += " "  + str(siteFrac_pc1_2)
                    out_string += " "  + str(siteFrac_pcVa_2)
                    out_string += " "  + str(gibbs_energy)
                    out_string += "\n"

                    fp.write(out_string)