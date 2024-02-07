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
    spA = str("Cr")
    spB = str("Fe")     # Balanced chemical species in the alloy
    sp1 = str("N")     # 1st chemical species in the alloy
    sp2 = str("C")     # 2nd chemical species in the alloy

    spV = str("Va")     # Vaccancy used in sublattice logic

    ph = str("BCC_A2")  # Phase for which gibbs free energy need to be extracted

    system = session.select_database_and_elements("TCFE12.0", [spA, spB, sp1, sp2]).get_system()
    calc = system.with_single_equilibrium_calculation().disable_global_minimization()
    
    dc = float(0.0375)
    maxminDiff = float(0.75)
    
    c1_start = float(0.00)
    c1_end   = float(c1_start + maxminDiff + dc)

    # 100 ppm, alloy composition for fixed mole fraction of Sr
    c2_start = float(0.00)
    c2_end   = float( c2_start + maxminDiff + dc)

    cA_start = float(0.00)
    cA_end   = float(cA_start + 1 - maxminDiff + dc)

    

    tempStart = int(813)
    tempStop = int(813)
    tempIntv = int(10)
    temp = np.arange(tempStart,tempStop+tempIntv,tempIntv)
    

    #print("[TempStart : TempIntv : TempStop] : ", tempStart,tempIntv,tempStop)
    print("For : ",sp1)
    print(c1_start, c1_end)
    print("For : ",sp2)
    print(c2_start, c2_end)
    print("For : ",spA)
    print(cA_start, cA_end)

    N_c1 = int((c1_end - c1_start)/dc)
    print("N_c1: ",N_c1)
    c1_loops = np.arange(0,N_c1,1)

    for tempI in temp:
        (calc
            .set_condition("T", tempI + 0.0)
            .set_phase_to_suspended(ALL_PHASES)
            .set_phase_to_entered(ph, 1.0))

        file_name = str("gibbsFree_") + str(tempI) + str(".dat")
        print("Temperature: ",tempI)

        with open(file_name,"w") as fp:

            index1 = 0
            for c1_index in c1_loops:

                c1_content = float(c1_start + c1_index * dc)
                calc.set_condition(ThermodynamicQuantity.mole_fraction_of_a_component(sp1), c1_content)
            
                N_c2 = int(N_c1 - index1)
                #print("N_c2: ",N_c2)
                c2_loops = np.arange(0,N_c2,1)
                index1 = index1 + 1
                
                for c2_index in c2_loops:               

                    c2_content = float(c2_start +  c2_index * dc)
                    calc.set_condition(ThermodynamicQuantity.mole_fraction_of_a_component(sp2), c2_content)
                    #print("iterations: ",ph, c1_content, c2_content)

                    N_cA = int((cA_end - cA_start)/dc - index1) 
                    #N_cA = int(((cA_end - cA_start)/dc) - index1) 
                    #print("N_c2: ",N_c2)
                    cA_loops = np.arange(0,N_cA,1)
                    

                    for cA_index in cA_loops:               

                        cA_content = float(cA_start + cA_index * dc)
                        calc.set_condition(ThermodynamicQuantity.mole_fraction_of_a_component(spA), cA_content)
                        #print("iterations: ",ph, c1_content, c2_content)


                        result = calc.calculate()

                        gibbs_energy = float(result.get_value_of(ThermodynamicQuantity.gibbs_energy_of_a_phase(ph)))
                        


                        # site fraction of species in phase p at sublattice 1
                        siteFrac_pcA_1 = result.get_value_of(SiteFractionOfComponentInPhase(ph,spA,1))

                        # Site fraction of blanced species in phase p at sublattice 1
                        siteFrac_pcB_1 = result.get_value_of(SiteFractionOfComponentInPhase(ph,spB,1))

                        # Site fraction of c1 in phase p at sublattice 2
                        siteFrac_pc1_2 = result.get_value_of(SiteFractionOfComponentInPhase(ph,sp1,2))
    
                        # Site fraction of c2 in phase p at sublattice 2
                        siteFrac_pc2_2 = result.get_value_of(SiteFractionOfComponentInPhase(ph,sp2,2))

                        # Site fraction of Va in phase p at sublattice 2
                        siteFrac_pcVa_2 = result.get_value_of(SiteFractionOfComponentInPhase(ph,spV,2))

                        #print("GM: ",c1_content, c2_content, gibbs_energy)
                        #print("-------------")

                        out_string  = ""
                        out_string += " "  + str(cA_content)
                        out_string += " "  + str(c1_content)
                        out_string += " "  + str(c2_content)
                        out_string += " "  + str(siteFrac_pcA_1)
                        out_string += " "  + str(siteFrac_pcB_1)
                        out_string += " "  + str(siteFrac_pc1_2)
                        out_string += " "  + str(siteFrac_pc2_2)
                        out_string += " "  + str(siteFrac_pcVa_2)
                        out_string += " "  + str(gibbs_energy)
                        out_string += "\n"

                        fp.write(out_string)