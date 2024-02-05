# real case make molefraction array
# Convert molefraction to siteFraction
# Should consider the mixing in all the sublattice
# Interstitial compnents are phase 2
# 1상이 단일일 경우 confidential entropy 가 없음
#(A,B)k(C,D,F)l k,l:stoichmetric coefficient
# per mole formula unit\
# 입실론 감마프라임 상 추가해서 완성하기

import os
import numpy as np
import matplotlib.pyplot as plt
import math

# Set Condition
T = 813 # temperature(K)

current_directory = rf"C:\Users\a2942\OneDrive\바탕 화면\thermo clc\DATABASE"  # 현재 경로 설정

def extract_values(file_path):  # 값 추출 함수
    x_axis_values1 = []
    y_axis_values1 = []

    try:
        with open(file_path, 'r') as file:
            for line in file:
                columns = line.split()
                x_axis_values1.append(float(columns[0]))   # 2번째 열 값 저장  N molefraction
                y_axis_values1.append(float(columns[4]))   # 6번째 열 값 저장  GibbsEnergy

    except FileNotFoundError:
        print(f'{file_path}을(를) 찾을 수 없습니다.')

    return x_axis_values1, y_axis_values1

# Magnetic contribution 
def Magnetic_contribution(T):
    Tc= 1043.0
    p = 0.4
    tau = T/Tc
    D = 518/1125 + 11692/15975 * (1/p-1)
    gtau = 1.0 - (79/140/(p*tau) + 474/497 * (1/p-1)*tau**3/6+tau**9/135 + tau**15/600)/D
    B0 = 2.22
    R = 8.3144
    Gmag = R*T*math.log(B0+1) * gtau
    return Gmag


def GHSERCC(T):
    result = -17368.441+170.73*T-24.3*T*math.log(T)-4.723E-04*T**2+2562600*T**(-1)-2.643E+08*T**(-2)+1.2E+10*T**(-3)
    return result

def GDIACC(T):
    result = -16359.441+175.61*T-24.31*T*math.log(T)-4.723E-04*T**2+2698000*T**(-1)-2.61E+08*T**(-2)+1.11E+10*T**(-3)
    return result

def GHSERFE(T):
    result = 1225.7+124.134*T-23.5143*T*math.log(T)-0.00439752*T**2-5.8927E-08*T**3+77359*T**(-1)
    return result

def GFCCFE(T):
    result = -1462.4+8.282*T-1.15*T*math.log(T)+6.4E-04*T**2+GHSERFE(T)
    return result

def GHCPFE(T):
    result = -3705.78+12.591*T-1.15*T*math.log(T)+6.4E-04*T**2+GHSERFE(T)
    return result

def GLIQFE(T):
    result = 12040.17-6.55843*T+GHSERFE(T)-3.67516E-21*T**7
    return result

def GHSERNN(T):
    result = -3750.675-9.45425*T-12.7819*T*math.log(T)-0.00176686*T**2+2.681E-09*T**3-32374*T**(-1)
    return result

phases = ['BCC_A2', 'HCP_A3', 'FE4N_LP1']
for phase in phases:

    if phase == 'BCC_A2':
        k = 1  # k and l are stoichimetric coefficient
        l = 3

        G_Fe_N = GHSERFE(T) + 3 * GHSERNN(T) + 93562 + 165.07 * T
        G_Fe_Va = GHSERFE(T)
        G_Ex = 0 #BCC_A2 has no excess term for N and C
        G_Mag = Magnetic_contribution(T)
        

    elif phase == 'HCP_A3':
        k = 1  # k and l are stoichimetric coefficient
        l = 0.5

        G_Fe_N = GHSERFE(T) + 0.5 * GHSERNN(T)-13863 + 40.2123*T
        G_Fe_Va = GHCPFE(T)
        G_Fe_Va_N0 =  8186-18.127*T   # 0 function of excess function
        G_Fe_Va_N1 = -24378+24.959*T  # 1           ''
        G_Ex =  (G_Fe_Va_N0 + G_Fe_Va_N1)    # 과잉 함수 계산시 site fraction 고려??(현재 미고려)
        G_Mag = 0
        

    elif phase == 'FE4N_LP1':
        k = 4
        l = 1

        G_Fe_N = 4*GHSERFE(T) + GHSERNN(T) -37744 + 72.786 * T
        G_Fe_Va = 4*GHSERFE(T) + 12066 + 3.691*T
        G_Ex = 0 # no excess function
        G_Mag = 0
        
    maxDiff = l / (k+l)
    
    
    file_pattern = f'gibbsFree_{phase}_813.dat'
    ascii_file_paths = [os.path.join(current_directory, file_pattern)]
        
    for file_path in ascii_file_paths:
        x_axis_list1,y_axis_list1 = extract_values(file_path)
        plt.plot(x_axis_list1,y_axis_list1, label=('DB'))


        C_Fe_start =0.00
        interval = 0.01
        C_Fe_end = C_Fe_start + maxDiff + interval

        x_axis = np.arange(C_Fe_start,C_Fe_end,interval)
        y_axis =[]

        k_f = k/(k+l) #k and l are stoichimetric coefficient
        k_f = round(k_f, 2) # 소수점 처리 넷째 자리에서 반올림
        l_f = l/(k+l)
        l_f = round(l_f, 2)


        for comp in x_axis:
            s_N = ((k_f*comp)/(l_f*(1-comp)))  # Site Fraction of N
            s_Va = 1 - s_N # Site Fraction of Va
            s_Va = abs(s_Va)

            if s_Va  == 0.0:
                s_Va = 0.005
                s_N = 0.995
            elif s_N  == 0.0:
                s_N = 0.005
                s_Va = 0.995

            s_Fe = float(1.0)           # Fe 단일로 ph1 의 sitefraction=1
            F_unit = k*s_Fe + l*s_N     # should consider formula-unit permole  1/(k*s_a+l*sitefraction) 
            G_result = (s_N * G_Fe_N) + (s_Va * G_Fe_Va) + (l*(8.3144)*T)*(s_N * math.log(s_N) + s_Va * math.log(s_Va)) + G_Ex + G_Mag
            G_result_f = G_result / F_unit
            y_axis.append(G_result_f)

        plt.plot(x_axis, y_axis, color ='g',label =('cal'))
        plt.title(f'{phase} Comparing Gibbs Energy ')
        plt.xlim(0,float(l_f)) #  stoichimetric coefficient로 정해진 고용한도로 x축 제한 
        plt.xlabel('N molefraction')
        plt.ylabel('GibbsEnergy')
        plt.legend()
        plt.show()
