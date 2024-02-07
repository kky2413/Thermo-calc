import os
import numpy as np
import matplotlib.pyplot as plt

# 수정된 현재 작업 디렉토리
comp = ['Cr', 'C', 'N']
for comp_name in comp:


    current_directory = rf"C:\Users\a2942\OneDrive\바탕 화면\thermo clc\FeCrCN 데이터\{comp_name}"

    def extract_values(file_path): 
        x_axis_values1 = []
        y_axis_values1 = []

        try:
            with open(file_path, 'r') as file:
                for line in file:
                    columns = line.split()
                    if comp_name == 'Cr':
                        x_axis_values1.append(float(columns[0])) # 4번째 열 값 저장  N molefraction
                    elif comp_name == 'N':
                        x_axis_values1.append(float(columns[2]))
                    elif comp_name == 'C':
                        x_axis_values1.append(float(columns[2])) #나중에 데이터 고치면 1로 바꿔야됨
                    y_axis_values1.append(float(columns[8]))   # 7번째 열 값 저장  GibbsEnergy


        except FileNotFoundError:
            print(f'{file_path}을(를) 찾을 수 없습니다.')

        return x_axis_values1, y_axis_values1


    phases = ["BCC_A2", "HCP_A3", "FE4N_LP1"]
    for ph in phases:
        file_pattern = f'gibbsFree_{ph}_813.dat'
        ascii_file_paths = [os.path.join(current_directory, file_pattern)]

        for file_path in ascii_file_paths:
            x_axis_list1,y_axis_list1 = extract_values(file_path)
            plt.plot(x_axis_list1,y_axis_list1, label=(f'{ph}'))

            

        # 플롯 생성 및 표시
    plt.xlabel(f'{comp_name} mole fraction')
    plt.ylabel('Gibbs Energy')
    plt.xlim(0,0.8)
    plt.title('extraction Gibbs Energy')
    plt.legend()
    output_image_path = os.path.join(current_directory, 'plot.png')

    plt.show()



    #site fraction은 4개씩 존재
    #하나의 조성을 고정하고 한 상에서의 site fraction을 구해야된다?
    #site fraction 나오는 그림이 9개가 나와야되나??
    #내일 와서 피피티 만들고 site fraction plot하는 코드 짜기
    #그림이 9개가 나와야 함.......
