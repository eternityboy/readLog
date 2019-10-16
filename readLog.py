version = "1.5.1"

## Version 1.5.1
## Bug Fixes
## To-Do

## Version 1.5
## support for both OF branches
## To-Do
## suggestions? QT? FFT analyze?

## Version 1.4
## reference point changed to reference range
## To-Do
## support for both OF versions

## Version 1.3
## OpenFOAM-v1806 support
## To-Do
## 

## Version 1.2
## added roll motion analize
## To-Do
## system argument - reference point 

## Version 1.1
## Whats New
## added possibilty to specify range of values to calculate mean
## added mean value to plot
## re-aranged part of the code with numpy
## SixDoF support
## Volume Phase
####################################################################
import sys
import numpy as np
import math
import re
import matplotlib.pyplot as plt
import time
#from scipy import signal
#from scipy.fftpack import fft
from functions import *

fileName = sys.argv[1]#.strip().split('/')
inputName = sys.argv[1]

#print(sys.argv[1])
#fileName = sys.argv[1].strip().split('.').split('/')
#print(fileName)
#
#inputName = fileName[-2]+'.'+fileName[-1]


try:
    ref_point_s
except NameError:
    rp = 10
else:
    rp = ref_point_s

def load_log():
    start_time = time.time()
    COL = 25
    global f_l, n_i
    f_l = file_len(inputName)
    n_i = num_iter(inputName,f_l)
    f = open(inputName, 'r').readlines()
    results = [[0 for x in range(n_i)] for y in range(COL)]
    print("--- Reading log file ---")
    print(sys.argv[1])
    i = 0
    loop = -1
    for i in range(f_l):
        if f[i].startswith('Time = '):
            loop += 1
            T  = f[i].strip().split(' ')
            results[0][loop] = T[2]
        if f[i].lstrip().startswith('Centre of mass:'):
            CoM = f[i].strip().split(' ')
            results[1][loop] = CoM[3].replace('(', '')
            results[2][loop] = CoM[4]
            results[3][loop] = CoM[5].replace(')', '')
        if f[i].lstrip().startswith('Linear velocity:'):
            LV = f[i].strip().split(' ')
            results[4][loop] = LV[2].replace('(', '')
            results[5][loop] = LV[3]
            results[6][loop] = LV[4].replace(')', '')
        if f[i].lstrip().startswith('Angular velocity:'):
            AV = f[i].strip().split(' ')
            results[7][loop] = AV[2].replace('(', '')
            results[8][loop] = AV[3]
            results[9][loop] = AV[4].replace(')', '')
        if f[i].lstrip().startswith('Execution time for'):
            MU = f[i].strip().split(' ')
            results[10][loop] = MU[5]
        if f[i].lstrip().startswith('ExecutionTime ='):
            Et = f[i].strip().split(' ')
            results[11][loop] = Et[2]
        if f[i].lower().lstrip().startswith('sum of forces'):
            if f[i+1].lower().lstrip().startswith('pressure'):
                Fp = re.sub(" +", " ", f[i+1]).strip()
                Fp = re.sub("[()]","",Fp).split(" ")
                results[12][loop] = Fp[2]
                results[13][loop] = Fp[3]
                results[14][loop] = Fp[4]
            if f[i+2].lower().lstrip().startswith('pressure'):
                Fp = re.sub(" +", " ", f[i+2]).strip()
                Fp = re.sub("[()]","",Fp).split(" ")
                results[12][loop] = Fp[2]
                results[13][loop] = Fp[3]
                results[14][loop] = Fp[4]
            if f[i+2].lower().lstrip().startswith('viscous'):
                Fv = re.sub(" +", " ", f[i+2]).strip()
                Fv = re.sub("[()]","",Fv).split(" ")
                results[15][loop] = Fv[2]
                results[16][loop] = Fv[3]
                results[17][loop] = Fv[4]
            if f[i+3].lower().lstrip().startswith('viscous'):
                Fv = re.sub(" +", " ", f[i+3]).strip()
                Fv = re.sub("[()]","",Fv).split(" ")
                results[15][loop] = Fv[2]
                results[16][loop] = Fv[3]
                results[17][loop] = Fv[4]
        if f[i].lower().lstrip().startswith('sum of moments'):
            if f[i+1].lower().lstrip().startswith('pressure'):
                Mp = re.sub(" +", " ", f[i+1]).strip()
                Mp = re.sub("[()]","",Mp).split(" ")
                results[18][loop] = Mp[2]
                results[19][loop] = Mp[3]
                results[20][loop] = Mp[4]
            if f[i+2].lower().lstrip().startswith('pressure'):
                Mp = re.sub(" +", " ", f[i+2]).strip()
                Mp = re.sub("[()]","",Mp).split(" ")
                results[18][loop] = Mp[2]
                results[19][loop] = Mp[3]
                results[20][loop] = Mp[4]
            if f[i+2].lower().lstrip().startswith('viscous'):
                Mv = re.sub(" +", " ", f[i+2]).strip()
                Mv = re.sub("[()]","",Mv).split(" ")
                results[21][loop] = Mv[2]
                results[22][loop] = Mv[3]
                results[23][loop] = Mv[4]
            if f[i+3].lower().lstrip().startswith('viscous'):
                Mv = re.sub(" +", " ", f[i+3]).strip()
                Mv = re.sub("[()]","",Mv).split(" ")
                results[21][loop] = Mv[2]
                results[22][loop] = Mv[3]
                results[23][loop] = Mv[4]
        if f[i].lstrip().startswith('Phase-1 volume fraction =') and f[i-1].lstrip().startswith('smoothSolver:'):
            vFr = f[i].strip().split(' ')
            results[24][loop] = vFr[4]
    print("--- Reading log file complete: %s seconds ---" % (time.time() - start_time))
    print("--- Creating NumPy array ---")
    results = np.array(results, dtype=np.float32) #.transpose()
    print("--- Creating NumPy array complete: %s seconds ---" % (time.time() - start_time))
    print("--- Creating global variables ---")
    global r_t,ETfMU,ET
    global CoM_x,CoM_y,CoM_z
    global LV_x,LV_y,LV_z
    global AV_x,AV_y,AV_z
    global Fp_x,Fp_y,Fp_z,Fv_x,Fv_y,Fv_z
    global Mp_x,Mp_y,Mp_z,Mv_x,Mv_y,Mv_z
    global VF
    global F_x,F_y,F_z,M_x,M_y,M_z
    global rp
    r_t = results[0]
    CoM_x = results[1]
    CoM_y = results[2]
    CoM_z = results[3]
    LV_x = results[4]
    LV_y = results[5]
    LV_z = results[6]
    AV_x = results[7]
    AV_y = results[8]
    AV_z = results[9]
    ETfMU = results[10]
    ET = results[11]
    Fp_x = results[12]
    Fp_y = results[13]
    Fp_z = results[14]
    Fv_x = results[15]
    Fv_y = results[16]
    Fv_z = results[17]
    Mp_x = results[18]
    Mp_y = results[19]
    Mp_z = results[20]
    Mv_x = results[21]
    Mv_y = results[22]
    Mv_z = results[23]
    VF = results[24]
    F_x = Fp_x + Fv_x
    F_y = Fp_y + Fv_y
    F_z = Fp_z + Fv_z
    M_x = Mp_x + Mv_x
    M_y = Mp_y + Mv_y
    M_z = Mp_z + Mv_z
    rp = [10,n_i-1]
    print("--- Creating global variables complete: %s seconds ---" % (time.time() - start_time))
    print("--- Reading motion values ---")
    if all(v1==0 for v1 in CoM_x) and all(v2==0 for v2 in CoM_z) :
        print("")
        print("No motions in the log file!")
    else:
        global sixDoF
        sixDoF = [[0 for x in range(n_i-1)] for y in range(5)]
        i = t = heave = pitch = roll = 0
        for i in range(len(r_t)-1):
            dT = float(r_t[i+1])-float(r_t[i])
            heave = heave + 0.5*(float(LV_z[i+1])+float(LV_z[i]))*dT
            pitch = pitch + 0.5*(float(AV_y[i+1])+float(AV_y[i]))*dT
            roll = roll + 0.5*(float(AV_x[i+1])+float(AV_x[i]))*dT
            t = t + dT
            sixDoF[0][i] = t
            sixDoF[1][i] = float(CoM_z[i])
            sixDoF[2][i] = heave
            sixDoF[3][i] = pitch
            sixDoF[4][i] = roll
            i += 1
        sixDoF = np.array(sixDoF, dtype=np.float32)
    print("--- Reading motion variables complete: %s seconds ---" % (time.time() - start_time))

#############################################################################
def time_consumption():
    N = len(rp)-2
    print("Total calculation time: %.2f m or %.1f h" % (ET[N]/60,ET[N]/3600))
    print("Calculation ration Execution Time / Real Time: %f" % (ET[N]/r_t[N]))
    print("Total time to update mesh: %f m and it is %.2f %% of total time" % (np.sum(ETfMU)/60,(ET[N]-np.sum(ETfMU))/ET[N]))

    sT = 1.0/3600.0
    plotgraph(r_t,ETfMU,0,0,'Execution time for mesh update','Time, s','Execution time for mesh update, s',rp)
    plotgraph(r_t,ET*sT,0,0,'Execution time, hours','Time, s','Execution Time, hours',rp)
    plotgraph(r_t,VF,0,0,'Phase-1 volume fraction','Time, s','Volume fraction',rp)

#############################################################################
def reference_point():
    global ref_point_s,ref_point_e
         
    plt
    plt.xlabel('Time, s')
    plt.ylabel('Total Force, N')
    plt.subplot(3,1,1)
    plt.title('F_x')
    plt.plot(r_t[50:-1],F_x[50:-1],'r-')
    plt.subplot(3,1,2)
    plt.title('F_y')
    plt.plot(r_t[50:-1],F_y[50:-1],'b-')
    plt.subplot(3,1,3)
    plt.title('F_z')
    plt.plot(r_t[50:-1],F_z[50:-1],'g-')  
    plt.show()
    
    ref_point_s = float(input("Enter reference point for time to START calculating Mean values:\n") or r_t[10])
    print ("Start time: %s" % ref_point_s)
    ref_point_s = find_nearest(r_t,ref_point_s)
    
    ref_point_e = float(input("Enter reference point for time to END calculating Mean values:\n") or r_t[-10])
    print ("End time: %s" % ref_point_e)
    ref_point_e = find_nearest(r_t,ref_point_e)  
    
############################################################################
def reference_pointOld():
    global ref_point
    while True:
        cs = input("Choose plane to set reference point (x default, y or z:\n").lower() or 'x'
        if cs == 'x':
           m_plane = F_x
           break
        elif cs == 'y':
           m_plane = F_y
           break
        elif cs == 'z':
           m_plane = F_z
           break
        else:
           print("Invalid command\n")

    
    plotgraph(r_t,m_plane,0,0,'F_{'+cs+'}','Time, s','Total Force, N')
    
    ref_point = input("Enter reference point for time to start calculating Mean values:\n") or 0.1
    ref_point = float(ref_point)
    ref_point = find_nearest(r_t,ref_point)

#############################################################################

#############################################################################
def forces(a):
    print('FORCES')
    m_F_x = calc_mean(r_t,F_x,rp)
    m_F_y = calc_mean(r_t,F_y,rp)
    m_F_z = calc_mean(r_t,F_z,rp)
    
    m_Fp_x = calc_mean(r_t,Fp_x,rp)
    m_Fp_y = calc_mean(r_t,Fp_y,rp)
    m_Fp_z = calc_mean(r_t,Fp_z,rp)
    
    m_Fv_x = calc_mean(r_t,Fv_x,rp)
    m_Fv_y = calc_mean(r_t,Fv_y,rp)
    m_Fv_z = calc_mean(r_t,Fv_z,rp)

    if a != 0:
        print("X: Mean total force: %s, Pressure force: %s, Viscous force: %s" % (m_F_x, m_Fp_x, m_Fv_x))
        print("Y: Mean total force: %s, Pressure force: %s, Viscous force: %s" % (m_F_y, m_Fp_y, m_Fv_y))
        print("Z: Mean total force: %s, Pressure force: %s, Viscous force: %s" % (m_F_z, m_Fp_z, m_Fv_z))
        print('############')
        print('Export view:')
        print("%s %s %s %s %s %s %s %s %s" % (m_F_x, m_F_y, m_F_z, m_Fp_x, m_Fp_y, m_Fp_z, m_Fv_x, m_Fv_y, m_Fv_z))
    else:
        print("%s %s %s %s %s %s %s %s %s" % (m_F_x, m_F_y, m_F_z, m_Fp_x, m_Fp_y, m_Fp_z, m_Fv_x, m_Fv_y, m_Fv_z))
    
    if a != 0:
        plotgraph(r_t,F_x,m_F_x,0,'F_{x}','Time, s','Total Force, N',rp)
        plotgraph(r_t,F_y,m_F_y,0,'F_{y}','Time, s','Total Force, N',rp)
        plotgraph(r_t,F_z,m_F_z,0,'F_{z}','Time, s','Total Force, N',rp)

#############################################################################
def moments(a):
    print('MOMENTS')
    m_M_x = calc_mean(r_t,M_x,rp)
    m_M_y = calc_mean(r_t,M_y,rp)
    m_M_z = calc_mean(r_t,M_z,rp)
    
    m_Mp_x = calc_mean(r_t,Mp_x,rp)
    m_Mp_y = calc_mean(r_t,Mp_y,rp)
    m_Mp_z = calc_mean(r_t,Mp_z,rp)
    
    m_Mv_x = calc_mean(r_t,Mv_x,rp)
    m_Mv_y = calc_mean(r_t,Mv_y,rp)
    m_Mv_z = calc_mean(r_t,Mv_z,rp)
    
    if a != 0:
        print("X: Mean total moment: %s, Pressure moment: %s, Viscous moment: %s" % (m_M_x, m_Mp_x, m_Mv_x))
        print("Y: Mean total moment: %s, Pressure moment: %s, Viscous moment: %s" % (m_M_y, m_Mp_y, m_Mv_y))
        print("Z: Mean total moment: %s, Pressure moment: %s, Viscous moment: %s" % (m_M_z, m_Mp_z, m_Mv_z))
        print('############')
        print('Export view:')
        print("%s %s %s %s %s %s %s %s %s" % (m_M_x, m_M_y, m_M_z, m_Mp_x, m_Mp_y, m_Mp_z, m_Mv_x, m_Mv_y, m_Mv_z))
    else:
        print("%s %s %s %s %s %s %s %s %s" % (m_M_x, m_M_y, m_M_z, m_Mp_x, m_Mp_y, m_Mp_z, m_Mv_x, m_Mv_y, m_Mv_z))
    
    if a != 0:
        plotgraph(r_t,M_x,m_M_x,0,'M_{x}','Time, s','Total Moment, N*m',rp)
        plotgraph(r_t,M_y,m_M_y,0,'M_{y}','Time, s','Total Moment, N*m',rp)
        plotgraph(r_t,M_z,m_M_z,0,'M_{z}','Time, s','Total Moment, N*m',rp)

#############################################################################
def phases():
    #rp = find_nearest(r_t,float(6.0))
    H_phase_x = calc_phase(sixDoF[0],sixDoF[2],rp)
    P_phase_x = calc_phase(sixDoF[0],sixDoF[3],rp)
    print("Heave amplitude: %s, Heave shift: %s," % (H_phase_x[0], H_phase_x[1]))
    print("Pitch amplitude: %s, Pitch shift: %s," % (P_phase_x[0], P_phase_x[1]))
    

#############################################################################
def motions(a):
    #N = np.size(sixDoF[0])
    #mean_period = sixDoF[0][-1]/N
    #FFT = fft_welch(r_t[rp:-1],sixDoF[2][rp:-1])
    #plotgraph(FFT[0],FFT[1],0,0,'FFT','Frequency','m**2/s')

    mean_heave = calc_mean_ampl(sixDoF[0],sixDoF[2],rp)
    mean_pitch = calc_mean_ampl(sixDoF[0],sixDoF[3],rp)
    mean_roll = calc_mean_ampl(sixDoF[0],sixDoF[4],rp)

    print("HEAVE mean value: %s, local maximum: %s, local minimum: %s" % (mean_heave[0],mean_heave[1],mean_heave[2]))
    print("PITCH mean value: %s, local maximum: %s, local minimum: %s" % (mean_pitch[0],mean_pitch[1],mean_pitch[2]))
    print("ROLL mean value: %s, local maximum: %s, local minimum: %s" % (mean_roll[0],mean_roll[1],mean_roll[2]))
    if a != 0:
        plotgraph(sixDoF[0],sixDoF[2]-mean_heave[3],mean_heave[0],1,'Heave','Time, s','Z_{g}, m',rp)
        plotgraph(sixDoF[0],sixDoF[3]-mean_pitch[3],mean_pitch[0],1,'Pitch','Time, s','Degree',rp)
        plotgraph(sixDoF[0],sixDoF[4]-mean_roll[3],mean_roll[0],1,'Roll','Time, s','Degree',rp)

#############################################################################
def wave_forces():
    print("--- This is Beta program! Ship values based on the Wigley III Delft experiments ---")

    while True:
        lL = input("Please set wave length ratio lambda/L = ") or '1'
        if is_number(lL) == True:
           lL = float(lL)
           break
        else:
           print("Invalid command\n")

    Fr = 0.3
    L = 3.0
    g = 9.81
    U = Fr*math.sqrt(g*L)
    vD = 0.078
    C33 = 6119.0
    C55 = 2874.0
    rho = 1000.0
    wHeight = 0.025
    wLambda = lL*L
    wK = 2*math.pi/wLambda
    wOmega = wK*U + math.sqrt(wK*g)
    
    i = t = X_1a = X_1b = X_3a = X_3b = X_5a = X_5b = 0
    for i in range(len(r_t)-1):
        dT = float(r_t[i+1])-float(r_t[i])
        X_1a = X_1a + 2.0/(r_t[-1]-r_t[rp])*(0.5*(float(F_x[i+1])+float(F_x[i])))*dT*math.cos(wOmega*r_t[i+1])
        X_1b = X_1b + 2.0/(r_t[-1]-r_t[rp])*(0.5*(float(F_x[i+1])+float(F_x[i])))*dT*math.sin(wOmega*r_t[i+1])
        X_3a = X_3a + 2.0/(r_t[-1]-r_t[rp])*(0.5*(float(F_z[i+1])+float(F_z[i])))*dT*math.cos(wOmega*r_t[i+1])
        X_3b = X_3b + 2.0/(r_t[-1]-r_t[rp])*(0.5*(float(F_z[i+1])+float(F_z[i])))*dT*math.sin(wOmega*r_t[i+1])
        X_5a = X_5a + 2.0/(r_t[-1]-r_t[rp])*(0.5*(float(M_y[i+1])+float(M_y[i])))*dT*math.cos(wOmega*r_t[i+1])
        X_5b = X_5b + 2.0/(r_t[-1]-r_t[rp])*(0.5*(float(M_y[i+1])+float(M_y[i])))*dT*math.sin(wOmega*r_t[i+1])
        t = t + dT
        i += 1
    X1 = 2*math.sqrt(X_1a**2+X_1b**2)
    X3 = 2*math.sqrt(X_3a**2+X_3b**2)
    X5 = 2*math.sqrt(X_5a**2+X_5b**2)
    eps_X1 = math.atan(X_1b/X_1a)
    eps_X3 = math.atan(X_3b/X_3a)
    eps_X5 = math.atan(X_5b/X_5a)

    print("Longitudial wave force X1: %f, Wave force phase eps_X1: %f, dim_X1: %f" % (X1,eps_X1*180.0/math.pi,X1/(wK*wHeight*rho*g*vD)))
    print("HEAVE: Wave force X3: %f, Wave force phase eps_X3: %f, dim_X3: %f" % (X3,eps_X3*180.0/math.pi,X3/(wHeight*C33)))
    print("PITCH: Wave moment X5: %f, Wave moment phase eps_X5: %f, dim_dX5: %f" % (X5,eps_X5*180.0/math.pi,X5/(wK*wHeight*C55)))

#############################################################################
def menu():

    global ref_point_s, ref_point_e
    global rp
    
    if len(sys.argv) > 2:
        ref_point_s = find_nearest(r_t,float(sys.argv[3]))
        ref_point_e = find_nearest(r_t,float(sys.argv[4]))
        rp = [ref_point_s,ref_point_e]
        if int(sys.argv[2]) == 1:
            forces(0)
        elif int(sys.argv[2]) == 2:
            moments(0)
        else:
            print("Valid args:\n1 - goes for extract forces only,\n2 - goes for extract moments only,\n next arg min reference point,\n and after max reference point. or just nothing and have fun.")
        exit()
    
    print ("\n### ReadLog v%s for OpenFOAM v1806 ###" % version)
    if 'ref_point_s' in globals():
        rp = [ref_point_s,ref_point_e]
        print ("\n### Reference range is set to %s - %s ###" % (r_t[rp[0]],r_t[rp[1]]))
        ####print("Nearest point is %s and it is N %s in the list" % (r_t[ref_point], ref_point))

    
    strs = ('\n1. Time\n'
            '2. Set reference point\n'
            '3. Forces\n'
            '4. Moments\n'
            '5. Calculate longitidual and vertical wave forces\n'
            '6. Solid motions\n'
            '7. Calculate Phase shift\n'
            '8. Exit\n'
            '0. Update\n'
            '--- Your choice : ')
    choice = input(strs) or "9"
    return int(choice) 

load_log()

while True:
    choice = menu()
    if choice == 1:
        time_consumption()
    elif choice == 2:
        reference_point()
    elif choice == 3:
        forces(1)
    elif choice == 4:
        moments(1)
    elif choice == 5:
        wave_forces()
    elif choice == 6:
        motions()
    elif choice == 7:
        phases()
    elif choice == 0:
        load_log()
    elif choice == 8:
        break
    else:
        menu()
        
