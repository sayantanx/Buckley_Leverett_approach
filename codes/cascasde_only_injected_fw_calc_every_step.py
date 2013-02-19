#-------------------------------------------------------------------------------
# Name:         Random Walker with BL velocities, single-phase flow
# Purpose:      Use Buckley-Leverett theory to estimate phase velocities, and use it in the transtion probability
#
# Author:       Sayantan
#
# Created:      04/02/2013
# Copyright:    (c) Sayantan Bhowmik 2013
#-------------------------------------------------------------------------------

import numpy as np
import math
import os

def rel_perm(fname):
    rel_perm_f = fname
    rp_f = open(rel_perm_f,'r')
    temp1 = rp_f.readline()
    temp2 = temp1.split('-')
    temp3 = temp2[0].split()
    krw0,krg0 = float(temp3[0]), float(temp3[1])
    temp1 = rp_f.readline()
    temp2 = temp1.split('-')
    temp3 = temp2[0].split()
    n1,n2 = float(temp3[0]), float(temp3[1])
    temp1 = rp_f.readline()
    temp2 = temp1.split('-')
    temp3 = temp2[0].split()
    Swr,Sgr = float(temp3[0]), float(temp3[1])
    rp_f.close()
    rel_perm_table = [[] for i in range(3)]
    for s in range(0,1001):
        sat = s/1000.0
        rel_perm_table[0].append(sat)
        if sat<Sgr:
            kr=0
        else:
            if sat>1-Swr:
                kr = krg0
            else:
                kr = krg0 * ((sat-Sgr)/(1-Sgr-Swr))**n2
        rel_perm_table[1].append(kr)
    for s in range(1,1001):
        sat=1-s/1000.0
        if sat<Swr:
            kr=0
        else:
            if sat>1-Sgr:
                kr = krw0
            else:
                kr = krw0 * ((sat-Swr)/(1-Sgr-Swr))**n1
        rel_perm_table[2].append(kr)
    return(rel_perm_table)

def front_saturation(rel_perm_table,del_rho,perm,v_macro,mu_g,mu_w,d_curr,d_tar,length):
    import math, numpy as np
    frac_flow = [[],[]]
    perm = perm*9.8*(10**(-16))
    del_rho = del_rho*16.018
    v_macro = v_macro/(3600*24*3.28084)
    mu_g = mu_g*0.001
    mu_w = mu_w*0.001
    for s in range(1000):
        satn = round(float(s)/1000.0,3)
        krg = rel_perm_table[1][s]
        krw = rel_perm_table[2][s]
        sin_alpha = math.sin(math.atan((d_curr-d_tar)/length))
        if krg==0:
            f_g=0
        else:
            f_g = (1-perm*krw*del_rho*9.8*sin_alpha/(v_macro*mu_w))/(1+krw*mu_g/krg/mu_w)
        frac_flow[0].append(satn)
        frac_flow[1].append(f_g)
    dfg_dSg = [[],[]]

    for s in range(1,1000):
        satn = round(float(s)/1000.0,3)
        f1,f2 = frac_flow[1][s-1], frac_flow[1][s+1]
        slope_actual = (f2-f1)/(0.002)
        slope_target = frac_flow[1][s]/satn
        if abs(slope_actual-slope_target)<0.1*slope_target:
            value = satn
            break
    return(value)

def BL_velocity(rel_perm_table,del_rho,perm,v_macro,mu_g,mu_w,d_curr,d_tar,length,sat):
    import math
    frac_flow = [[],[]]
    if sat>0.99:
        sat=0.99
    satn = sat+0.01
    krg = rel_perm_table[1][rel_perm_table[0].index(round(satn,3))]
    krw = rel_perm_table[2][rel_perm_table[0].index(round(1-satn,3))-1]
    sin_alpha = math.sin(math.atan((d_curr-d_tar)/length))
    f_g1 = (1-perm*krw*del_rho*9.8*sin_alpha)/(1+krw*mu_g/krg/mu_w)
    satn = sat-0.01
    krg = rel_perm_table[1][rel_perm_table[0].index(round(satn,3))]
    krw = rel_perm_table[2][rel_perm_table[0].index(round(1-satn,3))-1]
    sin_alpha = math.sin(math.atan((d_curr-d_tar)/length))
    f_g2 = (1-perm*krw*del_rho*9.8*sin_alpha)/(1+krw*mu_g/krg/mu_w)
    slope = (f_g1-f_g2)/0.02
    return(slope)

def monte_carlo_sampling(values):
    import random
    sample = random.random()
    probabilities = [val/sum(values) for val in values]
    for i in range(len(values)):
        if sample<sum(probabilities[:i+1]):
            break
    return(i)

def write_saturations(count,_max):
    import struct
    NZ = len(count)
    NY = len(count[0])
    NX = len(count[0][0])
    fp = open('../results/' + results_file +'/'+ 'temp_satn2.bin','ab')
    for i in range(NZ):
        for j in range(NY):
            for k in range(NX):
                a = struct.pack('f',float(count[i][j][k])/float(_max[i][j][k]))
                fp.write(a)
    fp.close()

def report(current_model,NX,NY,NZ):
    import struct
    import numpy as np
    fp = open('../results/' + results_file +'/'+ 'temp_satn2.bin','rb')
    a = fp.read()
    fp.close()
    print len(a)
    num_values = len(a)/4
    num_time = num_values/NX/NY/NZ
    values = struct.unpack('f'*num_values,a)
    report_list = np.reshape(values,(num_time,NX*NY*NZ))
    report_list = np.transpose(report_list)
    fname = 'Model' + str(current_model) + '.txt'
    fp = open('../results/' + results_file +'/'+ fname,'w')
    fp.write('%s\n%d\n' %('All times',num_time))
    for i in range(num_time):
        fp.write('%s\n' %('Time '+str(i+1)))
    for i in range(NX*NY*NZ):
        for l in range(num_time):
            fp.write('%.3f\t' %(report_list[i][l]))
        fp.write('\n')
    fp.close()
    import os
    os.remove('../results/' + results_file +'/'+ 'temp_satn2.bin')


## ******************* Reading parameter file ************************
# par_file = raw_input('Parameter file name:')
import time
par_temp = raw_input('Run parameters file name: ')
par_file = '../data/' + par_temp
print 'Reading from ',par_file
par_f = open(par_file,'r')

# Results file
temp1 = par_f.readline()
temp2 = temp1.split('-')
results_file = temp2[0].strip()
if os.path.isdir('../results/'+results_file):
    del_dir = raw_input('Results directory exists, replace? (y/n)')
    if del_dir=='y':
        import shutil
        shutil.rmtree('../results/'+results_file)
    else:
        exit()

os.mkdir('../results/' + results_file)
p_write = open('../results/' + results_file +'/parameters.txt','w')
del temp1, temp2

# NX, NY, NZ
temp1 = par_f.readline()
temp2 = temp1.split('-')
temp3 = [int(item) for item in temp2[0].split()]
NX,NY,NZ = temp3[0],temp3[1],temp3[2]
p_write.write('Grid size: %d  %d  %d\n' %(NX,NY,NZ))
del temp1, temp2, temp3

# dX, dY, dZ
temp1 = par_f.readline()
temp2 = temp1.split('-')
temp3 = [int(item) for item in temp2[0].split()]
dx,dy,dz = temp3[0],temp3[1],temp3[2]
p_write.write('Grid block size: %d  %d  %d\n' %(dx,dy,dz))
del temp1, temp2, temp3

# Number of injectors
temp1 = par_f.readline()
temp2 = temp1.split('-')
num_inj = int(temp2[0])
p_write.write('Number of injectors: %d\n' %(num_inj))
del temp1, temp2

# Injection locations
inj_x, inj_y, inj_z = [], [], []
for i in range(num_inj):
    temp1 = par_f.readline()
    temp2 = temp1.split()
    inj_x.append(int(temp2[0]))
    inj_y.append(int(temp2[1]))
    inj_z.append(int(temp2[2]))
    p_write.write('%d  %d  %d\n' %(int(temp2[0]),int(temp2[1]),int(temp2[2])))

# Permeability file
temp1 = par_f.readline()
temp2 = temp1.split('-')
perm_file = temp2[0].strip()
p_write.write('Perm file: %s\n' %(perm_file))
del temp1, temp2

# Porosity file
temp1 = par_f.readline()
temp2 = temp1.split('-')
por_file = temp2[0].strip()
p_write.write('Porosity file: %s\n' %(por_file))
del temp1, temp2

# Depth file
temp1 = par_f.readline()
temp2 = temp1.split('-')
pres_file = temp2[0].strip()
p_write.write('Depth file: %s\n' %(pres_file))
del temp1, temp2

# ************ Fluid Properties *******************
p_write.write(' *********** Fluid Properties ***********\n')
# Relative permeability file
temp1 = par_f.readline()
temp2 = temp1.split('-')
rel_perm_f = temp2[0].strip()
p_write.write('Rel perm file: %s\n' %(rel_perm_f))
rp_f = open(rel_perm_f,'r')
temp1 = rp_f.readline()
temp2 = temp1.split('-')
temp3 = temp2[0].split()
krw0,krg0 = float(temp3[0]), float(temp3[1])
temp1 = rp_f.readline()
temp2 = temp1.split('-')
temp3 = temp2[0].split()
n1,n2 = float(temp3[0]), float(temp3[1])
temp1 = rp_f.readline()
temp2 = temp1.split('-')
temp3 = temp2[0].split()
Swr,Sgr = float(temp3[0]), float(temp3[1])
rp_f.close()

rel_perms = rel_perm(rel_perm_f)


temp1 = par_f.readline()
temp2 = temp1.split('-')
brine_den = float(temp2[0])
p_write.write('Brine density: %.2f\n' %(brine_den))
temp1 = par_f.readline()
temp2 = temp1.split('-')
co2_den = float(temp2[0])
p_write.write('CO2 density: %.2f\n' %(co2_den))
temp1 = par_f.readline()
temp2 = temp1.split('-')
brine_visc = float(temp2[0])
p_write.write('Brine viscosity: %.2f\n' %(brine_visc))
temp1 = par_f.readline()
temp2 = temp1.split('-')
co2_visc = float(temp2[0])
p_write.write('CO2 viscosity: %.2f\n' %(co2_visc))

# **************** Recurrent data ****************
temp1 = par_f.readline()
temp2 = temp1.split('-')
total_days = int(temp2[0])
p_write.write('Total simulation time: %d days\n' %(total_days))
temp1 = par_f.readline()
temp2 = temp1.split('-')
vol_rate = float(temp2[0])
p_write.write('Volumetric injection rate: %.2f m3/day\n' %(vol_rate))
temp1 = par_f.readline()
temp2 = temp1.split('-')
num_models = int(temp2[0])
p_write.write('Number of models: %d\n' %(num_models))
temp1 = par_f.readline()
temp2 = temp1.split('-')
reporting_intervals = int(temp2[0])
p_write.write('Reporting intervals: %d days\n' %(reporting_intervals))
del temp1, temp2

par_f.close()
p_write.close()

## ************************** End of data read ********************************

# To conserve space, create only one list for perm, porosity etc,
# which can be used on every model while it is being run


# ***************************** Depth file *********************************
data_f = open(pres_file,'r')
depth = [[[0 for i in range(NX)] for j in range(NY)] for k in range(NZ)]
temp1 = data_f.readlines()
data_f.close()
temp_counter = 0
for j in range(NY):
    for k in range(NX):
        depth[1][j][k] = float(temp1[temp_counter].strip())
        for i in range(1,NZ):
            depth[i][j][k]=depth[1][j][k]+(i-1)*dz
        temp_counter = temp_counter+1
del temp_counter, temp1

for models in range(num_models):
    start_time1 = time.clock()
    reports = [[[[0 for i in range(NX)] for j in range(NY)] for k in range(NZ)] for l in range(total_days/reporting_intervals)]
    # Permeability values
    print 'Reading perms from ',perm_file
    data_f = open(perm_file,'r')
    perms = [[[0 for i in range(NX)] for j in range(NY)] for k in range(NZ)]
    for i in range(NZ):
        for j in range(NY):
            for k in range(NX):
                temp1 = data_f.readline().split()
                perms[i][j][k] = float(temp1[models])
    data_f.close()

    # Porosity values
    print 'Reading porosity from ',por_file
    data_f = open(por_file,'r')
    accomodation = [[[0 for i in range(NX)] for j in range(NY)] for k in range(NZ)]
    total_PV = 0
    total_particles = 0
    for i in range(NZ):
        for j in range(NY):
            for k in range(NX):
                temp1 = data_f.readline().split()
                accomodation[i][j][k] = int(20*float(temp1[models]))
                if accomodation[i][j][k]<1:
                    accomodation[i][j][k]=1
                total_particles = total_particles + accomodation[i][j][k]
                total_PV = total_PV + dx*dy*dz*float(temp1[models])
    data_f.close()

    print 'Max accomodation is ',max(max(max(accomodation[1:NZ+1][1:NY+1][1:NX+1])))

    # Volumetric calculations
    N_inj_rate = int(round(vol_rate * total_particles/total_PV))
    print 'Injection rate is ',N_inj_rate

    question = raw_input('Continue?')


    # Initialize particle counts
    carbon_count = [[[0 for i in range(NX)] for j in range(NY)] for k in range(NZ)]

    # Macroscopic velocity (used in fractional flow calculations)
    v_macro_x = vol_rate/(dy*dz)
    v_macro_y = vol_rate/(dx*dz)
    v_macro_z = vol_rate/(dx*dy)

    # Starting walk for current model
    start_time2 = time.clock()
    for curr_time in range(total_days):
        for particle in range(N_inj_rate):
            for inj in range(num_inj):
                carbon_count[inj_z[inj]-1][inj_y[inj]-1][inj_x[inj]-1] = carbon_count[inj_z[inj]-1][inj_y[inj]-1][inj_x[inj]-1] + 1
                x,y,z = inj_x[inj]-1,inj_y[inj]-1,inj_z[inj]-1
                go_back=1
                occupy = [[[0 for i in range(NX)] for j in range(NY)] for k in range(NZ)]
                while (go_back==1):
                    Tr = [0 for i in range(7)]
                    occupy[z][y][x]=1

                    # ********************************** TRANSITION PROBABILITIES *******************************************

##                    if carbon_count[z][y][x]/accomodation[z][y][x]>1-S1r:
##                        Tr_temp=[1000 for i in range(6)]
##                        Tr[1:7]=Tr_temp
##                    else:
                        # positive X-direction
                    if x==NX-1 or occupy[z][y][x+1]==1:
                        Tr[1]=0
                    else:
                        satn = (float(carbon_count[z][y][x])/accomodation[z][y][x]+float(carbon_count[z][y][x+1])/accomodation[z][y][x+1])/2
                        if satn>0.999:
                            satn=0.999
                        krg = rel_perms[1][rel_perms[0].index(round(satn,3))]
                        krw = rel_perms[2][rel_perms[0].index(round(1-satn,3))-1]
                        avg_perm =(perms[z][y][x]*perms[z][y][x+1])**0.5
                        front_sat = front_saturation(rel_perms,co2_den-brine_den,avg_perm,v_macro_x,co2_visc,brine_visc,depth[z][y][x],depth[z][y][x+1],dx)
                        if satn<front_sat:
                            Tr[1] = 0
                        else:
                            velocity = BL_velocity(rel_perms,co2_den-brine_den,avg_perm,v_macro_x,co2_visc,brine_visc,depth[z][y][x],depth[z][y][x+1],dx,satn)
                            Tr[1] = (carbon_count[z][y][x]-carbon_count[z][y][x+1])/accomodation[z][y][x+1] + velocity/v_macro_x

                    # negative X-direction
                    if x==0 or occupy[z][y][x-1]==1:
                        Tr[2]=0
                    else:
                        satn = (float(carbon_count[z][y][x])/accomodation[z][y][x]+float(carbon_count[z][y][x-1])/accomodation[z][y][x-1])/2
                        if satn>0.999:
                            satn=0.999
                        krg = rel_perms[1][rel_perms[0].index(round(satn,3))]
                        krw = rel_perms[2][rel_perms[0].index(round(1-satn,3))-1]
                        avg_perm =(perms[z][y][x]*perms[z][y][x-1])**0.5
                        front_sat = front_saturation(rel_perms,co2_den-brine_den,avg_perm,v_macro_x,co2_visc,brine_visc,depth[z][y][x],depth[z][y][x-1],dx)
                        if satn<front_sat:
                            Tr[2] = 0
                        else:
                            velocity = BL_velocity(rel_perms,co2_den-brine_den,avg_perm,v_macro_x,co2_visc,brine_visc,depth[z][y][x],depth[z][y][x-1],dx,satn)
                            Tr[2] = (carbon_count[z][y][x]-carbon_count[z][y][x-1])/accomodation[z][y][x-1] + velocity/v_macro_x

                    # positive Y-direction
                    if y==NY-1 or occupy[z][y+1][x]==1:
                        Tr[3]=0
                    else:
                        satn = (float(carbon_count[z][y][x])/accomodation[z][y][x]+float(carbon_count[z][y+1][x])/accomodation[z][y+1][x])/2
                        if satn>0.999:
                            satn=0.999
                        krg = rel_perms[1][rel_perms[0].index(round(satn,3))]
                        krw = rel_perms[2][rel_perms[0].index(round(1-satn,3))-1]
                        avg_perm =(perms[z][y][x]*perms[z][y+1][x])**0.5
                        front_sat = front_saturation(rel_perms,co2_den-brine_den,avg_perm,v_macro_x,co2_visc,brine_visc,depth[z][y][x],depth[z][y+1][x],dx)
                        if satn<front_sat:
                            Tr[3] = 0
                        else:
                            velocity = BL_velocity(rel_perms,co2_den-brine_den,avg_perm,v_macro_x,co2_visc,brine_visc,depth[z][y][x],depth[z][y+1][x],dx,satn)
                            Tr[3] = (carbon_count[z][y][x]-carbon_count[z][y+1][x])/accomodation[z][y+1][x] + velocity/v_macro_x

                    # negative Y-direction
                    if y==0 or occupy[z][y-1][x]==1:
                        Tr[4]=0
                    else:
                        satn = (float(carbon_count[z][y][x])/accomodation[z][y][x]+float(carbon_count[z][y-1][x])/accomodation[z][y-1][x])/2
                        if satn>0.999:
                            satn=0.999
                        krg = rel_perms[1][rel_perms[0].index(round(satn,3))]
                        krw = rel_perms[2][rel_perms[0].index(round(1-satn,3))-1]
                        avg_perm =(perms[z][y][x]*perms[z][y-1][x])**0.5
                        front_sat = front_saturation(rel_perms,co2_den-brine_den,avg_perm,v_macro_x,co2_visc,brine_visc,depth[z][y][x],depth[z][y-1][x],dx)
                        if satn<front_sat:
                            Tr[4] = 0
                        else:
                            velocity = BL_velocity(rel_perms,co2_den-brine_den,avg_perm,v_macro_x,co2_visc,brine_visc,depth[z][y][x],depth[z][y-1][x],dx,satn)
                            Tr[4] = (carbon_count[z][y][x]-carbon_count[z][y-1][x])/accomodation[z][y-1][x] + velocity/v_macro_x

                    # positive Z-direction
                    if z==NZ-1 or occupy[z+1][y][x]==1:
                        Tr[5]=0
                    else:
                        satn = (float(carbon_count[z][y][x])/accomodation[z][y][x]+float(carbon_count[z+1][y][x])/accomodation[z+1][y][x])/2
                        if satn>0.999:
                            satn=0.999
                        krg = rel_perms[1][rel_perms[0].index(round(satn,3))]
                        krw = rel_perms[2][rel_perms[0].index(round(1-satn,3))-1]
                        avg_perm =(perms[z][y][x]*perms[z+1][y][x])**0.5
                        front_sat = front_saturation(rel_perms,co2_den-brine_den,avg_perm,v_macro_x,co2_visc,brine_visc,depth[z][y][x],depth[z+1][y][x],dx)
                        if satn<front_sat:
                            Tr[5] = 0
                        else:
                            velocity = BL_velocity(rel_perms,co2_den-brine_den,avg_perm,v_macro_x,co2_visc,brine_visc,depth[z][y][x],depth[z+1][y][x],dx,satn)
                            Tr[5] = (carbon_count[z][y][x]-carbon_count[z+1][y][x])/accomodation[z+1][y][x] + velocity/v_macro_x


                    # negative Z-direction
                    if z==0 or occupy[z-1][y][x]==1:
                        Tr[6]=0
                    else:
                        satn = (float(carbon_count[z][y][x])/accomodation[z][y][x]+float(carbon_count[z-1][y][x])/accomodation[z-1][y][x])/2
                        if satn>0.999:
                            satn=0.999
                        krg = rel_perms[1][rel_perms[0].index(round(satn,3))]
                        krw = rel_perms[2][rel_perms[0].index(round(1-satn,3))-1]
                        avg_perm =(perms[z][y][x]*perms[z-1][y][x])**0.5
                        front_sat = front_saturation(rel_perms,co2_den-brine_den,avg_perm,v_macro_x,co2_visc,brine_visc,depth[z][y][x],depth[z-1][y][x],dx)
                        if satn<front_sat:
                            Tr[6] = 0
                        else:
                            velocity = BL_velocity(rel_perms,co2_den-brine_den,avg_perm,v_macro_x,co2_visc,brine_visc,depth[z][y][x],depth[z-1][y][x],dx,satn)
                            Tr[6] = (carbon_count[z][y][x]-carbon_count[z-1][y][x])/accomodation[z-1][y][x] + velocity/v_macro_x

                    for i in range(7):
                        if Tr[i]<0:
                            Tr[i]=0
                    if sum(Tr[1:7])==0:
                        Tr[0] = 1

                    Tr = [item/float(sum(Tr)) for item in Tr]
                    movement = monte_carlo_sampling(Tr)
                    carbon_count[z][y][x] = carbon_count[z][y][x] - 1
#                    print x,y,z, ' --> ', [round(item,3) for item in Tr], movement,
                    if movement==1:
                        x=x+1
                    if movement==2:
                        x=x-1
                    if movement==3:
                        y=y+1
                    if movement==4:
                        y=y-1
                    if movement==5:
                        z=z+1
                    if movement==6:
                        z=z-1
                    carbon_count[z][y][x] = carbon_count[z][y][x] + 1
#                    print ' -->',x,y,z

                    if movement==0:
                        go_back=0

        if curr_time%reporting_intervals==0:
            write_saturations(carbon_count,accomodation)
        if curr_time%50==0:
            print 'Done with step ',curr_time+1, ' in ',time.clock()-start_time2, ' secs.'
        print 'Done with step ',curr_time+1
    print 'Done with model ',models+1,' in',time.clock()-start_time1, ' secs.'
    report(models+1,NX,NY,NZ)





