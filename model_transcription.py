def simulate(all_t_real, all_OFF_real, all_ON_real, all_r_real, all_pd_real, all_pm_real, n):
    
    import random
    import math
    import numpy as np
    
    def find_nearest(list, value): #This function allows for aproximations and return the index (not the number!)
        array = np.asarray(list)
        index = (np.abs(array - value)).argmin()
        return index
    
    #starting conditions
    OFF = 1 #promoter starts in OFF state
    ON = 0
    r = 0
    pd = 0
    pm = 0

    
    #constants
    kon = 0.3 #promoter in the ON state
    koff = 1 #promoter in the OFF state
    kt = 10 #transcription rate
    ktl = 4000 #translation rate
    kd1 = 0.15 #degradation of mRNA
    kd2 = 1 #degradation of GFP protein 
    km = 2 #maturation of GFP protein 
    
    
    #time settings in hours
    t = 0
    T = 20
    
    #Store data of one simulation
    
    t_sim = []
    OFF_sim = []
    ON_sim = []
    r_sim = []
    pd_sim = []
    pm_sim = []
    
    t_sim.append(t)
    OFF_sim.append(OFF)
    ON_sim.append(ON)
    r_sim.append(r)
    pd_sim.append(pd)
    pm_sim.append(pm)
    
    while t < T:
        K_1 = OFF * kon
        K_2 = ON * kt
        K_3 = ON * koff
        K_4 = r * kd1 #rna is degraded
        K_5 = r * ktl #rna is translated 
        K_6 = pd * km #protein is matured
        K_7 = pm * kd2 #protein is degraded
        
        #SUM of K
        K_SUM = K_1 + K_2 + K_3 + K_4 + K_5 + K_6 + K_7
        
        #time step
        tau_r = random.random()
        tau = -(1/K_SUM)* math.log(tau_r)
        t = t + tau
        
        # Probability of each reaction
        p_1 = K_1 / K_SUM
        p_2 = K_2 / K_SUM
        p_3 = K_3 / K_SUM
        p_4 = K_4 / K_SUM
        p_5 = K_5 / K_SUM
        p_6 = K_6 / K_SUM
        p_7 = K_7 / K_SUM

        
        #choosing a reaction and effects of the reaction
        
        rnd = random.random()
        
        if rnd < p_1:
            OFF = 0
            ON = 1
        elif rnd < p_2 + p_1:
            r = r + 1
        elif rnd < p_3 + p_2 + p_1:
            OFF = 1
            ON = 0
        elif rnd < p_4 + p_3 + p_2 + p_1:
            r = r - 1
        elif rnd < p_5 + p_4 + p_3 + p_2 + p_1:   
            pd = pd + 1
        elif rnd < p_6 + p_5 + p_4 + p_3 + p_2 + p_1:
            pm = pm + 1
            pd = pd - 1
        elif rnd < p_7 + p_6 + p_5 + p_4 + p_3 + p_2 + p_1:
            pm = pm - 1
        
        #Store values and start again
        
        t_sim.append(t)
        OFF_sim.append(OFF)
        ON_sim.append(ON)
        r_sim.append(r)
        pd_sim.append(pd)
        pm_sim.append(pm)

        
    #from lists to arrays
    
    t_sima = np.array(t_sim)
    OFF_sima = np.array(OFF_sim)
    ON_sima = np.array(ON_sim)
    r_sima = np.array(r_sim)
    pd_sima = np.array(pd_sim)
    pm_sima = np.array(pm_sim)
    
    #Regenerate time-representative (real) lists of each specie
    
    t_real = []
    OFF_real = []
    ON_real = []
    r_real = []
    pd_real = []
    pm_real = []
    
    t = 0
    T = 20
    
    while(t < T):
        t_real.append(t)
        index = find_nearest(t_sima, t)
        OFF_real.append(OFF_sima[index])
        ON_real.append(ON_sima[index])
        r_real.append(r_sima[index])
        pd_real.append(pd_sima[index])
        pm_real.append(pm_sima[index])
        t = t + 15/60
    
    #Append the list of each simulation in the bigger list of all the simulations (defined in the execute script)
    all_t_real.append(t_real)
    all_OFF_real.append(OFF_real)
    all_ON_real.append(ON_real)
    all_r_real.append(r_real)
    all_pd_real.append(pd_real)
    all_pm_real.append(pm_real)