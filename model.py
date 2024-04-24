def simulate(all_t_real, all_r_real, all_pd_real, all_pm_real, n):
    
    import random
    import math
    import numpy as np
    
    def find_nearest(list, value): 
        array = np.asarray(list)
        index = (np.abs(array - value)).argmin()
        return index
    
    #starting conditions for random rna amount
    #num = np.random.normal(loc=150, scale=150, size=1) #initial distribution of mRNA per cell, obtained from smFISH experiments
    #abs_num = abs(num)
    #r = int(abs_num) #pick one random number from the distribution
    
    #starting conditions for fixed rna amount
    r = 150 
    pd = 0
    pm = 0

    
    #constants
    ktl = 22889.75 #translation rate
    kd1 = 0.363 #degradation of mRNA
    kd2 = 0.039 #degradation of GFP protein 
    km = 0.49 #maturation of GFP protein 
    
    
    #time settings in hours
    t = 0
    T = 10
    
    #Store data of one simulation
    
    t_sim = []
    r_sim = []
    pd_sim = []
    pm_sim = []
    
    t_sim.append(t)
    r_sim.append(r)
    pd_sim.append(pd)
    pm_sim.append(pm)
    
    while t < T:
        K_1 = r * kd1 #rna is degraded
        K_2 = r * ktl #rna is translated 
        K_3 = pd * km #protein is matured
        K_4 = pm * kd2 #protein is degraded
        
        #SUM of K
        K_SUM = K_1 + K_2 + K_3 + K_4
        
        #time step
        tau_r = random.random()
        tau = -(1/K_SUM)* math.log(tau_r)
        t = t + tau
        
        # Probability of each reaction
        p_1 = K_1 / K_SUM
        p_2 = K_2 / K_SUM
        p_3 = K_3 / K_SUM
        p_4 = K_4 / K_SUM

        
        #choosing a reaction and effects of the reaction
        
        rnd = random.random()
        
        if rnd < p_1:
            r = r - 1
        elif rnd < p_2 + p_1:
            pd = pd + 1
        elif rnd < p_3 + p_2 + p_1:
            pm = pm + 1
            pd = pd - 1
        elif rnd <= p_4 + p_3 + p_2 + p_1:
            pm = pm - 1
        
        #Store values and start again
        
        t_sim.append(t)
        r_sim.append(r)
        pd_sim.append(pd)
        pm_sim.append(pm)

        
    #from lists to arrays
    
    t_sima = np.array(t_sim)
    r_sima = np.array(r_sim)
    pd_sima = np.array(pd_sim)
    pm_sima = np.array(pm_sim)
    
    #Regenerate time-representative (real) lists of each specie
    
    t_real = []
    r_real = []
    pd_real = []
    pm_real = []
    
    t = 0
    T = 10
    
    while(t < T):
        t_real.append(t)
        index = find_nearest(t_sima, t)
        r_real.append(r_sima[index])
        pd_real.append(pd_sima[index])
        pm_real.append(pm_sima[index])
        t = t + 15/60 #every 15 minutes of simulated expression
    
    #Append the list of each simulation in the bigger list of all the simulations (defined in the execute script)
    all_t_real.append(t_real)
    all_r_real.append(r_real)
    all_pd_real.append(pd_real)
    all_pm_real.append(pm_real)