import numpy as np
from numpy.core.numeric import NaN
from scipy.io import loadmat as load

def calc_power_margin(ROP_dw_dBm_init, time, Req_MI, InterleavingTime = 50*1e-3, step_power_margin_dB = 0.1):
    
    tmp = load('/scratchm/eklotz/Lib/CAQTUS/perfo/IM_10Gbps_DPSK_EDFA_PIN_SOFT.mat')          
    # donne la correspondance entre la puissance injectée en entrée du LNOA et l'information mutuelle par paquet FEC
    ROP_MI = tmp['ROP_MI'].T[0]
    meanIM = tmp['meanIM'].T[0]
    power_margin_dB = 0
    
    while 1 : 
        ROP_dw_dBm = ROP_dw_dBm_init - power_margin_dB
        
        MI = np.zeros(len(ROP_dw_dBm))
        storageIndexROP  = np.zeros(len(ROP_dw_dBm))
        
        for idx_time in range(len(time)):
            tmp = np.abs(ROP_MI - ROP_dw_dBm[idx_time])
            indexROP = np.argmin(tmp)
            
            # En suivant, le but est de choisir la ROP la plus proche mais la plus
            # pessimiste.
            
            # On identifie les ROP les plus proches mais qui sont plus optimistes
            if (ROP_dw_dBm[idx_time] < ROP_MI[indexROP]):
            # On les transforme en plus pessimiste
                indexROP = indexROP + 1
            
            # Pour les indices qui depassent la taille max des ROP disponibles, on leur
            # affecte la derniere ROP existante.
            if (indexROP > (len(ROP_MI)-1)) : indexROP = len(ROP_MI)-1
            
            MI[idx_time] = meanIM[indexROP]
            
            storageIndexROP[idx_time] = indexROP
        
        
        ### application de l'entrelaceur
        time_step = time[1] - time[0]
        Nb_sample_in_interleaver_windows = int(round(InterleavingTime/time_step))
        Nb_sliding_windows = len(time) - Nb_sample_in_interleaver_windows
        meanIM_interleaved = np.zeros(Nb_sliding_windows)
        
        for index in range(Nb_sliding_windows) : # slidding windows
            MI_tmp = MI[index : index+Nb_sample_in_interleaver_windows+1]
            meanIM_interleaved[index]= sum(MI_tmp)/len(MI_tmp)

        
        ### calcul du pire cas de MI sur la série temporelle et comparaison au minimum requis par le modem de reception
        min_MI = np.min(meanIM_interleaved)
        
        if Req_MI> min_MI : # si le pire cas de MI de la série temporelle inferieur au minimum requis par le modem de reception: le lien n'est plus sans erreur. On récupère la marge système précédente pour laquelle le lien était encore sans erreur.
            if power_margin_dB >0: 
                return power_margin_dB - step_power_margin_dB
            else : return -1
            
        power_margin_dB = power_margin_dB + step_power_margin_dB



def get_limit_power(TC, time, Req_MI, init = (-36.3, -80,4), InterleavingTime = 50*1e-3, step_power_margin_dB = 0.1) :
    TC_dB = 10*np.log10(TC)
    sup = init[0]
    inf = init[1]
    
    if calc_power_margin(TC_dB + sup, time, Req_MI, InterleavingTime , step_power_margin_dB) < 0 : 
        return NaN
    
    while 1 : 
        dic = (inf + sup)/2
        power_margin = calc_power_margin(TC_dB + dic, time, Req_MI, InterleavingTime , step_power_margin_dB)
        
        if power_margin >= 0 : sup = dic
        else : inf = dic
        
        if sup - inf <= step_power_margin_dB : 
            return sup
        
    
    