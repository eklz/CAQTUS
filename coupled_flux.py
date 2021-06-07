# %%
from astropy.io import fits

from .Cn2 import Cn2
from tqdm import tqdm
import pandas as pd
from multiprocessing import Pool
import os
import shutil
from os import truncate, remove
from shutil import copyfile
from time import time, time_ns
import re
import logging
import numpy as np

# %%


def write_params(name_txt, name_cn2, name_ws,  path='/scratchm/eklotz/calcul_coupled_flux/Params/tmp/', default='/scratchm/eklotz/calcul_coupled_flux/Params/param_geo_oa_60cm.txt'):

    dest = path + name_txt + '.txt'

    copyfile(default, dest)

    with open(dest, 'r+') as f:
        text = f.read()
        path_to_Cn2 = path + name_cn2 + '.fits'
        path_to_ws = path + name_ws + '.fits'
        text = re.sub('/path/to/Cn2', path_to_Cn2, text)
        text = re.sub('/path/to/ws', path_to_ws, text)
        f.seek(0)
        f.write(text)
        f.truncate()


def make_fits(name_cn2, name_ws, Cn2, Ws, alt, path='/scratchm/eklotz/calcul_coupled_flux/Params/tmp/'):
    dz = alt[1:]-alt[:-1]
    profil_Cn2 = np.array([alt[:-1], Cn2[:-1]*dz])
    profil_ws = np.array([alt[:-1], Ws[:-1]])
    hdu_Cn2 = fits.PrimaryHDU(profil_Cn2)
    hdu_Cn2 = fits.HDUList([hdu_Cn2])
    hdu_Cn2.writeto(path + name_cn2 + '.fits', overwrite=1)

    hdu_ws = fits.PrimaryHDU(profil_ws)
    hdu_ws = fits.HDUList([hdu_ws])
    hdu_ws.writeto(path + name_ws + '.fits', overwrite=1)

# %%


def coupled_flux(data, date: str, params={'fech': 4.7e3, 'ttot': 1, 'delay': 2, 'dpup': 0.6, 'nmax': 12.0, 'nrad': 30.0, 'seed': 10, 'delay': 2},  
                 path='/scratchm/eklotz/calcul_coupled_flux/Params/tmp/', default='/scratchm/eklotz/calcul_coupled_flux/Params/param_geo_oa_60cm.txt', logFile = 'error.log'):
    
    logging.basicConfig(filename=logFile)
      
    ndate = date.replace(" ", "_")
    ndate = ndate.replace("-", "_")
    ndate = ndate.strip(':00:00')
    name_cn2 = ndate + '_Cn2'
    name_ws = ndate + '_wspeed'
    
        
    dir = path +ndate +'/'
    os.mkdir(dir)
    
    write_params(ndate, name_cn2, name_ws, path = dir, default= default)
    tmp = data[date]
    Cn2 = tmp.Cn2.values
    alt = tmp.alt.values
    WS = tmp.wspeed.values
    make_fits(name_cn2, name_ws, Cn2, WS, alt, path = dir)

    from idlpy import IDL
    IDL.run('.COMPILE -v /scratchm/eklotz/calcul_coupled_flux/SAOST_EK/open_param_turandot_ek.pro')
    IDL.run(
        '.COMPILE -v /scratchm/eklotz/calcul_coupled_flux/SAOST_EK/gen_profils_leo_ek.pro')
    IDL.run('.COMPILE -v /scratchm/eklotz/calcul_coupled_flux/SAOST_EK/simu_oa_simplifiee_ek.pro')
    
    for key in params : 
        IDL.run(f'{key} = {params[key]} ')

    IDL.run('tnmax = indgen(nmax)+1.')
    IDL.run('nmodes = round((nrad*(nrad+3))/2+1)')
    IDL.run('nocc = fech*ttot')
    IDL.dir = path +ndate +'/'
    IDL.file_param = dir + ndate + '.txt'
    cmd = 'results = simu_oa_simplifiee_ab2(file_param = file_param, nmax = nmax, fech = fech,  nmodes = nmodes, nocc = nocc, ttot = ttot, delay = delay, seed = seed, dir = dir, fast = 0, tempo = 1)'
    try:
        logging.info(f'{date} done')
        IDL.run(cmd)

        res = IDL.results
        df = pd.DataFrame(res['FC'], columns=['FC'])
        df['date'] = date
        logging.info(f'Ok pour {date}')
        #shutil.rmtree(dir)
        return df

    except:
        #shutil.rmtree(dir)
        print(f' Erreur rencontrée pour {date}')
        logging.error(f' Erreur rencontrée pour {date}')

def get_all_coupled_flux(data, nbCores, params={'fech': 4.7e3, 'ttot': 1, 'delay': 2, 'dpup': 0.6, 'nmax': 12.0, 'nrad': 30.0, 'seed': 10, 'delay': 2},  
                 path='/scratchm/eklotz/calcul_coupled_flux/Params/tmp/', default='/scratchm/eklotz/calcul_coupled_flux/Params/param_geo_oa_60cm.txt', logFile = 'error.log'):
    
    def foo(i): return coupled_flux(data, i, params=params, path=path, default=default, logFile = logFile)
    
    dates = data.dates
    print('0K')
    with Pool(nbCores) as P : 
        res = P.map(foo, list(dates))
        
    res = pd.concat(res)
    
    return res
    
# %%
