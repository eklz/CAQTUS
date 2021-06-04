# %%
from astropy.io import fits

from .Cn2 import *

from numpy import *
from tqdm import tqdm
import pandas as pd
from multiprocessing import Pool
import os
from os import truncate, remove
from shutil import copyfile
from time import time, time_ns
import re
import logging
import glob
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
    profil_Cn2 = array([alt[:-1], Cn2[:-1]*dz])
    profil_ws = array([alt[:-1], Ws[:-1]])
    hdu_Cn2 = fits.PrimaryHDU(profil_Cn2)
    hdu_Cn2 = fits.HDUList([hdu_Cn2])
    hdu_Cn2.writeto(path + name_cn2 + '.fits', overwrite=1)

    hdu_ws = fits.PrimaryHDU(profil_ws)
    hdu_ws = fits.HDUList([hdu_ws])
    hdu_ws.writeto(path + name_ws + '.fits', overwrite=1)

# %%


def coupled_flux(data, date: str, params={'fech': 4.7e3, 'ttot': 1, 'delay': 2, 'dpup': 0.6, 'nmax': 12.0, 'nrad': 30.0, seed: 10, 'delay': 2},  
                 path='/scratchm/eklotz/calcul_coupled_flux/Params/tmp/', default='/scratchm/eklotz/calcul_coupled_flux/Params/param_geo_oa_60cm.txt', logging = 'error.log'):
    
    logging.basicConfig(filename='error.log')
    ndate = date.replace(" ", "_")
    ndate = ndate.replace("-", "_")
    ndate = ndate.strip(':00:00')
    name_cn2 = ndate + '_Cn2'
    name_ws = ndate + '_wspeed'
    write_params(ndate, name_cn2, name_ws, path, default)
    tmp = data[date]
    Cn2 = tmp.Cn2.values
    alt = tmp.alt.values
    WS = tmp.wspeed.values
    make_fits(name_cn2, name_ws, Cn2, WS, alt, path)

    from idlpy import IDL
    IDL.run('.COMPILE -v /scratchm/eklotz/calcul_coupled_flux/SAOST_EK/open_param_turandot_ek.pro')
    IDL.run(
        '.COMPILE -v /scratchm/eklotz/calcul_coupled_flux/SAOST_EK/gen_profils_leo_ek.pro')
    IDL.run('.COMPILE -v /scratchm/eklotz/calcul_coupled_flux/SAOST_EK/simu_oa_simplifiee_ek.pro')
    IDL.dpup = 0.6
    IDL.nmax = 12.0
    IDL.run('tnmax = indgen(nmax)+1.')
    IDL.nrad = 30.0
    IDL.run('nmodes = round((nrad*(nrad+3))/2+1)')
    os.mkdir(path + str(i)+'/')
    IDL.dir = path + str(i)+'/'
    IDL.file_param = '/scratchm/eklotz/calcul_coupled_flux/Params/tmp/' + ndate + '.txt'
    cmd = 'results = simu_oa_simplifiee_ab2(file_param = file_param, nmax = nmax, fech = 4.7e3,  nmodes = nmodes, nocc = 4700, ttot = 1., delay = 2, seed = seed, dir = dir, fast = 0, tempo = 1)'
    try:
        logging.info(f'{date} done')
        IDL.run(cmd)

        res = IDL.results
        df = pd.DataFrame(res['FC'], columns=['FC'])
        df['date'] = date
        logging.info(f'Ok pour {date}')
        return df

    except:
        logging.error(f' Erreur rencontrée pour {date}')
