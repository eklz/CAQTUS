# %%
import numpy as np
import pandas as pd
# %%
from collections.abc import Iterable
#import seaborn as sns
from collections import namedtuple
from scipy.integrate._ivp.radau import T
from tqdm import tqdm
import warnings
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from .Dewan import *
from .Maths import *
from .Masciadri import *
from .decimation import *

import swifter 


class Cn2:

    warnings.simplefilter("ignore")
    tqdm.pandas()

    def __init__(self, dataFrame, date="date", alt="alt", Cn2="Cn2", wspeed="wspeed", **kwargs):
        """Initialize Cn2 class

        Args:
            dataFrame (pd.DataFrame): DataFrame containing the profiles
            date (str, optional): [description]. Defaults to "date".
            alt (str, optional): [description]. Defaults to "alt".
            Cn2 (str, optional): [description]. Defaults to "Cn2".
            wspeed (str, optional): [description]. Defaults to "wspeed".
        """
        ColTuple = namedtuple(
            "Columns", ["date", "alt", "Cn2", "wspeed", *kwargs.keys()])
        self.columns = ColTuple(date, alt, Cn2, wspeed, **kwargs)
        self._data = dataFrame.sort_values(by=[self.columns.date, self.columns.alt])
        self.dates = np.unique(self._data[self.columns.date])


    def filtre(self, size=11, type=1, inplace=True, column = 'Cn2'):
        """Filtre the given profile

        Args:
            size (int, optional): Size of the filtre to apply. Defaults to 11.
            type (int, optional): 1 --> moving_average, 2--> gaussianfilter1D. Defaults to 1.
            inplace (bool, optional): Defaults to True.
            column (str, optional): Column to filter. Defaults to 'Cn2'.
        """
        if type == 1:
            filtre = moving_average
        elif type == 2:
            filtre = gaussian_filter1d
        else:
            raise ValueError("Wrong value for type, please enter 1 for moving average and 2 for gaussian")

        print("Applying filter : ")
        col = getattr(self.columns, column)
        
        if len(self.dates) == 1 : 
            _filtered = self._data.groupby(self.columns.date, sort = False, dropna=False).progress_apply(
                        lambda x: pd.Series(filtre(x[col], size))).values[0]
        else :
            _filtered = self._data.groupby(self.columns.date, sort = False, dropna=False).progress_apply(
                        lambda x: pd.Series(filtre(x[col], size))).values
        if inplace:
            self._data[col] = _filtered
        else:
            copy = self._data.copy()
            copy[col] = _filtered
            return Cn2(copy, **self.columns._asdict())

    def rm_incomplete(self, alt_min, verb = 1, inplace = True):
        """Remove all profiles that don't have values after at least alt_min

        Args:
            alt_min (int)
        """

        invalid = abs((self._data.groupby('date').apply(lambda x : np.max(x.alt.values))>=alt_min)-1)
        invalid_date = invalid[invalid == 1].index.values
        res = self._data.loc[~self._data.date.isin(invalid_date)]
        
        if verb : 
            print(f'{len(invalid_date)} dates invalides sur {len(self.dates)}')
            print(invalid_date)
        
        if inplace:
            self._data = res
            self.dates = np.unique(self._data[self.columns.date])
        else:
            return Cn2(res, **self.columns._asdict())
               
    def decimate(self, nbsegments, nbpoints_per_segment, log = True):
        
        if nbsegments >1 : 
            def foo(x):
                try : 
                    roots, res = interp_gauss_legendre_pp(nbpoints_per_segment, np.array_split(x[self.columns.alt].values, nbsegments), np.array_split(x[self.columns.Cn2].values, nbsegments), norm = True)
                    f = interp1d(x.alt.values, x.wspeed.values)
                    return pd.DataFrame(np.transpose([roots, res, f(roots)]), columns=[self.columns.alt, self.columns.Cn2, self.columns.wspeed])
                except : 
                    print(f'Error for {x}')
        else :
            def foo(x):
                try : 
                    roots, res =interp_gauss_legendre(nbpoints_per_segment, x[self.columns.alt].values, x[self.columns.Cn2].values, norm=True)
                    f = interp1d(x.alt.values, x.wspeed.values)
                    return pd.DataFrame(np.transpose([roots, res, f(roots)]), columns=[self.columns.alt, self.columns.Cn2, self.columns.wspeed])
                except : 
                    print(f'Error for {x}')
                           
        df = self._data.groupby(self.columns.date).progress_apply(lambda x : foo(x))
        df = df.reset_index()
        df = df.drop(columns=['level_1'])
        return Cn2(df,date="date", alt="alt", Cn2="Cn2", wspeed="wspeed")


    def moments(self, lambda_m = 1.55e-6, zenithAngle = 0):
        f = (lambda x : calc_moments(x[self.columns.Cn2].values, x[self.columns.alt].values, x[self.columns.wspeed].values,
                                    lambda_m = lambda_m, zenithAngle = zenithAngle))
        res = self._data.groupby(self.columns.date, sort = False, dropna=False).apply(lambda x : pd.Series(f(x)))
        return res

    
    def set_ground_level(self, value = None, inplace =True):
        if not value:
            f = lambda x : pd.Series(x[self.columns.alt].values-x[self.columns.alt].values[0])
            _alt = self._data.groupby(self.columns.date, sort = False, dropna=False).apply(lambda x : f(x)).values
        else : 
            f = lambda x : pd.Series(x[self.columns.alt].values-value)
            _alt = self._data.groupby(self.columns.date, sort = False, dropna=False).apply(lambda x : f(x)).values
        if len(self.dates)==1:
            _alt = _alt[0]
                
        if inplace :
            self._data[self.columns.alt] = _alt
        else :
            copy = self._data.copy()
            copy[self.columns.alt] = _alt
            return Cn2(copy, **self.columns._asdict())

                        
    def rm_zeros(self, inplace = True):
        
        nnzeros = (self._data[self.columns.Cn2] != 0) & (~np.isnan((self._data[self.columns.Cn2])))
        if inplace:
            self._data = self._data[nnzeros]
        else : 
            return Cn2(self._data[nnzeros], **self.columns._asdict())

        
    def change_elevation(self, theta, inplace = True):
        """Changer l'élévation du profil en modifiant l'écart entre les valeurs de alt

        Args:
            theta ([int]): angle d'élévation en radians
        """            
        
        if inplace : 
            self._data.alt = self._data.alt.apply(lambda x : x/np.sin(theta))
        else : 
            copy = self._data.copy()
            copy.alt = copy.alt.apply(lambda x : x/np.sin(theta))
            return Cn2(copy, **self.columns._asdict())

        
    def to_csv(self, path):
        """Write the Cn2 object to a CSV file

        Args:
            path (str): Path of the CSV file to write, if file already exist it is being overwrite. 
        """
        self._data.to_csv(path, index = False)
        
    @classmethod
    def read_csv(cls, path, date="date", alt="alt", Cn2="Cn2", wspeed="wspeed", **kwargs):
        """Return an object of class Cn2 from a path to a CSV file

        Args:
            path (str): Path to the CSV file
            nameDate (str, optional): [description].Name of the column containing the dates. Defaults to "date".
            nameAlt (str, optional): [description]. Name of the column containing the altitutes. Defaults to "alt".
            nameCn2 (str, optional): [description]. Name of the column containing the Cn2 profiles. Defaults to "Cn2".
            nameWindSpeed (str, optional): [description]. Name of the column containing the wind speed profiles. Defaults to "wspeed".

        Returns:
            <Cn2>
        """
        return cls(pd.read_csv(path).reset_index(drop=True), date= date, alt= alt, Cn2= Cn2, wspeed= wspeed, **kwargs)

    @classmethod
    def from_soundings(cls, soundings, model="Dewan", alt="alt", date="date", press="press",
                       temp="temp", meridional="u", zonal="v", wspeed="wspeed", shear=None, gamma=None, **kwargs):
        """Generate Cn2 from radiosounding data.
        Args:
            soundings (str || pandas;DataFrame): Either a path to a CSV file (str) or a Pandas DataFrame.
            nbCPU (int, optional): Number of CPU used for the calcul of Cn2. Defaults to 1.
            model (str, optional): Model used to calculate Cn2 in ["Dewan", "Masciadri"]. Defaults to "Dewan". 
            nameAlt (str, optional): Name of the column containing the altitudes. Defaults to "alt".  
            nameDate (str, optional): Name of the column containing the dates. If there is only one measurement set nameDate = None. Defaults to "date". 
            namePress (str, optional): Name of the column containing the pressure profiles. Defaults to "press".
            nameTemp (str, optional): Name of the column containing the temperature profiles. Defaults to "temp". 
            nameU (str, optional): Name of the column containing the meridional wind profiles. Defaults to "u". 
            nameV (str, optional): Name of the column containing the zonal wind profiles. Defaults to "v". 
            nameShear (str, optional): Name of the column containing the wind shear profiles. If None it is being recalculated. Defaults to None. 
            gamma (float, optional): gamma used in the calcul of the Cn2 (Default 2.8 for Dewan and 1.5 for Masciadri). Default to None.
        """

        _gamma = {"Dewan": 2.8, "Masciadri": 1.5}

        ColTuple = namedtuple("Columns", [
                              "date", "alt", "press", "temp", "meridional", "zonal", "wspeed", "shear", *kwargs.keys()])
        columns = ColTuple(date, alt, press, temp, meridional,
                           zonal, wspeed, shear, **kwargs)
        # Check the variable soundings
        if isinstance(soundings, str):
            data = pd.read_csv(soundings).reset_index(drop=True)
            data.sort_values(by=[columns.date, columns.alt], inplace = True)
        elif isinstance(soundings, pd.DataFrame):
            data = soundings.sort_values(by=[columns.date, columns.alt])
        else:
            raise ValueError(
                "Invalid argument, enter either a path to a CSV file (str) or a Pandas DataFrame.")

        # Check the variable model
        if model not in _gamma.keys():
            raise ValueError(
                "Invalid model, select a model in [\"Dewan\", \"Masciadri\"]")

        # Check the variable gamma
        if not gamma:
            gamma = _gamma[model]

        # #Check the variable date
        # if not columns.date:
        #     dates = [1]
        # else :
        #     dates = np.unique(data[columns.date])

        # Check the wind shear
        if not columns.shear:  # If no shear we calculate it
            print("Calcul of the wind shear for the sounding data")
            data["S"] = data.groupby(columns.date, sort = False, dropna=False).progress_apply(lambda x: pd.Series(
                grad_wind(x[columns.meridional], x[columns.zonal], x[columns.alt]))).values
            columns = columns._replace(shear="S")

        if model == "Dewan":
            print("Calcul of the Cn2 for the sounding data")
            data['Cn2'] = data.groupby(columns.date, sort=False, dropna=False).progress_apply(lambda x: pd.Series(Cn2_Dewan(x[columns.press].values, x[columns.temp].values,
                                                                                                  x[columns.alt].values, x[columns.shear].values, gamma=gamma))).values

        # if model == "Masciadri":

        #     print("Calcul of the Cn2 for the sounding data")
        #     data['Cn2'] = data.groupby(columns.date).progress_apply(lambda x: pd.Series(Cn2_Dewan(x[columns.press].values, x[columns.temp].values,
        #                                                                                           x[columns.alt].values, x[columns.shear].values, gamma=gamma))).values

        return cls(data, **columns._asdict())


    def __getattr__(self, name: str):  # Attribute can be any column name
        try:
            return self._data[name]
        except :
            try:
                return self._data[getattr(self.columns, name)]
            except:
                raise AttributeError


    def __getitem__(self, name):

        if isinstance(name, int):  # If name is int we return the corresponding day
            try:
                return Cn2(self._data[self._data[self.columns.date] == self.dates[name]].reset_index(drop=True), **self.columns._asdict())
            except:
                raise IndexError(f'Index out of range {len(self.dates)}')

        # If name is a slice (ex : [0:100:2]) we make sure to return the right thing
        if isinstance(name, slice):
            try:
                return Cn2(self._data[self._data[self.columns.date].isin(self.dates[name])], **self.columns._asdict())
            except:
                raise IndexError(f'Index out of range {len(self.dates)}')
            
        if name in self.dates:  # If a date is given then we return the date
            return Cn2(self._data[self._data[self.columns.date] == name], **self.columns._asdict())

        try:  # If a column name is given, either the true name or the corresponding "key" in columns
            return self._data[name]
        except AttributeError:
            try:
                return self._data[getattr(self.columns, name)]
            except AttributeError:
                raise AttributeError
        
        if isinstance(name, Iterable):
            try:
                return Cn2(self._data[self._data[self.columns.date].isin(self.dates[name])], **self.columns._asdict())
            except:
                raise IndexError(f'Index out of range {len(self.dates)}')

        
        pass

    def __add__(self, other):
        data = pd.concat([self._data, other._data], axis = 0)
        return Cn2(data, **self.columns._asdict())

    def __len__(self):  # return the number of profiles (len(dates))
        return len(self.dates)


    def __repr__(self):
        return repr(self._data)


if __name__ == '__main__':
    alpha = Cn2.read_csv('Cn2_Tenerife_2020_fromSoundings.csv')

