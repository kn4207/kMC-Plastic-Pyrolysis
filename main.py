from multiprocessing import Pool
from rdkit import rdBase
from simulation import kmc_sim
from config import TEMPERATURE_VALUES


rdBase.DisableLog('rdApp.*') 

if __name__=='__main__': 
    pool = Pool(processes=8)
    results = pool.map(kmc_sim,TEMPERATURE_VALUES)
    pool.close()
    pool.join()