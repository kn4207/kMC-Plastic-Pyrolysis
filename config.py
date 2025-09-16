TEMPERATURE_VALUES = [693,713,723,733,773,873,973,1073] #Range of temperature values to run kmc_sim on

INIT_CARBON_CHAIN_LENGTH = 100 #Length of initial molecules in the system
INIT_FEEDSTOCK_SIZE = 100 #Initial number of molecules in the system
CHAIN_LENGTH_ENDPOINT = 17.83 #Average chain length of species in the system below which the simulation is terminated
MW_SWEEP_THRESHOLD = 90 #Mol wt. threshold below which species are removed from the system
SWEEP_FREQ = 10 #Frequency (Every x events) at which N2 Sweeping is done to remove species below MW_SWEEP_THRESHOLD
WRITE_FREQ = 1000 #Frequency at which data is written to output
PRINT_FREQ = 500 #Frequency at which event number is printed
RUN_FOLDER = 'kMC run 9 - ave len' ##########Change every run!