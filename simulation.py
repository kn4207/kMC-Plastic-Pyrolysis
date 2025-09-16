import os
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt #type:ignore
import time
import math
import numpy as np
import pandas as pd
import json
from utils import atom_metadata, n2_sweep, update_population
from tag import tag_polymer,tag_polymer_radical,tag_radical
from reaction_lib import ReactionLib
from config import WRITE_FREQ,SWEEP_FREQ,PRINT_FREQ,INIT_CARBON_CHAIN_LENGTH,INIT_FEEDSTOCK_SIZE,CHAIN_LENGTH_ENDPOINT,RUN_FOLDER


def update_species(species_df,species_pop,species): #Checks if species is new and adds to species_df. Also enumerates current population of all species and Total Reaction sites for each reaction
    s = Chem.MolFromSmiles(species)
    rad_list = ['[H][H]','[H]','[CH3]']
    rxn_list = ["rxn{}".format(i) for i in range(1, 64)]

    start_time = time.time()
    new_adata_time = start_time
    tag_time = start_time
    append_time = start_time
    
    existing_entry = species_df[species_df['Species SMILES'] == species]
    find_time = time.time()
    if existing_entry.empty:
        adata = atom_metadata(species)
        rad_bool = bool(adata['Rad'] or species in rad_list)
        sat_bool = bool(adata['Sp2'])
        chain_len = species.count('C')
        new_adata_time = time.time()
        if str(species)=='[H]' or str(species)=='[H][H]' or str(species)=='[CH3]':
            r_tags = tag_radical(species)
        elif not adata['Rad']:
            r_tags = tag_polymer(species)
        else:
            r_tags = tag_polymer_radical(species)
        tag_time = time.time()
        new_species = {
                    'Species SMILES':species,
                    'Current Population':1,
                    'Radical?': rad_bool,
                    'Double Bond?': sat_bool,
                    'Chain Length': chain_len,
                    'Molecular Weight':MolWt(s)
                    }
        new_species.update(r_tags)
        
        species_df = species_df._append(new_species, ignore_index=True)
        append_time = time.time()
    for species_smiles, population_history in species_pop.items():
        matching_rows = species_df[species_df['Species SMILES'] == species_smiles]
        
        if not matching_rows.empty:
            idx = matching_rows.index[0]
            current_population = population_history[-1]
            species_df.loc[idx, 'Current Population'] = current_population

    pop_time = time.time()

    for rxn in rxn_list:
        site_density = rxn+'-Reaction Site Density'
        total_sites = rxn+'-Total Reaction Sites'
        species_df[total_sites] = species_df['Current Population']*species_df[site_density]

    return(species_df)

def initialize_species(all_species_pop,species_pop,species_df,initial_species): #Pass a dictionary containing SMILES of initial species present at t=0 along with their populations to populate species_pop and species_df
    for species, population in initial_species.items():
        all_species_pop[species] = [population]
        species_pop[species] = [population] #Update Population first
        species_df = update_species(species_df,species_pop,species)
    return(all_species_pop,species_df,species_pop)

def calc_reaction_rate(species_df,reaction_rates,reaction_lib,T,event_index): #Calculates reaction rate constant and reaction rate for each reaction depending on species population at a certain event
    k_boltz = 3.297*10**-27
    R_gas = 1.985*10**-3
    h_planck = 1.58365200764818*10**-37

    #Initialize the reaction_rate dictionary
    reaction_rates[event_index] = {}
    h_transfer = ["rxn{}".format(i) for i in range(29, 42)]
    reaction_rate_constants = {}
    reaction_rate_constants[event_index] = {}

    #Loop through each reaction
    for reaction_name, reaction_data in reaction_lib.kinetic_params.items():
        enthalpy = reaction_data['act_enthalpy']
        entropy = reaction_data['act_entropy']
        num_reactants = reaction_data['reactants']

        #Calculate rate constant for specified temperature
        tst_prefactor = ((k_boltz*T/h_planck)*math.exp(entropy/R_gas))
        activation_energy = (math.exp((-1*enthalpy)/(R_gas*T)))
        rate_constant = tst_prefactor*activation_energy

        if num_reactants == 2:#Check Order of Reaction
            reactant1_col = species_df[reaction_name+'-Reactant 1']
            reactant2_col = species_df[reaction_name+'-Reactant 2']

            mask1 = reactant1_col ^ reactant2_col #XOR logic
            mask2 = reactant1_col & reactant2_col #AND logic to find species that can participate as both reactants 1 and 2

            r1_reaction_sites_partsum1 = (species_df[mask1][reaction_name+'-Reactant 1']*species_df[mask1][reaction_name+'-Total Reaction Sites']).sum()
            r2_reaction_sites_partsum1 = (species_df[mask1][reaction_name+'-Reactant 2']*species_df[mask1][reaction_name+'-Total Reaction Sites']).sum()
            r1_reaction_sites_partsum2 = (species_df[mask2][reaction_name+'-Reactant 1']*species_df[mask2][reaction_name+'-Total Reaction Sites']).sum()
            r2_reaction_sites_partsum2 = (species_df[mask2][reaction_name+'-Reactant 1']*species_df[mask2][reaction_name+'-Reaction Site Density']*(species_df[mask2]['Current Population']-1)).sum() 
            #Ensures species that can participate as both reactants have a population of >=2 in order to contribute to reaction rate

            r1_total = r1_reaction_sites_partsum1+r1_reaction_sites_partsum2
            r2_total = r2_reaction_sites_partsum1+r2_reaction_sites_partsum2

            reaction_rate = rate_constant*r1_total*r2_total

        else:
            reactant1_col = species_df[reaction_name+'-Reactant 1']
            reaction_sites_col = species_df[reaction_name+'-Total Reaction Sites']
            reaction_sites_sum = (reaction_sites_col*reactant1_col).sum()

            reaction_rate = rate_constant*reaction_sites_sum

        if reaction_rate<0:
            reaction_rate = 0

        #Store the total reactant density in the reactant_densities dictionary
        reaction_rates[event_index][reaction_name] = reaction_rate
        reaction_rate_constants[event_index][reaction_name] = rate_constant

    for i in h_transfer:
        reaction_rates[event_index][i] = 0
        
    return(reaction_rates)

def choose_reactants(species_df,rxn,rxn_lib): #Randomly chooses an appropriate reactants relevant to the specified reaction
    num_reactants = rxn_lib.kinetic_params[rxn]['reactants']

    if num_reactants==1:
        mask = species_df[rxn+'-Reactant 1']==True
        allowed_reactants = species_df[mask]['Species SMILES'].tolist()
        reactant_propensities = species_df[mask][rxn+'-Total Reaction Sites'].tolist()
        total_propensity = sum(reactant_propensities)
        probabilities = [p/total_propensity for p in reactant_propensities]

        selected_reactant = np.random.choice(allowed_reactants,p=probabilities)
        
        return([selected_reactant])
        
    else:
        mask1 = species_df[rxn+'-Reactant 1']==True
        mask2 = species_df[rxn+'-Reactant 2']==True
        allowed_reactant1 = species_df[mask1]['Species SMILES'].tolist()
        allowed_reactant2 = species_df[mask2]['Species SMILES'].tolist()
        reactant1_propensities = species_df[mask1][rxn+'-Total Reaction Sites'].tolist()
        total_propensity1 = sum(reactant1_propensities)
        probabilities1 = [p/total_propensity1 for p in reactant1_propensities]

        if allowed_reactant1:
            selected_reactant1 = np.random.choice(allowed_reactant1,p=probabilities1)
        else:
            return
        
        for species in allowed_reactant2:#Ensures that reactants that can participate as both reactants but have a population<2 do not get selected twice
            if species==selected_reactant1:
                mask_sp = species_df['Species SMILES']==species
                pop = species_df[mask_sp]['Current Population']
                if pop.iloc[0] >= 2:
                    continue
                else:
                    allowed_reactant2.remove(species)
            else:
                continue
                
        reactant2_propensities = species_df[mask2][rxn+'-Total Reaction Sites'].tolist()
        total_propensity2 = sum(reactant2_propensities)
        probabilities2 = [p/total_propensity2 for p in reactant2_propensities]

        if allowed_reactant2:
            selected_reactant2 = np.random.choice(allowed_reactant2,p=probabilities2)

            return([selected_reactant1, selected_reactant2])
        else:
            return

def choose_reaction(reaction_rates,event_index): #Chooses a random reaction based on current reaction rates and calculates reaction time
    rxn_list = ["rxn{}".format(i) for i in range(1, 64)]
    rates = list(reaction_rates[event_index].values())
    total_rate = sum(rates)
    probabilities = [p/total_rate for p in rates]

    selected_reaction = np.random.choice(rxn_list,p=probabilities)
    reaction_time = -np.log(np.random.rand())/total_rate
    
    return([selected_reaction,reaction_time])

def kmc_sim(temperature): #Main Algorithm
    compute_time_start = time.time()
    
    #Initialise dataframe, dictionaries and variables
    species_df = pd.DataFrame(columns=['Species SMILES','Current Population','Molecular Weight'])
    species_pop = {}
    swept_species = {}
    reaction_rates = {0:{'rxn1':0}}
    events={0:{'reaction':'','reactants':[],'products':[],'tau':0, 'time elapsed':0, 'conversion':0, 'average chain length':0}}
    timing = {'Reaction Rate Calc':[],'Choose Reaction':[],'Log Events':[],'Reaction':[],'Update Population':[],'Update DF':[]}
    rxn_lib = ReactionLib()
    current_time = 0
    event_num = 0
     
    #Initialise species present in system at t=0
    c_chain = 'C'*INIT_CARBON_CHAIN_LENGTH
    initial_species = {c_chain:INIT_FEEDSTOCK_SIZE,'[H]':1,}
    swept_species,species_df,species_pop = initialize_species(swept_species,species_pop,species_df,initial_species)

    #Define parameters for end point (either complete conversion or average product distribution can be used as endpoints)
    feedstock_pop = species_df.loc[species_df['Species SMILES']==c_chain].iloc[0]['Current Population']
    len_df = pd.DataFrame()
    len_df['Chain Length'] = species_df['Chain Length']
    pop_list=[]
    pop_list = [i[-1] for i in swept_species.values()]
    len_df['Population'] = pop_list
    weighted_len = (len_df['Chain Length']*len_df['Population']).sum()
    total_pop = len_df.loc[(len_df['Population']>0)&(len_df['Chain Length']>0)]['Population'].sum()
    average_len = weighted_len/total_pop

    #Name main folder for output
    temp_folder = f'Temp_{temperature-273}C'
    os.makedirs(RUN_FOLDER, exist_ok=True)
    
    #Start loop
    while average_len>CHAIN_LENGTH_ENDPOINT or feedstock_pop>0:#Endpoint condition
        
        start_time = time.time()
        reaction_rates = calc_reaction_rate(species_df,reaction_rates,rxn_lib,temperature,event_num)
        reaction_rate_time = time.time()
        chosen_reaction,tau=choose_reaction(reaction_rates,event_num)
        chosen_reactants = choose_reactants(species_df,chosen_reaction,rxn_lib)
        while not chosen_reactants:
            chosen_reaction,tau = choose_reaction(reaction_rates,event_num)
            chosen_reactants = choose_reactants(species_df,chosen_reaction,rxn_lib)
        choose_reaction_time = time.time()


        events[event_num+1]={}
        events[event_num+1]['reaction']=chosen_reaction
        events[event_num+1]['reactants']=chosen_reactants
        events[event_num+1]['tau']=tau
        current_time += tau
        events[event_num+1]['time elapsed']=current_time
        log_events_time = time.time()

        products = rxn_lib.reactions[chosen_reaction](chosen_reactants)

        rxn_time = time.time()
        event_num += 1
        

        events[event_num]['products']=products

        


        swept_species,species_pop=update_population(species_pop,swept_species,events,event_num)

        pop_time = time.time()
        species_time = time.time()
        for product in products:
            species_df = update_species(species_df,species_pop,product)
            species_time = time.time()

        feedstock_pop = species_df.loc[species_df['Species SMILES']==c_chain].iloc[0]['Current Population']
        len_df = pd.DataFrame()
        len_df['Chain Length'] = species_df['Chain Length']
        pop_list=[]
        pop_list = [i[-1] for i in swept_species.values()]
        len_df['Population'] = pop_list
        weighted_len = (len_df['Chain Length']*len_df['Population']).sum()
        total_pop = len_df.loc[(len_df['Population']>0)&(len_df['Chain Length']>0)]['Population'].sum()
        average_len = weighted_len/total_pop
        
        timing['Reaction Rate Calc'] += [reaction_rate_time-start_time]
        timing['Choose Reaction'] += [choose_reaction_time-reaction_rate_time]
        timing['Log Events'] += [log_events_time-choose_reaction_time]
        timing['Reaction'] += [rxn_time-log_events_time]
        timing['Update Population'] += [pop_time-rxn_time]
        timing['Update DF'] += [species_time-pop_time]

        if event_num%SWEEP_FREQ == 0:
            species_df,species_pop = n2_sweep(species_df,species_pop)

        if event_num%WRITE_FREQ == 0:
            iteration_folder = os.path.join(RUN_FOLDER, temp_folder, f'Iteration_{event_num}')
            os.makedirs(iteration_folder, exist_ok=True)

            # Write DataFrames and dictionaries to output files in the iteration folder
            species_df.to_csv(os.path.join(iteration_folder, 'species list.csv'), index=False)
            with open(os.path.join(iteration_folder, 'population history.json'), 'w') as json_file:
                json.dump(species_pop, json_file)
            with open(os.path.join(iteration_folder, 'reaction rates.json'), 'w') as json_file:
                json.dump(reaction_rates, json_file)
            with open(os.path.join(iteration_folder, 'events.json'), 'w') as json_file:
                json.dump(events, json_file)
            with open(os.path.join(iteration_folder, 'timing.json'), 'w') as json_file:
                json.dump(timing, json_file)
            with open(os.path.join(iteration_folder, 'swept species.json'), 'w') as json_file:
                json.dump(swept_species, json_file)
                
        current_compute_time = time.time()-compute_time_start
    
        if event_num%PRINT_FREQ==0 or event_num==100:
            print(f'For Temperature = {temperature} - Event num: {event_num}, Reaction time elapsed: {current_time}, Compute time elapsed: {current_compute_time}, Avg chain length: {average_len}, Feedstock conversion: {(INIT_FEEDSTOCK_SIZE-feedstock_pop)/INIT_FEEDSTOCK_SIZE}')

    iteration_folder = os.path.join(RUN_FOLDER, temp_folder, f'final_iteration_{event_num}')
    os.makedirs(iteration_folder, exist_ok=True)

    #Write DataFrames and dictionaries to output files at endpoint
    species_df.to_csv(os.path.join(iteration_folder, 'species list.csv'), index=False)
    with open(os.path.join(iteration_folder, 'population history.json'), 'w') as json_file:
        json.dump(species_pop, json_file)
    with open(os.path.join(iteration_folder, 'reaction rates.json'), 'w') as json_file:
        json.dump(reaction_rates, json_file)
    with open(os.path.join(iteration_folder, 'events.json'), 'w') as json_file:
        json.dump(events, json_file)
    with open(os.path.join(iteration_folder, 'timing.json'), 'w') as json_file:
        json.dump(timing, json_file)
    with open(os.path.join(iteration_folder, 'swept species.json'), 'w') as json_file:
        json.dump(swept_species, json_file)