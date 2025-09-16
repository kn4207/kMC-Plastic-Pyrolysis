from rdkit import Chem
from collections import Counter
from config import MW_SWEEP_THRESHOLD

def atom_metadata(species): #Used to obtain atomic metadata of a molecule such as hybridization
    mol = Chem.MolFromSmiles(species)
    metadata = {
        'atoms': [atom.GetSymbol() for atom in mol.GetAtoms()],
        'hybridization': [str(atom.GetHybridization()) for atom in mol.GetAtoms()],
        'hydrogens': [atom.GetTotalNumHs() for atom in mol.GetAtoms()],
        'degree': [atom.GetTotalDegree() for atom in mol.GetAtoms()],
        'radicals': [atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()],
    }

    sp3 = {i for i, (hyb, hyd, deg, rads) in enumerate(zip(metadata['hybridization'], metadata['hydrogens'], metadata['degree'], metadata['radicals'])) 
           if hyb == 'SP3' and hyd >= 0 and deg == 4 and rads == 0}
    sp2 = {i for i, (hyb, hyd, deg, rads) in enumerate(zip(metadata['hybridization'], metadata['hydrogens'], metadata['degree'], metadata['radicals'])) 
           if hyb == 'SP2' and hyd >= 0 and deg == 3 and rads == 0}
    rad = {i for i, (hyb, hyd, deg, rads) in enumerate(zip(metadata['hybridization'], metadata['hydrogens'], metadata['degree'], metadata['radicals'])) 
           if hyb == 'SP3' and hyd >= 0 and deg == 3 and rads == 1}

    metadata.update({
        'Sp3': sp3,
        'Sp2': sp2,
        'Rad': rad,
    })

    metadata.update({
        'Psp3': sp3 & {i for i, hyd in enumerate(metadata['hydrogens']) if hyd == 3},
        'Ssp3': sp3 & {i for i, hyd in enumerate(metadata['hydrogens']) if hyd == 2},
        'Tsp3': sp3 & {i for i, hyd in enumerate(metadata['hydrogens']) if hyd == 1},
        'Qsp3': sp3 & {i for i, hyd in enumerate(metadata['hydrogens']) if hyd == 0},
        'Usp2': sp2 & {i for i, hyd in enumerate(metadata['hydrogens']) if hyd == 2},
        'Msp2': sp2 & {i for i, hyd in enumerate(metadata['hydrogens']) if hyd == 1},
        'Dsp2': sp2 & {i for i, hyd in enumerate(metadata['hydrogens']) if hyd == 0},
    })

    metadata.update({
        'Prad': rad & {i for i, hyd in enumerate(metadata['hydrogens']) if hyd == 2},
        'Srad': rad & {i for i, hyd in enumerate(metadata['hydrogens']) if hyd == 1},
    })

    return metadata

def update_population(reactive_species_pop,all_species_pop,events,event_index): #Checks if species is new and adds it to species_pop and enumerates population history as 0 until the event it is added
    reactant_counts = Counter(events[event_index]['reactants'])
    product_counts = Counter(events[event_index]['products'])

    for species in product_counts: #Add population of newly produced species
        if species not in all_species_pop:
            all_species_pop[species] = [0]*(event_index)
            reactive_species_pop[species] = [0]*(event_index)

    for species in all_species_pop.keys():
        prev_all = all_species_pop.get(species,[0])[-1]  #Initialize with previous population
        prev_reactive = reactive_species_pop.get(species,[0])[-1]

        new_all = prev_all - reactant_counts.get(species, 0) + product_counts.get(species, 0)
        new_reactive = prev_reactive - reactant_counts.get(species, 0) + product_counts.get(species, 0)

        all_species_pop[species].append(new_all)
        reactive_species_pop[species].append(new_reactive)

    return all_species_pop, reactive_species_pop

def n2_sweep(species,reactive_pop): #Remove species below a certain threshold for mol. wt.

    mask = (species['Current Population'] > 0) & (species['Molecular Weight'] < MW_SWEEP_THRESHOLD) & (species['Radical?'] == False)

    #Filter the DataFrame using the mask
    filtered_species = species[mask]
    
    #Update the separate dictionary with the filtered species
    for _, row in filtered_species.iterrows():
        species_smiles = row['Species SMILES']
        reactive_pop[species_smiles][-1] = 0
    
    #Set the population to zero in the DataFrame using the mask
    species.loc[mask,'Current Population'] = 0

    # print(f"Swept {len(filtered_species)} species below {MW_SWEEP_THRESHOLD} g/mol")

    
    return(species,reactive_pop)