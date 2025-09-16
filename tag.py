from utils import atom_metadata
import re

def tag_polymer(pspecies): #Used to assign reaction tags to polymer species (without radicals present) to specify which reactions it can participate in 
    tags = {}
    tags_df = {}
    rxns = ["rxn{}".format(i) for i in range(1, 64)]
    for i in rxns:
        tags.update({i:{'Reactant 1':False, 'Reactant 2':False, 'Reaction Site Density':0, 'Total Reaction Sites':0 }})
    adata = atom_metadata(pspecies)
    plen = len(adata['atoms'])

    #rxn1
    if pspecies=='CC':
        tags['rxn1']['Reactant 1']=True
        tags['rxn1']['Reaction Site Density']+=1   
    #rxn2
    if not adata['Rad'] and plen>2 and adata['Psp3']: #If no radicals are present, polymer length is >2 and there are primary sp3 carbons present
        for i in adata['Psp3']:
            if i+1 in adata['Ssp3'] or i-1 in adata['Ssp3']: #If a Primary sp3 carbon is connected to a Secondary sp3 carbon
                tags['rxn2']['Reactant 1']=True
                tags['rxn2']['Reaction Site Density']+=1  #Counts number of possible reaction sites
    #rxn3
    if not adata['Rad'] and adata['Ssp3']:
        for i in adata['Ssp3']:
            tags['rxn7']['Reactant 1']=True
            tags['rxn7']['Reaction Site Density']+=1
            tags['rxn16']['Reactant 1']=True
            tags['rxn16']['Reaction Site Density']+=1 
            tags['rxn18']['Reactant 1']=True
            tags['rxn18']['Reaction Site Density']+=1
            tags['rxn20']['Reactant 1']=True
            tags['rxn20']['Reaction Site Density']+=1
            tags['rxn26']['Reactant 1']=True
            tags['rxn26']['Reaction Site Density']+=1
            if i+1 in adata['Ssp3']:
                tags['rxn3']['Reactant 1']=True
                tags['rxn3']['Reaction Site Density']+=1 
    #rxn5
    if pspecies=='C':
        tags['rxn5']['Reactant 1']=True
        tags['rxn5']['Reaction Site Density']+=4  #4 C-H bonds can be broken
        tags['rxn21']['Reactant 2']=True
        tags['rxn21']['Reaction Site Density']+=4 
        tags['rxn22']['Reactant 1']=True
        tags['rxn22']['Reaction Site Density']+=4
        tags['rxn27']['Reactant 2']=True
        tags['rxn27']['Reaction Site Density']+=4
    #rxn6
    if not adata['Rad'] and adata['Psp3']:
        for i in adata['Psp3']:
            tags['rxn6']['Reactant 1']=True
            tags['rxn6']['Reaction Site Density']+=1
            tags['rxn15']['Reactant 1']=True
            tags['rxn15']['Reaction Site Density']+=1 
            tags['rxn17']['Reactant 1']=True
            tags['rxn17']['Reaction Site Density']+=1 
            tags['rxn19']['Reactant 1']=True
            tags['rxn19']['Reaction Site Density']+=1
            tags['rxn25']['Reactant 1']=True
            tags['rxn25']['Reaction Site Density']+=1
    #rxn56
    if pspecies=='C=C':
        tags['rxn56']['Reactant 1']=True
        tags['rxn56']['Reaction Site Density']+=2
        tags['rxn60']['Reactant 1']=True
        tags['rxn60']['Reaction Site Density']+=2
        tags['rxn62']['Reactant 1']=True
        tags['rxn62']['Reaction Site Density']+=2
    #rxn57
    if not adata['Rad'] and plen>2 and adata['Usp2']: #If no radicals are present, polymer length is >2 and there are primary sp2 carbons present
        for i in adata['Usp2']:
            if i+1 in adata['Msp2'] or i-1 in adata['Msp2']: #If a Primary sp2 carbon is connected to a Secondary sp2 carbon
                tags['rxn57']['Reactant 1']=True
                tags['rxn57']['Reaction Site Density']+=1 
                tags['rxn58']['Reactant 1']=True
                tags['rxn58']['Reaction Site Density']+=1
                tags['rxn61']['Reactant 1']=True
                tags['rxn61']['Reaction Site Density']+=1
                tags['rxn63']['Reactant 1']=True
                tags['rxn63']['Reaction Site Density']+=1
    #rxn59
    if not adata['Rad'] and adata['Msp2']:
        for i in adata['Msp2']:
            if i+1 in adata['Msp2']:
                tags['rxn59']['Reactant 1']=True
                tags['rxn59']['Reaction Site Density']+=1       

    for rxn in tags.keys():
        for tag,value in tags[rxn].items():
            n = rxn+'-'+tag
            tags_df[n]=value

    return(tags_df)

def tag_polymer_radical(prspecies): #Used to assign reaction tags to polymer radical species
    tags = {}
    tags_df = {}
    rxns = ["rxn{}".format(i) for i in range(1, 64)]
    for i in rxns:
        tags.update({i:{'Reactant 1':False, 'Reactant 2':False, 'Reaction Site Density':0, 'Total Reaction Sites':0}})
    adata = atom_metadata(prspecies)
    plen = len(adata['atoms'])
    c_cchains=re.split('\[CH\]|C=C',prspecies) #type:ignore
    chainlength = len(max(c_cchains, key=len))

    #rxn10
    if adata['Prad']:
        for i in adata['Prad']:
            tags['rxn10']['Reactant 1']=True
            tags['rxn10']['Reaction Site Density']+=1
            tags['rxn13']['Reactant 1']=True
            tags['rxn13']['Reaction Site Density']+=1
            tags['rxn14']['Reactant 1']=True
            tags['rxn14']['Reactant 2']=2
            tags['rxn14']['Reaction Site Density']+=1
            tags['rxn15']['Reactant 2']=True
            tags['rxn15']['Reaction Site Density']+=1
            tags['rxn16']['Reactant 2']=True
            tags['rxn16']['Reaction Site Density']+=1
            tags['rxn21']['Reactant 1']=True
            tags['rxn21']['Reaction Site Density']+=1
            tags['rxn23']['Reactant 1']=True
            tags['rxn23']['Reaction Site Density']+=1
            tags['rxn47']['Reactant 2']=True
            tags['rxn47']['Reaction Site Density']+=1
            tags['rxn62']['Reactant 2']=True
            tags['rxn62']['Reaction Site Density']+=1
            tags['rxn63']['Reactant 2']=True
            tags['rxn63']['Reaction Site Density']+=1
            if i+1 in adata['Ssp3'] or i-1 in adata['Ssp3']:
                tags['rxn29']['Reactant 1']=True
                tags['rxn29']['Reaction Site Density']+=1
                tags['rxn42']['Reactant 1']=True
                tags['rxn42']['Reaction Site Density']+=1
                tags['rxn44']['Reactant 1']=True
                tags['rxn44']['Reaction Site Density']+=1
                tags['rxn48']['Reactant 1']=True
                tags['rxn48']['Reaction Site Density']+=1
                tags['rxn50']['Reactant 1']=True
                tags['rxn50']['Reaction Site Density']+=1
            if i+2 in adata['Ssp3'] or i-2 in adata['Ssp3']:
                tags['rxn30']['Reactant 1']=True
                tags['rxn30']['Reaction Site Density']+=1
            if (i+1 in adata['Ssp3'] and i+2 in adata['Ssp3']) or (i-1 in adata['Ssp3'] and i-2 in adata['Ssp3']):
                tags['rxn53']['Reactant 1']=True
                tags['rxn53']['Reaction Site Density']+=1
    #rxn11
    if adata['Srad']:
        for i in adata['Srad']:
            tags['rxn11']['Reactant 1']=True
            tags['rxn11']['Reaction Site Density']+=1
            tags['rxn25']['Reactant 2']=True
            tags['rxn25']['Reaction Site Density']+=1
            tags['rxn26']['Reactant 2']=True
            tags['rxn26']['Reaction Site Density']+=1
            tags['rxn27']['Reactant 1']=True
            tags['rxn27']['Reaction Site Density']+=1
            tags['rxn28']['Reactant 1']=True
            tags['rxn28']['Reaction Site Density']+=1
            tags['rxn48']['Reactant 2']=True
            tags['rxn48']['Reaction Site Density']+=1
            if i+1 in adata['Psp3'] or i-1 in adata['Psp3']:
                tags['rxn35']['Reactant 1']=True
                tags['rxn35']['Reaction Site Density']+=1  
            if i+1 in adata['Ssp3'] or i-1 in adata['Ssp3']:
                tags['rxn36']['Reactant 1']=True
                tags['rxn36']['Reaction Site Density']+=1
            if i+2 in adata['Ssp3'] or i-2 in adata['Ssp3']:
                tags['rxn37']['Reactant 1']=True
                tags['rxn37']['Reaction Site Density']+=1
            if i+1 in adata['Ssp3']:
                tags['rxn43']['Reactant 1']=True
                tags['rxn43']['Reaction Site Density']+=1
                tags['rxn45']['Reactant 1']=True
                tags['rxn45']['Reaction Site Density']+=1
                tags['rxn47']['Reactant 1']=True
                tags['rxn47']['Reaction Site Density']+=1
                tags['rxn51']['Reactant 1']=True
                tags['rxn51']['Reaction Site Density']+=1
            if i-1 in adata['Ssp3']:
                tags['rxn43']['Reactant 1']=True
                tags['rxn43']['Reaction Site Density']+=1
                tags['rxn45']['Reactant 1']=True
                tags['rxn45']['Reaction Site Density']+=1
                tags['rxn47']['Reactant 1']=True
                tags['rxn47']['Reaction Site Density']+=1
                tags['rxn51']['Reactant 1']=True
                tags['rxn51']['Reaction Site Density']+=1
            if (i+1 in adata['Ssp3'] and i+2 in adata['Ssp3']) or (i-1 in adata['Ssp3'] and i-2 in adata['Ssp3']):
                tags['rxn55']['Reactant 1']=True
                tags['rxn55']['Reaction Site Density']+=1
    #rxn31
    if adata['Prad']and not adata['Msp2'] and not adata['Dsp2']:
        for i in adata['Prad']:
            if i+3 in adata['Ssp3'] or i-3 in adata['Ssp3']:
                tags['rxn31']['Reactant 1']=True
                tags['rxn31']['Reaction Site Density']+=1
            if i+4 in adata['Ssp3'] or i-4 in adata['Ssp3']:
                tags['rxn32']['Reactant 1']=True
                tags['rxn32']['Reaction Site Density']+=1
            if i+5 in adata['Ssp3'] or i-5 in adata['Ssp3']:
                tags['rxn33']['Reactant 1']=True
                tags['rxn33']['Reaction Site Density']+=1
            if i+6 in adata['Ssp3'] or i-6 in adata['Ssp3']:
                tags['rxn34']['Reactant 1']=True
                tags['rxn34']['Reaction Site Density']+=1                        
    #rxn38
    if adata['Srad']and not adata['Msp2'] and not adata['Dsp2']:
        for i in adata['Srad']:
            if i+3 in adata['Ssp3'] or i-3 in adata['Ssp3']:
                tags['rxn38']['Reactant 1']=True
                tags['rxn38']['Reaction Site Density']+=1
            if i+4 in adata['Ssp3'] or i-4 in adata['Ssp3']:
                tags['rxn39']['Reactant 1']=True
                tags['rxn39']['Reaction Site Density']+=1
            if i+5 in adata['Ssp3'] or i-5 in adata['Ssp3']:
                tags['rxn40']['Reactant 1']=True
                tags['rxn40']['Reaction Site Density']+=1
            if i+6 in adata['Ssp3'] or i-6 in adata['Ssp3']:
                tags['rxn41']['Reactant 1']=True
                tags['rxn41']['Reaction Site Density']+=1
    #rxn46
    if adata['Prad']:
        for i in adata['Prad']:
            if i+1 in adata['Ssp3'] or i-1 in adata['Ssp3']:
                tags['rxn46']['Reactant 1']=True
                tags['rxn46']['Reactant 2']=True
                tags['rxn46']['Reaction Site Density']+=1
            else:
                tags['rxn46']['Reactant 2']=True
                tags['rxn46']['Reaction Site Density']+=1           
    #rxn49
    if adata['Srad']:
        for i in adata['Srad']:
            if i+1 in adata['Ssp3']:
                tags['rxn49']['Reactant 1']=True
                tags['rxn49']['Reaction Site Density']+=1
            if i-1 in adata['Ssp3']:
                tags['rxn49']['Reactant 1']=True
                tags['rxn49']['Reaction Site Density']+=1
            else:
                tags['rxn49']['Reactant 2']=True
                tags['rxn49']['Reaction Site Density']+=1                
    #rxn52
    if prspecies=='CC[CH2]' or prspecies=='[CH2]CC':
        tags['rxn52']['Reactant 1']=True
        tags['rxn52']['Reaction Site Density']+=1            
    #rxn54
    if prspecies=='CC[CH]C' or prspecies=='C[CH]CC':
        tags['rxn54']['Reactant 1']=True
        tags['rxn54']['Reaction Site Density']+=1

    for rxn in tags.keys():
        for tag,value in tags[rxn].items():
            n = rxn+'-'+tag
            tags_df[n]=value

    return(tags_df)

def tag_radical(rspecies): #Used to assign reaction tags to radical species
    tags = {}
    tags_df = {}
    rxns = ["rxn{}".format(i) for i in range(1, 64)]
    for i in rxns:
        tags.update({i:{'Reactant 1':False, 'Reactant 2':False, 'Reaction Site Density':0, 'Total Reaction Sites':0}})
    hh = ['rxn23','rxn24','rxn28'] #reactions that h2 particpates in
    h = ['rxn9','rxn10','rxn11','rxn19','rxn20','rxn22','rxn42','rxn43','rxn56','rxn57','rxn58','rxn59'] #Reaction that h' participates in
    ch3 = ['rxn13','rxn17','rxn18','rxn44','rxn45','rxn60','rxn61']
    if rspecies=='[H][H]':
        tags['rxn4']['Reactant 1']=True
        tags['rxn4']['Reaction Site Density']+=1
        for i in hh:
            tags[i]['Reactant 2']=True
            tags[i]['Reaction Site Density']+=1
    elif rspecies=='[H]':
        tags['rxn8']['Reactant 1']=True
        tags['rxn8']['Reactant 2']=True
        tags['rxn8']['Reaction Site Density']+=1
        for i in h:
            tags[i]['Reactant 2']=True
            tags[i]['Reaction Site Density']+=1
    elif rspecies=='[CH3]':
        tags['rxn9']['Reactant 1']=True
        tags['rxn9']['Reaction Site Density']+=1
        tags['rxn12']['Reactant 1']=True
        tags['rxn12']['Reactant 2']=True
        tags['rxn12']['Reaction Site Density']+=1
        tags['rxn24']['Reactant 1']=True
        tags['rxn24']['Reaction Site Density']+=1
        for i in ch3:
            tags[i]['Reactant 2']=True
            tags[i]['Reaction Site Density']+=1

    

    for rxn in tags.keys():
        for tag,value in tags[rxn].items():
            n = rxn+'-'+tag
            tags_df[n]=value

    return(tags_df)