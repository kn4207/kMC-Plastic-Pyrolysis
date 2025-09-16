import random
import json
from rdkit import Chem
from rdkit.Chem.rdChemReactions import ReactionFromSmarts

class ReactionLib:
    def __init__(self):
        with open('kinetic_params.json') as f:
            self.kinetic_params = json.load(f)
        
        self.reactions = {
            name:getattr(self,name) for name in dir(self)
            if callable(getattr(self,name)) and name.startswith('rxn')
        }
        
    #Initiation Reactions

    #Reaction 1 - C-C bond activation (Primary and Primary)
    def rxn1(self, reactant): #Pass a list containing SMILES of reactants
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn =  ReactionFromSmarts('[C;H3:1][CH3:2]>>[CH3:1].[CH3:2]')
        rad1 = rxn.RunReactants((rlist[0], ))[0][0]
        rad2 = rxn.RunReactants((rlist[0], ))[0][1]
        Chem.SanitizeMol(rad1) #Necessary to represent as radicals
        Chem.SanitizeMol(rad2)
        products = [Chem.MolToSmiles(rad1),Chem.MolToSmiles(rad2)]
        return (products) #Returns a list containing SMILES of products [CH3',CH3']
    
    #Reaction 2 - C-C bond activation (Primary and Secondary)
    def rxn2(self, reactant): #[RCH3]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn =  ReactionFromSmarts('[X4;H2:1][X4;H3:2]>>[CH2:1].[CH3:2]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        rad1 = prod[idx][0]
        rad2 = prod[idx][1]
        Chem.SanitizeMol(rad1)
        Chem.SanitizeMol(rad2)
        products = [Chem.MolToSmiles(rad1),Chem.MolToSmiles(rad2)]
        return (products) #[RCH2',CH3']
    
    #Reaction 3 - C-C bond activation (Secondary and Secondary)
    def rxn3(self, reactant):
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn =  ReactionFromSmarts('[X4;H2:1][X4;H2:2]>>[CH2:1].[CH2:2]')
        prod = rxn.RunReactants((rlist[0], )) #Generates all possible combinations of secondary C-C bond activation
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1) #Chooses a random outcome of secondary C-C bond activation
        rad1 = prod[idx][0]
        rad2 = prod[idx][1]
        Chem.SanitizeMol(rad1)
        Chem.SanitizeMol(rad2)
        products = [Chem.MolToSmiles(rad1),Chem.MolToSmiles(rad2)]
        return (products) #[RCH2',RCH2'] where R can be of varying lengths
    
    #Reaction 4 - H-H bond activation
    def rxn4(self, reactant):
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn =  ReactionFromSmarts('[H:1][H:2]>>[H:1].[H:2]')
        rad1 = rxn.RunReactants((rlist[0], ))[0][0]
        rad2 = rxn.RunReactants((rlist[0], ))[0][1]
        rad1.GetAtomWithIdx(0).SetNoImplicit(True) #Allows for H to be represented as single atom instead of H2
        rad2.GetAtomWithIdx(0).SetNoImplicit(True)
        Chem.SanitizeMol(rad1)
        Chem.SanitizeMol(rad2)
        products = [Chem.MolToSmiles(rad1),Chem.MolToSmiles(rad2)]
        return (products) #[H,H]
    
    #Reaction 5 - C-H bond activation (Methyl Carbon)
    def rxn5(self, reactant):
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn =  ReactionFromSmarts('[C;H4:1]>>[CH3:1].[H:2]')
        rad1 = rxn.RunReactants((rlist[0], ))[0][0]
        rad2 = rxn.RunReactants((rlist[0], ))[0][1]
        rad2.GetAtomWithIdx(0).SetNoImplicit(True)
        Chem.SanitizeMol(rad1)
        Chem.SanitizeMol(rad2)
        products = [Chem.MolToSmiles(rad1),Chem.MolToSmiles(rad2)]
        return (products) #[CH3',H]
    
    #Reaction 6 - C-H bond activation (Primary)
    def rxn6(self, reactant):
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn =  ReactionFromSmarts('[C;H3:1][CH2:2]>>[CH2:2].[H:1]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        rad1 = prod[idx][0]
        rad2 = prod[idx][1]
        rad2.GetAtomWithIdx(0).SetNoImplicit(True)
        Chem.SanitizeMol(rad1)
        Chem.SanitizeMol(rad2)
        products = [Chem.MolToSmiles(rad1),Chem.MolToSmiles(rad2)]
        return (products) #[RCH2',H]
        
    #Reaction 7 - C-H bond activation (Secondary)
    def rxn7(self, reactant):
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn =  ReactionFromSmarts('[X4;H2:1]>>[CH:1].[H:2]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        rad1 = prod[idx][0]
        rad2 = prod[idx][1]
        rad2.GetAtomWithIdx(0).SetNoImplicit(True)
        Chem.SanitizeMol(rad1)
        Chem.SanitizeMol(rad2)
        products = [Chem.MolToSmiles(rad1),Chem.MolToSmiles(rad2)]
        return (products) #[RCH'R,H]
    
    #Radical Recombination Reactions
    
    #Reaction 8 - H-H bond fusion
    def rxn8(self, reactants): #List should contain [H,H]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[H:1].[H:2]>>[H:1][H:2]')
        prod = rxn.RunReactants((rlist[0],rlist[1]))[0][0]
        Chem.SanitizeMol(prod)
        product = [Chem.MolToSmiles(prod)]
        return (product) #[H2]
    
    #Reaction 9 - C-H bond fusion (Methyl Carbon)
    def rxn9(self, reactants): #List should contain [CH3',H]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[CH3:1]>>[CH4:1]')
        prod = rxn.RunReactants((rlist[0], ))[0][0]
        Chem.SanitizeMol(prod)
        product = [Chem.MolToSmiles(prod)]
        return (product) #[CH4]
    
    
    #Reaction 10 - C-H bond fusion (Primary Carbon)
    def rxn10(self, reactants): #[RCH2',H] - radical must be primary carbon
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[v3;H2:1]>>[CH3:1]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return (product) #[RCH3]
    
    #Reaction 11 - C-H bond fusion (Secondary Carbon)
    def rxn11(self, reactants): #[RCH'R, H] - radical must be secondary carbon
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[v3;H1:1]>>[CH2:1]') #identifies carbon atoms with a bond order of 3 to distinguish between R[CH]R and RCH=CH2
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return (product) #[RCH2R]
    
    #Reaction 12 - C-C bond fusion (Primary and Primary)
    def rxn12(self, reactants): #[CH3',CH3']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[CH3:1].[CH3:2]>>[CH3:1][CH3:2]')
        prod = rxn.RunReactants((rlist[0],rlist[1] ))[0][0]
        Chem.SanitizeMol(prod)
        product = [Chem.MolToSmiles(prod)]
        return (product) #[CH3-CH3]
    
    #Reaction 13 - C-C bond fusion (Primary and Secondary)
    def rxn13(self, reactants): #[RCH2',CH3']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[v3;H2:1].[CH3:2]>>[CH2:1][CH3:2]')
        prod = rxn.RunReactants((rlist[0], rlist[1]))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return (product) #[RCH2CH3]
    
    #Reaction 14 - C-C bond fusion (Secondary and Secondary)
    def rxn14(self, reactants): #[RCH2',RCH2']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[v3;H2:1].[v3;H2:2]>>[CH2:1][CH2:2]')
        prod = rxn.RunReactants((rlist[0], rlist[1]))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return (product) #[RCH2CH2R]
    
    #Hydrogen Abstraction Reactions
    
    #Reaction 15 - H abstraction (By Primary Radical from Primary Carbon)
    def rxn15(self, reactants): #[RCH3,RCH2']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[X4;H3:1].[v3;H2:2]>>[CH2:1].[CH3:2]')
        prod = rxn.RunReactants((rlist[0],rlist[1] ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1),Chem.MolToSmiles(prod2)]
        return (product) #[RCH2',RCH3]
    
    #Reaction 16 - H abstraction (By Primary Radical from Secondary Carbon)
    def rxn16(self, reactants): #[RCH2R,RCH2']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[X4;H2:1].[v3;H2:2]>>[CH:1].[CH3:2]')
        prod = rxn.RunReactants((rlist[0],rlist[1] ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1),Chem.MolToSmiles(prod2)]
        return (product) #[RCH'R,RCH3]
    
    #Reaction 17 - H abstraction (Methyl radical from Primary Carbon)
    def rxn17(self, reactants): #[RCH3,CH3']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[X4;H3:1].[X3;H3:2]>>[CH2:1].[CH4:2]')
        prod = rxn.RunReactants((rlist[0],rlist[1] ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1),Chem.MolToSmiles(prod2)]
        return (product) #[RCH2',CH4]
    
    #Reaction 18 - H abstraction (Methyl radical from Secondary Carbon)
    def rxn18(self, reactants): #[RCH2R,CH3']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[X4;H2:1].[X3;H3:2]>>[CH:1].[CH4:2]')
        prod = rxn.RunReactants((rlist[0],rlist[1] ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1),Chem.MolToSmiles(prod2)]
        return(product)#[RCH'R,CH4]
    
    #Reaction 19 - H abstraction (H from Primary Carbon)
    def rxn19(self, reactants): #[RCH3, H]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn  = ReactionFromSmarts('[X4;H3:1]>>[CH2:1]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1), '[H][H]']
        return(product) #[RCH2',H2]
    
    #Reaction 20 - H abstraction (H from Secondary Carbon)
    def rxn20(self, reactants): #[RCH2R, H]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn  = ReactionFromSmarts('[X4;H2:1]>>[CH:1]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1), '[H][H]']
        return(product) #[RCH'R,H2]
    
    #Reaction 21 - H abstraction (Primary radical from Methane)
    def rxn21(self, reactants): #[RCH2', CH4]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn  = ReactionFromSmarts('[v3;H2:1].[CH4:2]>>[CH3:1].[CH3:2]')
        prod = rxn.RunReactants((rlist[0], rlist[1]))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1), Chem.MolToSmiles(prod2)]
        return(product) #[RCH3,CH3']
    
    #Reaction 22 - H abstraction (H radical from methane)
    def rxn22(self, reactants): #[CH4, H]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn  = ReactionFromSmarts('[CH4:1]>>[CH3:1]')
        prod = rxn.RunReactants((rlist[0], ))[0][0]
        Chem.SanitizeMol(prod)
        product = [Chem.MolToSmiles(prod), '[H][H]']
        return(product) #[CH3',H2]
    
    #Reaction 23 - H abstraction (Primary radical from H2)
    def rxn23(self, reactants): #[RCH2',H2]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn  = ReactionFromSmarts('[v3;H2:1]>>[CH3:1]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1), '[H]']
        return(product) #[RCH3,H]
    
    #Reaction 24 - H abstraction (Methyl radical from H2)
    def rxn24(self, reactants): #[CH3',H2]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn  = ReactionFromSmarts('[X3;H3:1]>>[CH4:1]')
        prod = rxn.RunReactants((rlist[0], ))[0][0]
        Chem.SanitizeMol(prod)
        product = [Chem.MolToSmiles(prod), '[H]']
        return(product) #[CH4,H]
    
    #Reaction 25 - H abstraction (Secondary radical from Primary Carbon)
    def rxn25(self, reactants): #[RCH3, RCH'R]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[X4;H3:1].[v3;H1:2]>>[CH2:1].[CH2:2]')
        prod = rxn.RunReactants((rlist[0],rlist[1]))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1), Chem.MolToSmiles(prod2)]
        return(product) #[RCH2',RCH2R]
    
    #Reaction 26 - H abstraction (Secondary Radical from Secondary Carbon)
    def rxn26(self, reactants): #[RCH2R, RCH'R]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[X4;H2:1].[v3;H1:2]>>[CH:1].[CH2:2]')
        prod = rxn.RunReactants((rlist[0],rlist[1]))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1), Chem.MolToSmiles(prod2)]
        return(product) #[RCH'R,RCH2R]
    
    #Reaction 27 - H abstraction (Secondary Radical from Methane)
    def rxn27(self, reactants): #[RCH'R,CH4]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[v3;H1:1].[X4;H4:2]>>[CH2:1].[CH3:2]')
        prod = rxn.RunReactants((rlist[0],rlist[1]))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1), Chem.MolToSmiles(prod2)]
        return(product) #[RCH2R,CH3]
    
    #Reaction 28 - H abstraction (Secondary Radical from H2)
    def rxn28(self, reactants): #[RCH'R,H2]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[v3;H1:1]>>[CH2:1]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1), '[H]']
        return(product) #[RCH2R,H]
    
    #Backbiting Reactions
    
    #Reaction 29 - Intramolecular H transfer (Primary Carbon C1 from Secondary Carbon C2)
    def rxn29(self, reactant): #[RCH2']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[X4;H2:1][v3;H2:2]>>[CH:1][CH3:2]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH'CH3]
    
    #Reaction 30 - Intramolecular H transfer (Primary Carbon C1 from Secondary Carbon C3)
    def rxn30(self, reactant): #[RCH2']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[X4;H2:1][X4;H2:2][v3;H2:3]>>[CH:1][CH2:2][CH3:3]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH'CH2CH3]
    
    #Reaction 31 - Intramolecular H transfer (Primary Carbon C1 from Secondary Carbon C4)
    def rxn31(self, reactant): #[RCH2']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[X4;H2:1][X4;H2:2][X4;H2:3][v3;H2:4]>>[CH:1][CH2:2][CH2:3][CH3:4]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH'CH2CH2CH3]
    
    #Reaction 32 - Intramolecular H transfer (Primary Carbon C1 from Secondary Carbon C5)
    def rxn32(self, reactant): #[RCH2']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[X4;H2:1][X4;H2:2][X4;H2:3][X4;H2:4][v3;H2:5]>>[CH:1][CH2:2][CH2:3][CH2:4][CH3:5]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH'CH2CH2CH2CH3]
    
    #Reaction 33 - Intramolecular H transfer (Primary Carbon C1 from Secondary Carbon C6)
    def rxn33(self, reactant): #[RCH2']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[X4;H2:1][X4;H2:2][X4;H2:3][X4;H2:4][X4;H2:5][v3;H2:6]>>[CH:1][CH2:2][CH2:3][CH2:4][CH2:5][CH3:6]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH'CH2CH2CH2CH3]
    
    #Reaction 34 - Intramolecular H transfer (Primary Carbon C1 from Secondary Carbon C7)
    def rxn34(self, reactant): #[RCH2']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[X4;H2:1][X4;H2:2][X4;H2:3][X4;H2:4][X4;H2:5][X4;H2:6][v3;H2:7]>>[CH:1][CH2:2][CH2:3][CH2:4][CH2:5][CH2:6][CH3:7]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH'CH2CH2CH2CH2CH3]
    
    #Reaction 35 - Intramolecular H transfer (Secondary Carbon C2 from Primary Carbon C1)
    def rxn35(self, reactant): #[RCH'CH3]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[v3;H1:1][X4;H3:2]>>[CH2:1][CH2:2]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product)
    
    #Reaction 36 - Intramolecular H transfer (Secondary Carbon x from Secondary Carbon x+1)
    def rxn36(self, reactant): #[RCH'R]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[v3;H1:1][X4;H2:2]>>[CH2:1][CH:2]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH'CH2R]
    
    #Reaction 37 - Intramolecular H transfer (Secondary Carbon x from Secondary Carbon x+2)
    def rxn37(self, reactant): #[RCH'R]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[v3;H1:1][X4;H2:2][X4;H2:3]>>[CH2:1][CH2:2][CH:3]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH'CH2CH2R]
    
    #Reaction 38 - Intramolecular H transfer (Secondary Carbon x from Secondary Carbon x+3)
    def rxn38(self, reactant): #[RCH'R]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[v3;H1:1][X4;H2:2][X4;H2:3][X4;H2:4]>>[CH2:1][CH2:2][CH2:3][CH:4]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH'CH2CH2CH2R]
    
    #Reaction 39 - Intramolecular H transfer (Secondary Carbon x from Secondary Carbon x+4)
    def rxn39(self, reactant): #[RCH'R]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[v3;H1:1][X4;H2:2][X4;H2:3][X4;H2:4][X4;H2:5]>>[CH2:1][CH2:2][CH2:3][CH2:4][CH:5]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH'CH2CH2CH2CH2R]
    
    #Reaction 40 - Intramolecular H transfer (Secondary Carbon x from Secondary Carbon x+5)
    def rxn40(self, reactant): #[RCH'R]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[v3;H1:1][X4;H2:2][X4;H2:3][X4;H2:4][X4;H2:5][X4;H2:6]>>[CH2:1][CH2:2][CH2:3][CH2:4][CH2:5][CH:6]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH'CH2CH2CH2CH2CH2R]
    
    #Reaction 41 - Intramolecular H transfer (Secondary Carbon x from Secondary Carbon x+6)
    def rxn41(self, reactant): #[RCH'R]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[v3;H1:1][X4;H2:2][X4;H2:3][X4;H2:4][X4;H2:5][X4;H2:6][X4;H2:7]>>[CH2:1][CH2:2][CH2:3][CH2:4][CH2:5][CH2:6][CH:7]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH'CH2CH2CH2CH2CH2CH2R]
    
    
    #Radical Disproportionation Reactions
    
    #Reaction 42 - Disproportionation (Primary Radical and H)
    def rxn42(self, reactants): #[RCH2',H]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[X4;H2:1][v3;H2:2]>>[CH:1]=[CH2:2]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1), '[H][H]']
        return(product) #[RCH=CH2,H2]
    
    #Reaction 43 - Disproportionation (Secondary Radical and H)
    def rxn43(self, reactants): #[RCH'R,H]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[X4;H2:1][v3;H1:2]>>[CH:1]=[CH:2]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1), '[H][H]']
        return(product) #[RCH=CHR,H2]
    
    #Reaction 44 - Disproportionation (Primary Radical and Methyl Radical)
    def rxn44(self, reactants): #[RCH2',CH3]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[X4;H2:1][v3;H2:2].[X3;H3:3]>>[CH:1]=[CH2:2].[CH4:3]')
        prod = rxn.RunReactants((rlist[0], rlist[1]))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1), Chem.MolToSmiles(prod2)]
        return(product) #[RCH=CH2,CH4]
    
    #Reaction 45 - Disproportionation (Secondary Radical and Methyl Radical)
    def rxn45(self, reactants): #[RCH'R,CH3]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[X4;H2:1][v3;H1:2].[X3;H3:3]>>[CH:1]=[CH:2].[CH4:3]')
        prod = rxn.RunReactants((rlist[0], rlist[1]))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1), Chem.MolToSmiles(prod2)]
        return(product) #[RCH=CHR,CH4]
    
    #Reaction 46 - Disproportionation (Primary Radical and Primary Radical)
    def rxn46(self, reactants): #[RCH2',RCH2']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[X4;H2:1][v3;H2:2].[v3;H2:3]>>[CH:1]=[CH2:2].[CH3:3]')
        prod = rxn.RunReactants((rlist[0], rlist[1]))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1), Chem.MolToSmiles(prod2)]
        return(product) #[RCH=CH2,RCH3]
    
    #Reaction 47 - Disproportionation (Secondary Radical and Primary Radical)
    def rxn47(self, reactants): #[RCH'R,RCH2']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[X4;H2:1][v3;H1:2].[v3;H2:3]>>[CH:1]=[CH:2].[CH3:3]')
        prod = rxn.RunReactants((rlist[0], rlist[1]))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1), Chem.MolToSmiles(prod2)]
        return(product) #[RCH=CHR,RCH3]
    
    #Reaction 48 - Disproportionation (Primary Radical and Secondary Radical)
    def rxn48(self, reactants): #[RCH2',RCH'R]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[X4;H2:1][v3;H2:2].[v3;H1:3]>>[CH:1]=[CH2:2].[CH2:3]')
        prod = rxn.RunReactants((rlist[0], rlist[1]))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1), Chem.MolToSmiles(prod2)]
        return(product) #[RCH=CH2,RCH2R]
    
    #Reaction 49 - Disproportionation (Secondary Radical and Secondary Radical)
    def rxn49(self, reactants): #[RCH'R,RCH'R]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[X4;H2:1][v3;H1:2].[v3;H1:3]>>[CH:1]=[CH:2].[CH2:3]')
        prod = rxn.RunReactants((rlist[0], rlist[1]))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1), Chem.MolToSmiles(prod2)]
        return(product) #[RCH=CHR,RCH2R]
    
    #Beta Scission Reactions
    
    #Reaction 50 - Beta scission of H and Primary Radical
    def rxn50(self, reactant): #[RCH2']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[X4;H2:1][v3;H2:2]>>[CH:1]=[CH2:2]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1), '[H]']
        return(product) #[RCH=CH2,H]
    
    #Reaction 51 - Beta Scission of H and Secondary Radical
    def rxn51(self, reactant): #[RCH'R]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[X4;H2:1][v3;H1:2]>>[CH:1]=[CH:2]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1), '[H]']
        return(product) #[RCH=CHR,H]
    
    #Reaction 52 - Beta Scission of Primary Carbon and Primary Radical 
    def rxn52(self, reactant): #[CH3CH2CH2']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[v3;H2:1][X4;H2:2][X4;H3:3]>>[CH2:1]=[CH2:2].[CH3:3]')
        prod = rxn.RunReactants((rlist[0], ))
        prod1 = prod[0][0]
        prod2 = prod[0][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1), Chem.MolToSmiles(prod2)]
        return(product) #[CH2=CH2,CH3]
    
    #Reaction 53 - Beta Scission of Secondary Carbon from Primary Radical
    def rxn53(self, reactant): #[RCH2']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[v3;H2:1][X4;H2:2][X4;H2:3]>>[CH2:1]=[CH2:2].[CH2:3]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1), Chem.MolToSmiles(prod2)]
        return(product) #[RCH=CHR,H]
    
    #Reaction 54 - Beta Scission of Primary Carbon from Secondary Radical
    def rxn54(self, reactant): #[CH3CH2CH'CH3]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[X4;H3:1][v3;H1:2][X4;H2:3][X4;H3:4]>>[CH3:1][CH:2]=[CH2:3].[CH3:4]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1), Chem.MolToSmiles(prod2)]
        return(product) #[CH3CH=CH2,CH3]
    
    #Reaction 55 - Beta Scission of Secondary Carbon from Secondary Radical
    def rxn55(self, reactant): #[RCH'R]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactant]
        rxn = ReactionFromSmarts('[v3;H1:1][X4;H2:2][X4;H2:3]>>[CH:1]=[CH2:2].[CH2:3]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        prod2 = prod[idx][1]
        Chem.SanitizeMol(prod1)
        Chem.SanitizeMol(prod2)
        product = [Chem.MolToSmiles(prod1), Chem.MolToSmiles(prod2)]
        return(product) #[CH2=CHR,RCH2']
    
    #Addition Reactions
    
    #Reaction 56 - Addition of H to Unsubstituted Ethene
    def rxn56(self, reactants): #[CH2=CH2,H]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[v4;H2:1]=[v4;H2:2]>>[CH3:1][CH2:2]')
        prod = rxn.RunReactants((rlist[0], ))[0][0]
        Chem.SanitizeMol(prod)
        product = [Chem.MolToSmiles(prod)]
        return(product) #[CH3CH2']
    
    #Reaction 57 - Addition of H to Monosubstitued Alkene (Non-Substituted Carbon)
    def rxn57(self, reactants): #[RCH=CH2,H]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[v4;H1:1]=[v4;H2:2]>>[CH2:1][CH2:2]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH2CH2']
    
    #Reaction 58 - Addition of H to Monosubstituted Alkene (Substituted Carbon)
    def rxn58(self, reactants): #[RCH=CH2,H]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[v4;H1:1]=[v4;H2:2]>>[CH:1][CH3:2]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH'CH3]
    
    #Reaction 59 - Addition of H to Disubstituted Alkene
    def rxn59(self, reactants): #[RCH=CHR,H]
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[v4;H1:1]=[v4;H1:2]>>[CH2:1][CH:2]')
        prod = rxn.RunReactants((rlist[0], ))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH2CH'R]
    
    #Reaction 60 - Addition of Methyl Radical to Unsubstituted Ethene
    def rxn60(self, reactants): #[CH2=CH2,CH3']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[v4;H2:1]=[v4;H2:2].[X3;H3:3]>>[CH2:1][CH2:2][CH3:3]')
        prod = rxn.RunReactants((rlist[0],rlist[1]))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[CH3CH2CH2']
    
    #Reaction 61 - Addition of Methyl Radical to Monosubstituted Alkene
    def rxn61(self, reactants): #[RCH=CH2,CH3']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[v4;H1:1]=[v4;H2:2].[X3;H3:3]>>[CH:1][CH2:2][CH3:3]')
        prod = rxn.RunReactants((rlist[0],rlist[1]))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH'CH2CH3]
    
    #Reaction 62 - Addition of Primary Radical to Unsubstituted Ethene
    def rxn62(self, reactants): #[CH2=CH2, RCH2']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[v4;H2:1]=[v4;H2:2].[v3;H2:3]>>[CH2:1][CH2:2][CH2:3]')
        prod = rxn.RunReactants((rlist[0],rlist[1]))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH2CH2CH2']
    
    #Reaction 63 - Addition of Primary Radical to Monosubstituted Alkene
    def rxn63(self, reactants): #[RCH=CH2,RCH2']
        rlist = [Chem.MolFromSmiles(mol) for mol in reactants]
        rxn = ReactionFromSmarts('[v4;H1:1]=[v4;H2:2].[v3;H2:3]>>[CH:1][CH2:2][CH2:3]')
        prod = rxn.RunReactants((rlist[0],rlist[1]))
        if not prod:
            return ([])
        if len(prod)==1:
            idx = 0
        else:
            idx = random.randint(0,len(prod)-1)
        prod1 = prod[idx][0]
        Chem.SanitizeMol(prod1)
        product = [Chem.MolToSmiles(prod1)]
        return(product) #[RCH'CH2CH2R]