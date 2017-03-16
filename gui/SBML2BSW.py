#!/usr/bin/env python
import re, libsbml, numpy
from libsbml import SBMLNamespaces
import sys, os
import operator
import numpy
import re
from numpy import savetxt
from math import isnan

sbmlns = SBMLNamespaces(3,1,"fbc",1);


class ObjectiveFunction:
        def __init__(self):
            self.type=None
            self.name=None
            self.ID=None
            self.fluxes=[]


"""
SBML to BioSimWare conversion.
TODO:   - assignment rules
        - non mass-action kinetics
        - compartments and re-labeling
        - feeds loaded from SBML
        - FBA
"""
class SBMLloader(object):

    def __init__(self, file_name):
        super(SBMLloader, self).__init__()
        self.model = None
        self.ISDET = False
        if file_name!="":
            self.reader     = libsbml.SBMLReader()
            self.sbml       = self.reader.readSBML(file_name)    
            self.model = self.sbml.getModel()            
            if self.model==None: 
                print "ERROR, empty model or file not valid:", file_name
                sys.exit(-2)
            else:
                print " * File", file_name,"loaded correctly"
        else:
            print " * Please specify a SBML file"

    def create_folder(self, OUTPATH):
        try: os.mkdir(OUTPATH)
        except: print "WARNING: directory", OUTPATH, "already exists"


    def get_initial_amounts(self, verbose=False):        
        
        ini = []

        # print " * Detected", len(self.model.getListOfSpecies()), "species"

        n = 0
        for species in self.model.getListOfSpecies(): 

            compartment = species.getCompartment()

            if n==0:
                if species.getSubstanceUnits()=="mole":
                    self.ISDET=True

        
            name = ""
            if species.getName()!="": name = species.getName()
            else:                     name = species.getId()            
            if compartment!="": 
                name = name + "_in_"+compartment

            if verbose:
                print " * Parsing species", n, name

            amount = 0
            if isnan(species.getInitialAmount()):
                amount = float(species.getInitialConcentration())
            else:
                amount = species.getInitialAmount()

            # print "Initial concentration:", species.getInitialConcentration()
            # print "Initial amount:", species.getInitialAmount()

            try:            
                if int(amount)==amount:
                    amount = int(amount)
            except:
                print "ERROR: cannot test amount of species", name
                exit(1)
            
            
            try:
                self.species_dict[name]                            
            except:
                self.species_dict[name]=n
                ini.append(amount)
                n += 1

        #print " * Found initial amounts for", len(ini), "species"
        #print " * Dictionary created for", len(self.species_dict.items()), "species"
        #print sorted(self.species_dict.items(),  key=operator.itemgetter(1))
        return ini



    def get_reactions(self, use_fba=False):

        tots_species = len(self.initial_amounts)
        self.reaction_names = []
        self.fluxes_boundaries = []
        self.reactants = []
        self.products  = []
        self.objective_functions=[]

        if use_fba:
            # FBA stuff
            fbc = self.model.getPlugin("fbc")
            print " * Detected objectives:", fbc.getNumObjectives()
            for o in xrange(fbc.getNumObjectives()):
                OF = ObjectiveFunction()            
                SBOf = fbc.getObjective(o)
                OF.name = SBOf.getName()
                OF.ID = SBOf.getId()
                OF.type = SBOf.getType()
                print " * Objective function detected:", OF.type, OF.ID
                
                # for ofl_ in range(SBOf.getNumFluxObjectives()):
                #    print ofl_
                for fo in xrange(SBOf.getNumFluxObjectives()):
                    fluxobj = SBOf.getFluxObjective(fo)
                    print " * Coefficient", fo, "for reaction", fluxobj.getReaction(), "is equal to", fluxobj.getCoefficient()
                

        # get dictionary of parameters
        PARAM_D = {}
        for p_ in range(self.model.getNumParameters()):
            P = self.model.getParameter(p_)
            pid = P.getId()
            pdict = {'id' : pid,
                     'value' : P.getValue(),
                     'constant' : P.getConstant(),
                     'sbo' : P.getSBOTermID(),
                     'name' : P.getName(),
                     'annotation': None,
                     'miriam' : None,
                     'association' : []
                     }
            PARAM_D[pid] = pdict

        """
        # get fluxes
        for r in range(self.model.getNumReactions()):
            SBRe = self.model.getReaction(r)
            RFBCplg = SBRe.getPlugin('fbc')
            lfbid = RFBCplg.getLowerFluxBound()
            ufbid = RFBCplg.getUpperFluxBound()
            #print lfbid, ufbid, PARAM_D[lfbid]['value'], PARAM_D[ufbid]['value']
        """
        
        for react in self.model.getListOfReactions():

            reactants_vector = numpy.zeros(tots_species)
            products_vector = numpy.zeros(tots_species)
            reaction_name = react.getId()
            create_reverse = False

            if reaction_name=="": reaction_name = react.getName()

            if use_fba:
                # getting fluxes
                RFBCplg = react.getPlugin('fbc')
                lfbid = RFBCplg.getLowerFluxBound()
                ufbid = RFBCplg.getUpperFluxBound()
                self.fluxes_boundaries.append( [PARAM_D[lfbid]['value'], PARAM_D[ufbid]['value']] )

            if react.getKineticLaw() != None:
            
                # gestire piu di un parametro!
                if len(react.getKineticLaw().getListOfParameters())==0:
                    #print "WARNING: cannot find local kinetic parameters"
                    nomi_parametri = [name.getId() for name in self.model.getListOfParameters() ]
                    parameters_in_kineticaw = re.findall(r"[\w']+", react.getKineticLaw().getFormula())
                    parameters_in_kineticaw = [x.strip() for x in parameters_in_kineticaw ]
                    parameters_in_kineticaw = filter( lambda x: x in nomi_parametri, parameters_in_kineticaw )
                    # print react.getName(), "PIKL:", parameters_in_kineticaw
                    if len(parameters_in_kineticaw)==0:
                        print "ERROR: can't find any kinetic parameters for reaction", reaction_name
                        exit(-3)
                    elif react.getReversible(): 
                        #print "WARNING: detected reversible reaction by getReversible", reaction_name
                        create_reverse = True                    
                    elif len(parameters_in_kineticaw)==2:
                        print "WARNING: detected two parameters in kinetic law of reaction", reaction_name, ", assuming reversible reaction"
                        create_reverse = True                    
                    elif len(parameters_in_kineticaw)==1:
                        pass
                    else:
                        print "ERROR: too many parameters in kinetic law, aborting"
                        exit(-3)

                    for el in parameters_in_kineticaw:                   

                        p = self.model.getParameter(el)

                        # print len(parameters_in_kineticaw), p
                        
                        if p.getValue()==0:
                            if not p.constant:
                                # print "WARNING: non constant parameter, assignment rule?"                            
                                if self.model.getListOfRules().get(p.getName()).isParameter():
                                    # print " * Rule for parameter", p.getName(), "detected"
                                    # print " * Rule implemented as", self.model.getListOfRules().get("k1").getFormula()
                                    tokenized_rule = self.model.getListOfRules().get("k1").getFormula()
                                    if tokenized_rule[0:8] == 'stepfunc':
                                        tokenized_rule = tokenized_rule.replace("stepfunc(", "")
                                        tokenized_rule = tokenized_rule.replace(")", "")
                                    tokenized_rule = tokenized_rule.replace(",", "")
                                    tokenized_rule =  tokenized_rule.split()
                                    temp = 0
                                    for token in tokenized_rule:                                       
                                        try:
                                            temp = float(token)
                                            if temp>0:
                                                break
                                        except: 
                                            pass
                                            #print token, "is NAN"

                                    self.PARAMS.append(temp)                
                                    #print p.getName(), float(token)
                            else:
                                #print "WARNING: constant value set to 0, parameter:", p.getName()
                                self.PARAMS.append(temp)           
                        else:
                            self.PARAMS.append(p.getValue())
                else:
                    for p in react.getKineticLaw().getListOfParameters():
                        self.PARAMS.append(p.getValue())

            # append reaction name
            self.reaction_names.append(reaction_name)                
            if create_reverse:
                self.reaction_names.append(reaction_name+" (reverse)")

            # append reactants
            temp_reactants = []
            for reactant in react.getListOfReactants():
                stoichiometry = reactant.getStoichiometry()
                alias = self.id2name[reactant.getSpecies()]            
                if alias == "": 
                    alias = reactant.getSpecies()                
                index = self.species_dict[alias]
                temp_reactants.append( str(stoichiometry) + "*" + alias )
                reactants_vector[index]=stoichiometry
            self.reactants.append( "+".join(temp_reactants) )

            # append products
            temp_products = []
            for product in react.getListOfProducts():
                alias = self.id2name[product.getSpecies()]      
                if alias == "": 
                    alias = reactant.getSpecies()
                stoichiometry = product.getStoichiometry()
                index = self.species_dict[alias]
                temp_products.append( str(stoichiometry) + "*" + alias)
                products_vector[index]=stoichiometry
            self.products.append( "+".join(temp_products) )

            self.LEFT.append(reactants_vector)
            self.RIGHT.append(products_vector)

            # create reverse reaction
            if create_reverse:
                self.LEFT.append(products_vector)
                self.RIGHT.append(reactants_vector)
                self.products.append( "+".join(temp_reactants) )
                self.reactants.append( "+".join(temp_products) )

            #exit()




    def convert_to_biosimware(self, OUTPATH, verbose=False):

        self.create_folder(OUTPATH)
        self.process_sbml()

        tots_species = len(self.initial_amounts)
        feeds = numpy.zeros(tots_species)

        cwd = os.getcwd()
        os.chdir(OUTPATH)

        if verbose: print " * Creating 'M_0' with initial state...",
        savetxt("M_0", [self.initial_amounts], fmt="%e", delimiter="\t")
        if verbose: print "DONE"
        
        if verbose: print " * Creating 'M_feed' file...",        
        for n,species in enumerate(self.model.getListOfSpecies()):
            if species.constant: feeds[n]=1
        numpy.savetxt("M_feed", [feeds], fmt="%d", delimiter="\t")
        if verbose: print "DONE"

        if verbose: print " * Creating 'alphabet' with names of species...", 
        with open("alphabet", "w") as fo:
            for (k,v) in self.sorted_species_dict: fo.write(k+"\t")
        if verbose: print "DONE"
        
        if verbose: print " * Creating 'left_side' matrix with reactants",
        numpy.savetxt("left_side", self.LEFT, fmt="%d", delimiter="\t")
        if verbose: print "DONE"; print " * Creating 'right_side' matrix with products...",
        numpy.savetxt("right_side", self.RIGHT, fmt="%d", delimiter="\t")
        if verbose: print "DONE"; print " * Creating 'c_vector' file for parameters...",
        numpy.savetxt("c_vector", numpy.array([self.PARAMS]), fmt="%e", delimiter="\n")
        if verbose: print " * Creating 'boundaries' matrix for FBA fluxes limits...", 
        numpy.savetxt("boundaries", self.fluxes_boundaries, fmt="%e", delimiter="\t")
        if verbose: print "DONE"

        os.chdir(cwd)
        

    def process_sbml(self):

        # create dictionary of species IDs
        self.id2name = {}
        for species in self.model.getListOfSpecies():            
            fullname = species.getName()
            if fullname == "":
                fullname = species.getId()
            if species.getCompartment()!="":
                fullname = fullname + "_in_"+species.getCompartment()
            #print fullname
            self.id2name[species.getId()] = fullname

        # create reverse dictionary of species indices
        self.species_dict = {}        
        
        # container of initial species amounts (sorted)
        self.initial_amounts = self.get_initial_amounts()
        tots_species = len(self.initial_amounts)
        
        self.sorted_species_dict = sorted(self.species_dict.items(),  key=operator.itemgetter(1))
       
        self.LEFT, self.RIGHT, self.PARAMS = [], [], []
        self.get_reactions()    


def separator():
    print
    print "*"*100
    print

if __name__ == '__main__':

    INPUT_FILE = "Breast_core.xml"
    OUTPUT_FOLDER = "./output_marzia"

    if len(sys.argv)>1: INPUT_FILE = argv[1]
    if len(sys.argv)>2: INPUT_FILE = argv[2]

    SL = SBMLloader(INPUT_FILE)  
    SL.convert_to_biosimware(OUTPUT_FOLDER)
    print "Reaction names", SL.reaction_names
    separator()
    print "Reaction parameters", SL.PARAMS
    separator()
    print "PSA kinetic parameters", map(lambda x: "K_"+x, SL.reaction_names)
    separator()
    print "PSA chemical species", map(lambda x: x[0], SL.sorted_species_dict)
    separator()
    print "FBA boundaries", SL.fluxes_boundaries
    separator()
    print "Reactants", SL.reactants
    separator()
    print "Products", SL.products
    separator()
    