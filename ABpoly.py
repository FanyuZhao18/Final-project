# this script is to build oligomers with different chain length and orientation in stk, built in Jelfs' group
# details of stk can be found in https://github.com/lukasturcani/stk
# surface area is calculated using Freesasa: https://freesasa.github.io/

from rdkit import Chem
from rdkit.Chem import Descriptors
import stk
import os
import freesasa


def poly(file_1, func_1, file_2, func_2, units):

    # turn off cache
    stk.OPTIONS['cache'] = False

    # make polymer
    name_base =  '_diol_difluorene_poly'
    name_1 = file_1.replace('.mol', '')
    name_2 = file_2.replace('.mol', '')

    global NAME
    NAME = name_1+'_'+name_2+name_base

    unit_1 = stk.StructUnit2(file_1, func_1)
    unit_2 = stk.StructUnit2(file_2, func_2)
    polymer = stk.Polymer([unit_1, unit_2], stk.Linear('AB', [0, 0], n=units, ends='h'))
    print(f'Creating polymer: {NAME}')
    polymer.write(NAME+'.mol')
    mol_polymer = Chem.MolFromMolFile(NAME + '.mol')
    
    # optimization
    print(f'Optimizing {NAME}')
    macromodel_dir = 'pathMacroModel/'
    rff = stk.MacroModelForceField(
    macromodel_path=macromodel_dir,
    restricted=True
    )

    uff = stk.MacroModelForceField(
    macromodel_path=macromodel_dir,
    restricted=False
    )

    md = stk.MacroModelMD(
    macromodel_path=macromodel_dir,
    temperature=700,
    simulation_time=2000,
    time_step=1,
    eq_time=100
    )
    
    macromodel = stk.OptimizerSequence(rff, uff, md)
    macromodel.optimize(polymer)
    print (f'Optimization completed: {NAME}')
    
    # save files
    # make different directories
    if name_base == '_anhydride_poly':
        new_dir_1 = file_dir+'Dianhydride/'
        if not os.path.exists(new_dir_1):
            os.makedirs(new_dir_1)
        else:
        	pass
        polymer.write(new_dir_1+NAME+'.mol')
        print (f'{NAME} has been saved as dianhydride.')
        return (new_dir_1+NAME+'.mol')

    else:
        new_dir_2 = file_dir+'Polybenzodioxane/'
        if not os.path.exists(new_dir_2):
            os.makedirs(new_dir_2)
        else:
        	pass
        polymer.write(new_dir_2+NAME+'.mol')
        print (f'{NAME} has been saved as polybenzodioxane.')
        return (new_dir_2+NAME+'.mol')


def sa_calc(polymer_pdb, radius):
    # pdb files are needed for calculation surface area
    mol_file = Chem.MolFromMolFile(polymer_pdb)
    # hydrogens are removed in the mol file
    pdb_file = Chem.AddHs(mol_file, addCoords = True)
    # convert mol file to pdb file in rdkit
    Chem.MolToPDBFile(pdb_file, out_dir+NAME+'_new.pdb')

	# hydrogens are removed in the default option
    option_with_Hs =  {    'hetatm' : True,
                           'hydrogen' : True,
                           'join-models' : False,
                           'skip-unknown' : False,
                           'halt-at-unknown' : False    }

    # calculate solvent accessible surface area(probe radius = 1.4 Å or 3.6 Å)
    para = freesasa.Parameters()
    freesasa.Parameters.setProbeRadius(para, radius)
    # calculate sa for different type of polymers
    free_struct = freesasa.Structure(out_dir+NAME+'_new.pdb', options = option_with_Hs)
    free_calc = freesasa.calc(free_struct, para)
    total = free_calc.totalArea()
    # round to 4 decimals
    decimal = round(total, 4)
    print (f'Total SASA is {decimal} Å^2 when probe radius is {radius} Å.')
    atom_number = mol_file.GetNumAtoms()
    normalized_sa = round(decimal / atom_number, 4)

    # save data to a txt file
    with open (out_dir + 'Average surface area.txt', 'a+') as Asa:
       Asa.write(f'The normalized surface area of {NAME} is ' + str(normalized_sa) + ' Å^2 with the probe size of ' + str(radius) + 'Å.\n'
        )
    print ('Nomalized solvent accessible surface area is '+ str(normalized_sa) + ' Å^2 with the probe size of ' + str(radius) + 'Å.\n')

def properties(polymer_pdb):
    mol_file = Chem.MolFromMolFile(polymer_pdb)
    weight = round(Chem.Descriptors.MolWt(mol_file),4 )
    rotatable_bond = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol_file)
    bridgehead = Chem.rdMolDescriptors.CalcNumBridgeheadAtoms(mol_file)
    spiro = Chem.rdMolDescriptors.CalcNumSpiroAtoms(mol_file)
    with open (out_dir + 'Average surface area.txt', 'a+') as Asa:
       Asa.write(       f'The molecular weight of {NAME} is {weight}.\n'
                        f'The number of rotatable bonds of {NAME} is {rotatable_bond}.\n'
                        f'The number of bridgehead atoms of {NAME} is {bridgehead}.\n'
                        f'The number of spiro atoms of {NAME} is {spiro}.\n'       )

if (__name__ == "__main__"):

    file_dir = 'path1/'
    out_dir = 'path2/'
    """
    functional groups:
    ['diol'] and ['dibromine']/['difluorene']
    or
    ['bromine'] and ['bromine']/['iodine']

    name_base_is_anhydride = '_anhydride_poly'
    name_base_is_diol_difluorene = '_diol_difluorene_poly'
 
    units = chain length
    radius = probe size(1.4Å for H2 or 1.8Å for N2)
    """
    
    polymers = poly('A.mol', ['diol'], 'B.mol', ['difluorene'], 5)
    properties(polymers)

    radii = [0, 1.4, 1.8]
    for radius in radii:
        print(f'Doing radius {radius}')
        sa_calc(polymers, radius)
