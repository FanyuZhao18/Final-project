# this script is to search for different conformers with ETKDG in RDKit
# the construction of short linear polymers is performed in stk (which detalis can be found in https://github.com/lukasturcani/stk)
# the surface area is calculated with Freesasa (which details can be found in https://freesasa.github.io)

from rdkit.Chem import AllChem as rdkit
import stk
import os
import freesasa
import itertools

# make conformers of the polymer
def sa_conformers(file_1, func_1, file_2, func_2, units, radius):
    # turn off cache
    stk.OPTIONS['cache'] = False
    
    # number of conformers
    N = 10
    """
    functional groups:
       ['diol'] and ['dibromine']/['difluorene']
       or
       ['bromine'] and ['bromine']/['iodine']
    """
    name_1 = file_1.replace('.mol', '')
    unit_1 = stk.StructUnit2(file_1, func_1)

    name_2 = file_2.replace('.mol', '')
    unit_2 = stk.StructUnit2(file_2, func_2)

    # make polymer
    NAME = name_1+'_'+name_2+'_AB_poly'
    print(f'Creating polymer: {NAME}')
    polymer = stk.Polymer([unit_1, unit_2], stk.Linear('AB', [0, 0], n=units, ends='h'))
    # write unoptimized structure
    polymer.write(NAME+'.mol')
    mol_polymer = rdkit.MolFromMolFile(NAME + '.mol')
    #print(f'{NAME} has {polymer.mol.get_no_atoms()} atoms!')
    print(f'Optimizing polymer {NAME} and saving {N} conformers')
    # clean molecule with ETKDG
    embedder = stk.UFF(use_cache=False)
    embedder.optimize(polymer, conformer=-1)
    # write optimized polymer to json
    polymer.dump(NAME+'_opt.json')
    polymer.write(NAME+'_opt.mol')
    # make N conformers of the polymer molecule
    etkdg = rdkit.ETKDGv2()
    etkdg.randomSeed = 1000
    etkdg.verbose = True
    etkdg.maxIterations = 200000
    cids = rdkit.EmbedMultipleConfs(
        mol=polymer.mol, 
        numConfs=N,
        params=etkdg
    )
    print(f'Made {len(cids)} conformers...')
    print(f'Warning! I have not implemented an optimization of the ETKDG cleaned polymers!')

    # iterate over conformers and save structure
    file_dir = '/home/fanyuzhao/Monomers/OH+F/dimer/conformers/'
    new_dir = file_dir+NAME+'_'+str(units)+'_'+str(radius)+'/'
    for cid in cids:
        # build directories
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
        # write optimized polymer to mol
        polymer.write(new_dir+NAME+'_'+str(cid)+'_opt.mol', conformer=cid)
        # write optimized polymer to pdb
        polymer.write(new_dir+NAME+'_'+str(cid)+'_opt.pdb', conformer=cid)
        print(f'Done! {N} ETKDG conformers of polymer written to {NAME}_{N}_opt.mol/pdb')

    # pdb file from stk can not be read in freesasa
    # save the new pdb file in rdkit from mol files
    for item in os.listdir(new_dir):
        if item.endswith('.mol'):
            file_pdb = item.replace('.mol', '')
            a = rdkit.MolFromMolFile(os.path.join(new_dir, item))
            # hydrogens are removed when converting the file in rdkit
            b = rdkit.AddHs(a, addCoords = True)
            rdkit.MolToPDBFile(b, new_dir + file_pdb + '_new.pdb')

    # calculate solvent accessible surface area(probe radius = 1.4Å and 3.6Å)
    # hydrogens are removed in the default option
    # hetatm are ignored in the default option
    options_with_Hs =  {    'hetatm' : True,
                            'hydrogen' : True,
                            'join-models' : False,
                            'skip-unknown' : False,
                            'halt-at-unknown' : False    }

    sa_list = []
    pdb_list = []
    # loop all new pdb files
    for pdb in os.listdir(new_dir):
        if pdb.endswith("_new.pdb"):
            # use freesasa to calculate SASA
            para = freesasa.Parameters()
            freesasa.Parameters.setProbeRadius(para, radius)
            free_struct = freesasa.Structure(os.path.join(new_dir, pdb), options = options_with_Hs)
            free_calc = freesasa.calc(free_struct, para)
            total = free_calc.totalArea()
            # keep 3 decimals
            decimal = round(total, 4)
            sa_list.append(decimal)
            name_pdb = pdb.replace('.pdb', '')
            pdb_list.append(name_pdb)
    # calculate average SASA(probe radius = 1.4Å)
    sa_average = round(sum(sa_list) / len(sa_list), 4)
    atom_number = mol_polymer.GetNumAtoms()
    normalized_sa = round(sa_average / atom_number, 4)
    with open (file_dir + 'Average surface area of conformers.txt', 'a+') as Asa:
        Asa.write(f'The normalized surface area of {NAME}_{units} is ' + str(normalized_sa) + ' Å^2 with the probe size of ' + str(radius) + f'Å and chain length of {units}.\n')
    print ('The avarage surface area of the conformers is ' + str(sa_average) + ' Å^2 with the probe size of ' + str(radius) + 'Å.')

    # save data to a csv table
    # save pdb file and surface area to a directory
    dic = {p: s for p, s in zip(pdb_list, sa_list)}
    download_dict = new_dir + 'Solvent accessible surface area of ' + NAME +'.csv'
    csv = open(download_dict, 'w')
    columnTitleRow = "Polymer_name, SASA\n"
    csv.write(columnTitleRow)

    for key in dic.keys():
        Polymer_name = key
        SASA = dic[key]
        row = Polymer_name + "," + str(SASA) + "\n"
        csv.write(row)
    print ('Nomalized solvent accessible surface area is '+ str(normalized_sa) + ' Å^2 with the probe size of ' + str(radius) + 'Å.')


if __name__ == "__main__":
    units = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    radii = [0, 1.4, 1.8]
    for unit, radius in itertools.product(units, radii):
        print(f'Doing unit {unit} and radius {radius}')
        sa_conformers('A.mol', ['diol'], 'B.mol', ['difluorene'], unit, radius)
