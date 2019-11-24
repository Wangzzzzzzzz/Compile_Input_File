from Bio.PDB import PDBParser, PDBList
from Bio.PDB.Polypeptide import three_to_one, PPBuilder
from Bio import SeqIO
import pandas as pd
import numpy as np
import os, warnings
import time

CACHE = {'pdb_id':None,'Model':None,'PP_Chain':None,'Num_Uni':None, 'Uni_Chain':None}
FLAG_TABL = {0:'Normal',1:'2 unique seq with duplicates',
             2: 'More than 2 unique seq', 3: 'Less than 2 seq'}


def read_in_file(file_dir):
    file_dir = './dataset/' + file_dir
    dataset = pd.read_csv(file_dir,sep='\t')
    return dataset
    
def obtain_pdb(dataset):
    # download the pdb set
    pdb_set = list(set(dataset['pdb']))
    if not os.path.exists('./PDB'):
        os.mkdir('./PDB')
    PDB_dl_handle = PDBList(verbose=False)
    PDB_dl_handle.download_pdb_files(pdb_codes=pdb_set,file_format='pdb',pdir='./PDB')
    # download the fasta set
    if not os.path.exists('./FASTA'):
        os.mkdir('./FASTA')
    for item in pdb_set:
        if not os.path.exists("./FASTA/{}.fasta".format(item)):
            os.system("""wget -q -O - "https://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList={}&compressionType=uncompressed" > ./FASTA/{}.fasta""".format(item,item))
            time.sleep(0.1)
        assert (os.stat("./FASTA/{}.fasta".format(item)).st_size != 0), 'Download Failed, Empty File generated.'

# find the unique chains 
def find_unique_seq(mutation):
    global CACHE
    Chains = CACHE['PP_Chain']
    mut_chain_id = mutation[1]
    # first add in the mutation chain id
    seq_set = {Chains[mut_chain_id]}
    unique_id = [mut_chain_id]
    # add any other different seq
    for chain_id, seq in Chains.items():
        if not seq in seq_set:
            unique_id.append(chain_id)
            seq_set.add(seq)
    # number of unique seq
    num_unique = len(unique_id)
    unique_chain = {uni_id:Chains[uni_id] for uni_id in unique_id}

    return num_unique, unique_chain

def output_mutant_seq(raw_dataset, output_file):
    global CACHE
    chain_tracker = set()
    chain_store = dict()
    output_file = './output/' + output_file
    for _, raw_series in raw_dataset.iterrows():
        pdb_id = raw_series['pdb']
        mutation = raw_series['mutation']
        WD = mutation[0]
        MU = mutation[-1]
        chain_pos = mutation[2:-1]
        mu_chain_id = mutation[1]
        if chain_pos[-1].isdigit():
            mu_residue_id = (" ", int(chain_pos), " ")
        else:
            mu_residue_id = (" ", int(chain_pos[:-1]), chain_pos[-1].upper())
        parser = PDBParser()
        builder = PPBuilder()
        if not CACHE['pdb_id'] == pdb_id:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                pdb_struc = parser.get_structure(id=pdb_id,file='./PDB/pdb'+pdb_id.lower()+'.ent')
            # store to CACHE if encounter a new pdb
            model = pdb_struc[0]
            CACHE['pdb_id'], CACHE['Model'], CACHE['Num_Uni'], CACHE['Uni_Chain'] = pdb_id, model, None, None
            with open('./FASTA/' + pdb_id + '.fasta','r') as hdl:
                CACHE['PP_Chain'] = {rec.id[5]:str(rec.seq) for rec in SeqIO.parse(hdl,"fasta")}
        else:
            model = CACHE['Model']

        num_unique_chains = CACHE['Num_Uni']
        unique_chains = CACHE['Uni_Chain']
        if num_unique_chains is None:
            num_unique_chains, unique_chains = find_unique_seq(mutation)
            CACHE['Num_Uni'] = num_unique_chains
            CACHE['Uni_Chain'] = unique_chains

        # handle wild type
        for unique_chain_id in unique_chains:
            protein_chain = model[unique_chain_id]
            seq_identifier = pdb_id + '_' + unique_chain_id
            if not seq_identifier in chain_tracker:
                chain = "".join([str(pp.get_sequence()).replace('\n', '')
                                 for pp in builder.build_peptides(protein_chain)])
                chain_store[seq_identifier] = chain
                chain_tracker.add(seq_identifier)
        
        # find the mutant chain
        protein_chain = model[mu_chain_id]
        seq_identifier = pdb_id + '_' + mu_chain_id + '_' + mutation
        if not seq_identifier in chain_tracker:
            mu_chain = []
            for residue in protein_chain:
                if residue.get_id()[0] != ' ':
                    continue
                else:
                    if residue.get_id()[1] == mu_residue_id[1] and residue.get_id()[2] == mu_residue_id[2]:
                        assert three_to_one(residue.get_resname()) == WD
                        mu_chain.append(MU)
                    else:
                        try:
                            mu_chain.append(three_to_one(residue.get_resname()))
                        except:
                            pass
            chain = "".join(mu_chain)
            chain_store[seq_identifier] = chain
            chain_tracker.add(seq_identifier)

    with open(output_file, 'w') as output_hl:
        for k, v in chain_store.items():
            output_hl.write(k+'\t'+v+'\n')

# handle cases where there are more than 2 unique chaines
def multi_unique_chain_handler(mutation, ddG, unique_chains):
    global CACHE
    model = CACHE['Model']
    pdb_id = CACHE['pdb_id']
    # take the unique chains as representative
    # also gather all of other chains as one group
    WD = mutation[0]
    MU = mutation[-1]
    chain_pos = mutation[2:-1]
    mu_chain_id = mutation[1]
    if chain_pos[-1].isdigit():
        residue_id = (" ", int(chain_pos), " ")
    else:
        residue_id = (" ", int(chain_pos[:-1]), chain_pos[-1].upper())
    assert three_to_one(model[mu_chain_id][residue_id].get_resname()) == WD, 'The mutation position is incorrect.'
    Wild = {'Mutation':"", 'Wild':""}
    Mut = {'Mutation': "", 'Wild': ""}
    for item in unique_chains.keys():
        chain = model[item]
        if chain.id == mu_chain_id:
            Wild['Mutation'] += (pdb_id + '_' + chain.id)
            Mut['Mutation'] += (pdb_id + '_' + chain.id + '_' + mutation)
        else:
            Wild['Wild'] += (pdb_id + '_' + chain.id + ', ')
            Mut['Wild'] += (pdb_id + '_' + chain.id + ', ')
    return (Wild['Mutation'], Wild['Wild'][:-2], Mut['Mutation'], Mut['Wild'][:-2], ddG, 2)


# handle cases where there are more than 2 chains in the protein
def multi_chain_handler(mutation, ddG):
    global CACHE
    model = CACHE['Model']
    pdb_id = CACHE['pdb_id']
    num_unique_chains = CACHE['Num_Uni']
    unique_chains = CACHE['Uni_Chain']
    if num_unique_chains is None:
        num_unique_chains, unique_chains = find_unique_seq(mutation)
        CACHE['Num_Uni'] = num_unique_chains
        CACHE['Uni_Chain'] = unique_chains
    if num_unique_chains > 2:
        return multi_unique_chain_handler(mutation, ddG, unique_chains)
    else:
        # take the unique chains as representative
        WD = mutation[0]
        MU = mutation[-1]
        chain_pos = mutation[2:-1]
        mu_chain_id = mutation[1]
        if chain_pos[-1].isdigit():
            residue_id = (" ", int(chain_pos), " ")
        else:
            residue_id = (" ", int(chain_pos[:-1]),chain_pos[-1].upper())
        assert three_to_one(model[mu_chain_id][residue_id].get_resname()) == WD, 'The mutation position is incorrect.'
        Wild = []
        Mut = []
        for item in unique_chains.keys():
            chain = model[item]
            Wild.append(pdb_id + '_' + chain.id)
            if chain.id == mu_chain_id:
                Mut.append(pdb_id + '_' + chain.id + '_' + mutation)
            else:
                Mut.append(pdb_id + '_' + chain.id)
        return (Wild[0], Wild[1], Mut[0], Mut[1], ddG, 1)
    
def generate_mut_row(raw_series):
    global CACHE
    pdb_id = raw_series['pdb']
    mutation = raw_series['mutation']
    ddG = raw_series['actual']
    parser = PDBParser()
    if not CACHE['pdb_id'] == pdb_id:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pdb_struc = parser.get_structure(id=pdb_id,file='./PDB/pdb'+pdb_id.lower()+'.ent')
        # store to CACHE if encounter a new pdb
        model = pdb_struc[0]
        CACHE['pdb_id'], CACHE['Model'], CACHE['Num_Uni'], CACHE['Uni_Chain'] = pdb_id, model, None, None
        with open('./FASTA/' + pdb_id + '.fasta','r') as hdl:
            CACHE['PP_Chain'] = {rec.id[5]:str(rec.seq) for rec in SeqIO.parse(hdl,"fasta")}
    else:
        model = CACHE['Model']

    if len(CACHE['PP_Chain'].keys()) > 2:
        return multi_chain_handler(mutation, ddG)
    elif len(CACHE['PP_Chain'].keys()) < 2:
        UserWarning(pdb_id + 'Protein has less than two chains and will be omitted')
        return (CACHE['pdb_id'], np.nan, np.nan, np.nan, np.nan, 3)
    else:
        WD = mutation[0]
        MU = mutation[-1]
        chain_pos = mutation[2:-1]
        mu_chain_id = mutation[1]
        if chain_pos[-1].isdigit():
            residue_id = (" ", int(chain_pos), " ")
        else:
            residue_id = (" ", int(chain_pos[:-1]),chain_pos[-1].upper())
        assert three_to_one(model[mu_chain_id][residue_id].get_resname()) == WD, 'The mutation position is incorrect.'
        Wild = []
        Mut = []
        for chain in model:
            Wild.append(pdb_id + '_' + chain.id)
            if chain.id == mu_chain_id:
                Mut.append(pdb_id + '_' + chain.id + '_' + mutation)
            else:
                Mut.append(pdb_id + '_' + chain.id)
        return (Wild[0], Wild[1], Mut[0], Mut[1], ddG, 0)

def summarize_flag(raw_dataset, output):
    # generate a new dataset with protein chain
    flag_dict = dict()
    protein_tracker = set()
    for _, row in raw_dataset.iterrows():
        pdb_id = row['Wild1'][0:4]
        if pdb_id in protein_tracker:
            assert row['Flag'] == flag_dict[pdb_id], 'One Protein will not have more than 1 flag.'
        else:
            flag_dict[pdb_id] = row['Flag']
            protein_tracker.add(pdb_id)
    protein_flag_df = pd.DataFrame(list(flag_dict.items()),columns=['Protein','Flag'])
    protein_flag_df['Flag_Meaning'] = [FLAG_TABL[i] for i in protein_flag_df['Flag']]
    flag_summary = protein_flag_df['Flag'].value_counts(sort=False)
    protein_flag_summary = pd.DataFrame({
        'Flag_Label': flag_summary.index,
        'Count': flag_summary.values,
        'Flag_Meaning': [FLAG_TABL[i] for i in flag_summary.index]
    })
    with pd.ExcelWriter(output) as output_wr:
        protein_flag_summary.to_excel(output_wr, 'Summary')
        protein_flag_df.to_excel(output_wr,'Protein_Flag')

def main():
    if not os.path.exists('./dataset'):
        os.mkdir('./dataset')

    if not os.path.exists('./output'):
        os.mkdir('./output')
    
    dataset = read_in_file('training_10t10f_S4169.csv')
    obtain_pdb(dataset)
    output = []
    for _, row in dataset.iterrows():
        output.append(generate_mut_row(row))

    raw_data = pd.DataFrame(output, columns=['Wild1', 'Wild2', 'Mut1', 'Mut2', 'ddG', 'Flag'])
    raw_data.to_csv('./output/raw.csv',sep='\t',index=False)
    formatted_data = raw_data[['Wild1', 'Wild2', 'Mut1', 'Mut2', 'ddG']]
    formatted_data.to_csv('./output/formatted.csv',sep='\t', index=False, header=False)
    Res_summary = raw_data['Flag'].value_counts(sort=False)
    print(Res_summary)
    Summary = pd.DataFrame({'Flag_Label':Res_summary.index, 
                            'Count':Res_summary.values,
                            'Flag_Meaning':[FLAG_TABL[i] for i in Res_summary.index]})
    Summary.to_csv('./output/Summary.csv', sep='\t', index=False)

    # output summary for each protein
    summarize_flag(raw_data,'./output/Protein_summary.xlsx')

    # output mutant chains
    output_mutant_seq(dataset, 'MU_WD_seq.txt')
    


if __name__ == "__main__":
    main()
