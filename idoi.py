from multiprocessing import Pool
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.svm import SVC
import glob
import os
import prody

def filter_by_sequence_identity(output_name, identity):
    df = pd.read_csv(output_name, index_col=0)
    df_short = df[df < identity].dropna()
    return df, df_short

def get_sequences_identity(df_short):
    selection = df_short.index.tolist()
    selection_split = set([x.split("_")[0] for i in selection for x in i.split(";")])
    return selection_split

def select_interactome_complexes():
    output_name = "/home/pepamengual/uep_mutations/alignment/merged_alignments.csv"
    identity = 30
    df, df_short = filter_by_sequence_identity(output_name, identity)
    selection = get_sequences_identity(df_short)
    path_folders="/home/pepamengual/UEPPi/ueppi_script/training/all_complexes/interactome_*"
    return selection, path_folders




def get_protein_protein_interactions(pdb_path):
    pdb = prody.parsePDB(pdb_path)
    pdb_chains = read_pdb_chains(pdb_path)
    interface_residues = get_interface_residues(pdb, pdb_chains)
    contact_data = get_contacts(pdb, interface_residues)
    vector_data = get_vector(contact_data)
    return vector_data

def get_contacts(pdb, interface_residues):
    letters = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
    contact_data = {}
    for chain, position, name in interface_residues:
        position_contacts = pdb.select("(protein within 5.00 of (noh resid {0} and chain {1})) and not chain {1}".format(position, chain))
        if position_contacts:
            for atom_contact, residue_name in zip(position_contacts.getNames().tolist(), position_contacts.getResnames().tolist()):
                if not atom_contact.startswith("H") and residue_name in letters:
                    properties = atom_mapping(residue_name, atom_contact)
                    if properties:
                        for property_ in properties:
                            contact_data.setdefault(chain, {}).setdefault(name, {}).setdefault(property_, 0)
                            contact_data[chain][name][property_] += 1
    return contact_data

    vector_data = get_vector(data)

def get_interface_residues(pdb, pdb_chains):
    interface_residues = []
    for chain in pdb_chains:
        interface_contacts = pdb.select("(protein ca same residue as within 5.00 of (noh chain {0})) and not chain {0}".format(chain))
        for chain_contact, position_contact, name_contact in zip(interface_contacts.getChids().tolist(), interface_contacts.getResnums().tolist(), interface_contacts.getResnames().tolist()):
            interface_residues.append((chain_contact, position_contact, name_contact))
    return interface_residues

def read_pdb_chains(pdb_path):
    pdb_chains = set()
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21]
                pdb_chains.add(chain)
    return list(pdb_chains)




def get_vector(data):
    aa_code = {}
    letters = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
    for i, letter in enumerate(letters):
        aa_code.setdefault(letter, i)

    vector_data = []
    for chain, position_dict in data.items():
        for position, property_dict in position_dict.items():
            amino_acid = position.split("_")[0]
            vector = {"P": 0, "N": 0, "A": 0, "D": 0, "H": 0, "R": 0, "S": 0, "Class": 0}
            if amino_acid in letters:
                for property_, value in property_dict.items():
                    vector[property_] += value
                vector["Class"] += aa_code[amino_acid]
            vector_data.append(vector)
    return vector_data

def atom_mapping(residue_name, atom_contact):
    #https://cdn.rcsb.org/wwpdb/docs/documentation/file-format/PDB_format_1992.pdf
    #http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/charge/
    #POSITIVE: P, NEGATIVE: N, HB ACCEPTOR: A,
    #HB DONOR: D, HYDROPHOBIC: H, AROMATIC: R, SULPHUR: S
    mapping = {"ALA": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H"},
               "ARG": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "CG": "H", "CD": "H", "NE": "PD", "CZ": "H", "NH1": "PD", "NH2": "PD"},
               "ASN": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "CG": "H", "OD1": "A", "ND2": "D"},
               "ASP": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "CG": "H", "OD1": "NA", "OD2": "NA"},
               "CYS": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "SG": "S"},
               "GLU": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "CG": "H", "CD": "H", "OE1": "NA", "OE2": "NA"},
               "GLN": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "CG": "H", "CD": "H", "OE1": "A", "NE2": "D"},
               "GLY": {"N": "D", "CA": "H", "C": "H", "O": "A"},
               "HIS": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "CG": "H", "ND1": "DA", "CD2": "H", "CE1": "H", "NE2": "DA"},
               "ILE": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "CG1": "H", "CG2": "H", "CD1": "H"},
               "LEU": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "CG": "H", "CD1": "H", "CD2": "H"},
               "LYS": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "CG": "H", "CD": "H", "CE": "H", "NZ": "PD"},
               "MET": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "CG": "H", "SD": "S", "CE": "H"},
               "PHE": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "CG": "R", "CD1": "R", "CD2": "R", "CE1": "R", "CE2": "R", "CZ": "R"},
               "PRO": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "CG": "H", "CD": "H"},
               "SER": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "OG": "AD"},
               "THR": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "OG1": "AD", "CG2": "H"},
               "TRP": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "CG": "R", "CD1": "R", "NE1": "D", "CE2": "R", "CZ2": "R", "CH2": "R", "CZ3": "R", "CE3": "R", "CD2": "R"},
               "TYR": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "CG": "R", "CD1": "R", "CE1": "R", "CZ": "R", "OH": "DA", "CE2": "R", "CD2": "R"},
               "VAL": {"N": "D", "CA": "H", "C": "H", "O": "A", "CB": "H", "CG1": "H", "CG2": "H"}}
    if residue_name in mapping and atom_contact in mapping[residue_name]:
        return mapping[residue_name][atom_contact]
    else:
        return None

def make_svm(df):
    column_names = list(df.columns.values)
    column_names.remove("Class")
    X = df[column_names].to_numpy()
    y = df[["Class"]].to_numpy()
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state = 0)
    svm_model_linear = SVC(kernel = 'linear', C = 1).fit(X_train, y_train)

    svm_predictions = svm_model_linear.predict(X_test)
    accuracy = svm_model_linear.score(X_test, y_test)
    cm = confusion_matrix(y_test, svm_predictions)
    print(accuracy)
    print(cm)

def export_data(all_vector_data, export_name):
    df = pd.DataFrame(all_vector_data)
    df.to_csv(export_name)

def main():
    selection, path_folders = select_interactome_complexes()

    all_vector_data = []
    export_name = "data_contacts_interactome3d.csv"
    pool = Pool(20) 
    multiple_results = []
    
    for folder in glob.glob(path_folders):
        pdb_list = glob.glob(os.path.join(folder, "*.pdb"))
        for pdb in pdb_list:
            if pdb.split("/")[-1] in selection:
                multiple_results.append(pool.apply_async(get_protein_protein_interactions, (pdb,)))
    for result in multiple_results:
        try:
            vector_data = result.get()
            print(vector_data)
            all_vector_data.extend(vector_data)
        except:
            continue
    export_data(all_vector_data, export_name)
    
    df = pd.read_csv(export_name, index_col=0)
    make_svm(df)

main()




