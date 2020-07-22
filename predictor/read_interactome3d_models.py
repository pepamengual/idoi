import prody
import mdtraj
from predictor.vectors import get_vector
from predictor.vectors import atom_mapping
from predictor.vectors import get_mapping
import copy

def get_protein_protein_interactions(pdb_path):
    pdb = prody.parsePDB(pdb_path)
    pdb_chains, ordered_residues = read_pdb(pdb_path)
    mapping = get_mapping()
    interface_residues = get_interface_residues(pdb, pdb_chains)
    if interface_residues:
        print(pdb_path)
        sasa_all, dssp_all = mdtraj_analysis(pdb_path)
        contact_data = get_contacts(pdb, interface_residues, mapping, ordered_residues, sasa_all, dssp_all)
        vector_data = get_vector(contact_data)
        return vector_data

def mdtraj_analysis(pdb_path):
    pdb = mdtraj.load_pdb(pdb_path)
    dssp = mdtraj.compute_dssp(pdb)
    sasa = mdtraj.shrake_rupley(pdb, mode="residue")
    return sasa[0], dssp[0]

def find_interactions(residue, chain, position, position_contacts, mapping, pdb):
    interaction_data = {}
    residue_selection = pdb.select("noh protein chain {0} and resid {1}".format(chain, position))
    distance_matrix = prody.buildDistMatrix(residue_selection, position_contacts)
    for i, (a1, r1) in enumerate(zip(residue_selection.getNames(), residue_selection.getResnames())):
        for j, (a2, r2) in enumerate(zip(position_contacts.getNames(), position_contacts.getResnames())):
            distance = distance_matrix[i][j]
            if r1 in mapping and r2 in mapping and a1 in mapping[r1] and a2 in mapping[r2]:
                a1_properties = mapping[r1][a1]
                a2_properties = mapping[r2][a2]
                
                if ("D" in a1_properties and "A" in a2_properties and distance < 3.5) or ("A" in a1_properties and "D" in a2_properties and distance < 3.5):
                    interaction_data.setdefault("h", 0)
                    interaction_data["h"] += 1 #hydrogen bonds
                if ("P" in a1_properties and "N" in a2_properties and distance < 5.0) or ("N" in a1_properties and "P" in a2_properties and distance < 5.0):
                    interaction_data.setdefault("b", 0)
                    interaction_data["b"] += 1 #salt bridge
                if ("R" in a1_properties and "R" in a2_properties and distance < 5.0):
                    interaction_data.setdefault("p", 0)
                    interaction_data["p"] += 1 #pi-pi stacking
                if ("P" in a1_properties and "A" in a2_properties and distance < 5.0) or ("A" in a1_properties and "P" in a2_properties and distance < 5.0):
                    interaction_data.setdefault("k", 0)
                    interaction_data["k"] += 1 #pi-cation
    return interaction_data


def sasa_classifier(sasa):
    if sasa < 0.25:
        return 0
    if sasa >= 0.25 and sasa < 0.50:
        return 1
    if sasa >= 0.50 and sasa < 0.75:
        return 2
    if sasa >= 0.75:
        return 3

def get_contacts(pdb, interface_residues, mapping, ordered_residues, sasa_all, dssp_all):
    contact_data = {}
    amino_acid_list = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
    for chain, position, name in interface_residues:
        position_contacts = pdb.select("(noh protein within 5.00 of (noh resid {0} and chain {1})) and not chain {1}".format(position, chain))
        ## here get interaction data between target and environment residues
        interaction_data = find_interactions(name, chain, position, position_contacts, mapping, pdb)
        if interaction_data:
            for property_ in interaction_data:
                contact_data.setdefault(chain, {}).setdefault(name, {}).setdefault(property_, 0)
                contact_data[chain][name][property_] += 1
        
        ### here get global index in the PDB for MDtraj: sasa and dssp
        mdtraj_identifier = "{}_{}_{}".format(name, chain, position)
        index_mdtraj = ordered_residues.index(mdtraj_identifier)
        sasa = sasa_all[index_mdtraj] #
        dssp = dssp_all[index_mdtraj] #
        dssp_dict = {"C": 0, "H": 1, "E": 2}
        sasa_class = sasa_classifier(sasa)
        contact_data.setdefault(chain, {}).setdefault(name, {}).setdefault("a", sasa_class)
        contact_data.setdefault(chain, {}).setdefault(name, {}).setdefault("s", dssp_dict[dssp])
        
        ## here get atom descriptors of the environment
        for atom_contact, residue_name in zip(position_contacts.getNames().tolist(), position_contacts.getResnames().tolist()):
            if residue_name in amino_acid_list:
                atom_properties = atom_mapping(residue_name, atom_contact)
                if atom_properties:
                    for property_ in atom_properties:
                        contact_data.setdefault(chain, {}).setdefault(name, {}).setdefault(property_, 0)
                        contact_data[chain][name][property_] += 1
    return contact_data

def get_interface_residues(pdb, pdb_chains):
    amino_acid_list = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
    interface_residues = []
    for chain in pdb_chains:
        interface_contacts = pdb.select("(noh protein ca same residue as within 5.00 of (noh chain {0})) and not chain {0}".format(chain))
        for chain_contact, position_contact, name_contact in zip(interface_contacts.getChids().tolist(), interface_contacts.getResnums().tolist(), interface_contacts.getResnames().tolist()):
            position_contacts = pdb.select("(noh protein within 5.00 of (noh resid {0} and chain {1})) and not chain {1}".format(position_contact, chain_contact))
            if position_contacts:
                resnames = position_contacts.getResnames().tolist()
                if position_contacts and len(resnames) > 2 and set(resnames).issubset(amino_acid_list) and name_contact in amino_acid_list:
                    interface_residues.append((chain_contact, position_contact, name_contact))
    return interface_residues

def read_pdb(pdb_path):
    pdb_chains, ordered_residues = set(), list()
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                residue = line[17:20].replace(" ", "")
                chain = line[21]
                num = line[22:26].replace(" ", "")
                name = "{}_{}_{}".format(residue, chain, num)
                pdb_chains.add(chain)
                if not name in ordered_residues:
                    ordered_residues.append(name)
    return list(pdb_chains), ordered_residues
