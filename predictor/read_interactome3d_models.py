import prody
import predictor.vectors as vectors

def get_protein_protein_interactions(pdb_path):
    pdb = prody.parsePDB(pdb_path)
    pdb_chains = read_pdb_chains(pdb_path)
    interface_residues = get_interface_residues(pdb, pdb_chains)
    contact_data = get_contacts(pdb, interface_residues)
    vector_data = vectors.get_vector(contact_data)
    return vector_data

def get_contacts(pdb, interface_residues):
    contact_data = {}
    for chain, position, name in interface_residues:
        position_contacts = pdb.select("(within 5.00 of (noh resid {0} and chain {1})) and not chain {1}".format(position, chain))
        if position_contacts:
            for atom_contact, residue_name in zip(position_contacts.getNames().tolist(), position_contacts.getResnames().tolist()):
                if not atom_contact.startswith("H"):
                    properties = vectors.atom_mapping(residue_name, atom_contact)
                    for property_ in properties:
                        contact_data.setdefault(chain, {}).setdefault(name, {}).setdefault(property_, 0)
                        contact_data[chain][name][property_] += 1
    return contact_data

    vector_data = get_vector(data)

def get_interface_residues(pdb, pdb_chains):
    interface_residues = []
    for chain in pdb_chains:
        interface_contacts = pdb.select("(ca same residue as within 5.00 of (noh chain {0})) and not chain {0}".format(chain))
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
