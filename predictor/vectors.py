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
