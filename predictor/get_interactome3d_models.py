import pandas as pd

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

