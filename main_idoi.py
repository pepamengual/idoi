from predictor.read_interactome3d_models import get_protein_protein_interactions
from multiprocessing import Pool
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.svm import SVC

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
    all_vector_data = []
    export_name = "data_contacts_interactome3d.csv"
    pdb_list = ["../1ao7.pdb"]
    
    pool = Pool(2) 
    multiple_results = []
    for pdb in pdb_list:
        multiple_results.append(pool.apply_async(get_protein_protein_interactions, (pdb,)))
    for result in multiple_results:
        vector_data = result.get()
        all_vector_data.extend(vector_data)
    
    export_data(all_vector_data, export_name)
    
    df = pd.read_csv(export_name, index_col=0)
    print(df)
    make_svm(df)

main()













