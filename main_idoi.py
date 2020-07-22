from predictor.read_interactome3d_models import get_protein_protein_interactions
from multiprocessing import Pool, TimeoutError
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.svm import SVC
from predictor.get_interactome3d_models import select_interactome_complexes
import glob
import os
import pickle


def make_nn(df):
    from keras.models import Model, Sequential
    from keras.layers import Input, Dense, Dropout, Flatten, Activation, Reshape
    from keras.callbacks import EarlyStopping
    import tensorflow as tf
    from keras import regularizers
    from keras import backend as K
    from sklearn.preprocessing import OneHotEncoder, StandardScaler, LabelEncoder
    from keras.utils import np_utils

    column_names = list(df.columns.values)
    column_names.remove("Class")
    X = df[column_names].to_numpy()
    y = df[["Class"]].to_numpy()
    encoder = LabelEncoder()
    encoder.fit(y)
    encoded_y = encoder.transform(y)
    y = np_utils.to_categorical(encoded_y)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state = 20, test_size=0.25, shuffle=True)

    neurons = len(list(df.drop(['Class'], axis=1)))
    b_size = 100
    model = Sequential()
    model.add(Dense(500, input_dim=neurons, activation='relu'))#, kernel_initializer='he_normal', kernel_regularizer=regularizers.l2(0.01)))
    model.add(Dense(250, activation='relu'))
    model.add(Dense(100, activation='relu'))
    #model.add(Dropout(0.3))
    model.add(Dense(20, activation='softmax'))
    model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
    es = EarlyStopping(monitor='val_loss', mode='min', patience=5, verbose=1)
    #history = model.fit(X_train, y_train, epochs=100, batch_size=b_size, validation_data=(X_test, y_test), callbacks=[es], verbose=1)
    history = model.fit(X_train, y_train, epochs=20, batch_size=b_size, validation_data=(X_test, y_test), verbose=1)

def make_svm(df):
    column_names = list(df.columns.values)
    column_names.remove("Class")
    X = df[column_names].to_numpy()
    y = df[["Class"]].to_numpy()
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state = 20, test_size=0.90)
    print("Train and test splitted")
    svm_model_linear = SVC(kernel = 'linear', C = 1).fit(X_train, y_train)
    print("Model finished...")

    filename = "SVC_linear_model.sav"
    pickle.dump(svm_model_linear, open(filename, "wb"))
    #svm_model_linear = pickle.load(open(filename, "rb"))
    svm_predictions = svm_model_linear.predict(X_test)
    accuracy = svm_model_linear.score(X_test, y_test)
    cm = confusion_matrix(y_test, svm_predictions)
    print(accuracy)
    print(cm)

def export_data(all_vector_data, export_name):
    df = pd.DataFrame(all_vector_data)
    df.to_csv(export_name)


def generate_data(export_name, selection, path_folders):
    all_vector_data = []
    pool = Pool(20)
    multiple_results = []

    for folder in glob.glob(path_folders):
        pdb_list = glob.glob(os.path.join(folder, "*.pdb"))
        for pdb in pdb_list:
            if pdb.split("/")[-1] in selection:
                #vector_data = get_protein_protein_interactions(pdb)
                multiple_results.append(pool.apply_async(get_protein_protein_interactions, (pdb,)))
    print("Getting results")
    for result in multiple_results:
        try:
            vector_data = result.get(timeout=30)
            if vector_data:
                #print(vector_data)
                all_vector_data.extend(vector_data)
        except Exception as e:
            print(e)
            continue
        except TimeoutError:
            continue
    pool.close()
    pool.terminate()
    return all_vector_data


def main():
    all_vector_data = []
    export_name = "data_contacts_interactome3d.csv"
    selection, path_folders = select_interactome_complexes()
    print("Interactome3D complexes selected")
    all_vector_data = generate_data(export_name, selection, path_folders)
    export_data(all_vector_data, export_name)
    
    df = pd.read_csv(export_name, index_col=0)
    print(df)
    make_nn(df)

main()




