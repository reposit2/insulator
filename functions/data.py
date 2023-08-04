# This code was created based on codes in 'Tasaki, S., Gaiteri, C., Mostafavi, S. & Wang, Y. Deep learning decodes the principles of differential gene expression. Nature Machine Intelligence (2020)'

def batch_iter2(
               X_promoter_train,
               Y_train,
               batch_size,
               shuffle=True):

    import numpy as np

    # X_promoter_train: promoter features
    # Y_train: target data
    # batch_size: The number of data in each batch
    # shuffle=True: Shuffling the order of data in each epoch

    # Number of data in each batch
    data_size = len(Y_train)
    num_batches_per_epoch = int((data_size - 1) / batch_size) + 1

    # Obtain size of promoter feature matrix
    n_feature_promoter=X_promoter_train[0].shape[0]

    def data_generator():
        while True:

            # Shuffle the data at each epoch
            if shuffle:
                shuffle_indices = np.random.permutation(np.arange(data_size))
                shuffled_data_promoter = X_promoter_train[shuffle_indices]
                shuffled_labels = Y_train[shuffle_indices]
            else:
                shuffled_data_promoter = X_promoter_train
                shuffled_labels = Y_train

            # Generate data for each batch
            for batch_num in range(num_batches_per_epoch):
                start_index = batch_num * batch_size
                end_index = min((batch_num + 1) * batch_size, data_size)

                # Prepare promoter feature matrix
                # Obtain max length of promoter sequence in this batch
                seq_len_promoter=[np.shape(x)[1] for x in shuffled_data_promoter[start_index: end_index]]
                max_seq_promoter=np.max(seq_len_promoter)

                # Initialize promoter feature matrix
                X_promoter=np.zeros((end_index-start_index, max_seq_promoter, n_feature_promoter))
                k=0
                for i in range(start_index,end_index):
                    X_promoter[k,0:seq_len_promoter[k],:]=shuffled_data_promoter[i].transpose()
                    k=k+1
                X_promoter = np.asarray(X_promoter, dtype=np.float32)

                # Prepare target data
                Y = shuffled_labels[start_index: end_index]
                Y = np.asarray(Y, dtype=np.float32)

                yield X_promoter, Y

    return num_batches_per_epoch, data_generator()


def batch_iter_DeepLIFT2(
                        X_promoter_train,
                        Y_train,
                        batch_size,
                        med_promoter_len,
                        shuffle=True):
    import numpy as np

    # X_promoter_train: promoter features
    # Y_train: target data
    # batch_size 
    # shuffle=True

    # Number of data in each batch
    data_size = len(Y_train)
    num_batches_per_epoch = int((data_size - 1) / batch_size) + 1

    # Obtain the size of promotor feature matrix
    n_feature_promoter=X_promoter_train[0].shape[0]

    def data_generator():
        while True:
            # Shuffle the data at each epoch
            if shuffle:
                shuffle_indices = np.random.permutation(np.arange(data_size))
                shuffled_data_promoter = X_promoter_train[shuffle_indices]
                shuffled_labels = Y_train[shuffle_indices]
            else:
                shuffled_data_promoter = X_promoter_train
                shuffled_labels = Y_train

            # Generate data for each batch
            for batch_num in range(num_batches_per_epoch):
                start_index = batch_num * batch_size
                end_index = min((batch_num + 1) * batch_size, data_size)

                # Prepare promoter feature matrix
                # Obtain max length of promoter sequence in this batch
                seq_len_promoter=[np.shape(x)[1] for x in shuffled_data_promoter[start_index: end_index]]
                max_seq_promoter=np.max(seq_len_promoter)
                max_seq_promoter=np.max([max_seq_promoter,med_promoter_len])

                # Initializing promoter feature matrix
                X_promoter=np.zeros((end_index-start_index, max_seq_promoter, n_feature_promoter))
                k=0
                for i in range(start_index,end_index):
                    X_promoter[k,0:seq_len_promoter[k],:]=shuffled_data_promoter[i].transpose()
                    k=k+1

                # Prepare target data
                Y = shuffled_labels[start_index: end_index]

                yield X_promoter, Y

    return num_batches_per_epoch, data_generator()


def load_ml_data(deg_data_file,
                 mRNA_data_loc,
                 mRNA_annotation_data,
                 promoter_data_loc,
                 promoter_annotation_data):
    
    # deg_data_file: a full path to transcriptome data
    # mRNA_data_loc: a full path to a directory where mRNA (promoter in this study) annotation data is stored
    # mRNA_annotation_data: a file name of mRNA (promoter in this study) annotation data
    # promoter_data_loc:  a full path to a directory where promoter (enhancer in this study) annotation data is stored
    # promoter_annotation_data: a file name of promoter (enhancer in this study) annotation data
                    
    import pandas as pd
    import numpy as np
    from scipy.sparse import vstack
    
    # Subset DEG data
    import os.path
    root, ext = os.path.splitext(deg_data_file)
    if ext==".txt":
        deg_data = pd.read_csv(deg_data_file,sep="\t")
        if deg_data.columns[0] != 'Name':
            print('The first column must be Name')
            sys.exit(1) 

        med_exp = np.median(deg_data.values[:,1:],axis=1)
        for i in range(deg_data.shape[0]):
            deg_data.iloc[i,1:]  = deg_data.values[i,1:] - med_exp[i]
        deg_data['MedianExp'] = med_exp
    elif ext==".gz":
        deg_data = pd.read_csv(deg_data_file,sep="\t",compression='gzip')
        if deg_data.columns[0] != 'Name':
            print('The first column must be Name')
            sys.exit(1) 

        med_exp = np.median(deg_data.values[:,1:],axis=1)
        deg_data['MedianExp'] = med_exp
        print('deg_data', deg_data) 
    elif ext==".pkl":
        deg_data = pd.read_pickle(deg_data_file)
        if deg_data.columns[-1] != "MedianExp":
            print('The first column must be MedianExp')
            sys.exit(1) 
    else:
        sys.exit(1) 

    # Cocatenate mRNA data
#    mRNA_data = pd.read_pickle(mRNA_data_loc+"/"+mRNA_annotation_data[0]+".pkl")
    mRNA_data = pd.read_pickle(mRNA_data_loc+"/"+mRNA_annotation_data[0]+".pkl.gz")
    mRNA_data = pd.merge(deg_data[['Name']], mRNA_data,
        how='inner', on='Name')
    mRNA_feature_name = pd.read_csv(mRNA_data_loc+"/"+mRNA_annotation_data[0]+"_feature_name.txt.gz",
                                      header=None,compression='gzip')
    mRNA_feature_name = mRNA_feature_name.values

    for mRNA_Annot in mRNA_annotation_data[1:]:
        range_data_temp = pd.read_pickle(mRNA_data_loc+"/"+mRNA_Annot+".pkl")
        
        # Filter genes
        range_data_temp = pd.merge(deg_data[['Name']], range_data_temp,
            how='inner', on='Name')
        # Read feature names
        feature_name_temp = pd.read_csv(mRNA_data_loc+"/"+mRNA_Annot+"_feature_name.txt.gz",
                                          header=None,compression='gzip')
        feature_name_temp = feature_name_temp.values

        # Check gene order is identical
        if sum(mRNA_data.Name != range_data_temp.Name) >0:
            raise Exception("gene name does not math!")

#        # Remove exon
        indx=feature_name_temp[:,0]!="promotor_annot_mask.rds"
        for i in range(mRNA_data.shape[0]):
            mRNA_data.values[i,1]=np.vstack((mRNA_data.values[i,1],range_data_temp.values[i,1][indx,]))

        # Update feature name
        mRNA_feature_name=np.concatenate((mRNA_feature_name,feature_name_temp[indx]))

    # Convert to dense matrix
    for i in range(mRNA_data.shape[0]):
        mRNA_data.values[i,1]=np.array(mRNA_data.values[i,1].todense())

    # Filter features
    has_annot=[np.sum(x,axis=1) for x in mRNA_data.values[:,1]]
    has_annot=np.vstack(has_annot)
    indx = np.sum(has_annot>0,axis=0) < 0
    for i in range(mRNA_data.shape[0]):
        mRNA_data.values[i,1]=mRNA_data.values[i,1][~indx,:]
    mRNA_feature_name=mRNA_feature_name[~indx,0]
#    print('mRNA_data.values[0,1]_shape', mRNA_data.values[0,1].shape)
#    print('mRNA_feature_name_len', len(mRNA_feature_name))

    # Cocatenate promoter data
#    promoter_data = pd.read_pickle(promoter_data_loc+"/"+promoter_annotation_data[0]+".pkl")
    promoter_data = pd.read_pickle(promoter_data_loc+"/"+promoter_annotation_data[0]+".pkl.gz")
    promoter_data = pd.merge(deg_data[['Name']], promoter_data,
        how='inner', on='Name')
    promoter_feature_name = pd.read_csv(promoter_data_loc+"/"+promoter_annotation_data[0]+"_feature_name.txt.gz",
                                          header=None,compression='gzip')
    promoter_feature_name = promoter_feature_name.values

    for mRNA_Annot in promoter_annotation_data[1:]:
        range_data_temp = pd.read_pickle(promoter_data_loc+"/"+mRNA_Annot+".pkl")
        
        # Filter genes
        range_data_temp = pd.merge(deg_data[['Name']], range_data_temp,
            how='inner', on='Name')
        # Read feature names
        feature_name_temp = pd.read_table(promoter_data_loc+"/"+mRNA_Annot+"_feature_name.txt.gz",
                                          header=None,compression='gzip')
        feature_name_temp = feature_name_temp.values

        # Check gene order is identical
        if sum(promoter_data.Name != range_data_temp.Name) >0:
            raise Exception("gene name does not math!")

        # Remove mask
        indx=feature_name_temp[:,0]!="promotor_annot_mask.rds"
        for i in range(promoter_data.shape[0]):
            promoter_data.values[i,1]=np.vstack((promoter_data.values[i,1],range_data_temp.values[i,1][indx,]))

        # update feature name
        promoter_feature_name=np.concatenate((promoter_feature_name,feature_name_temp[indx]))
    
    # Convert to dense matrix
    for i in range(promoter_data.shape[0]):
        promoter_data.values[i,1]=np.array(promoter_data.values[i,1].todense())

    # Filter features
    has_annot=[np.sum(x,axis=1) for x in promoter_data.values[:,1]]
    has_annot=np.vstack(has_annot)
##    indx = np.sum(has_annot>0,axis=0) < 30
    indx = np.sum(has_annot>0,axis=0) < 0
    for i in range(promoter_data.shape[0]):
        promoter_data.values[i,1]=promoter_data.values[i,1][~indx,:]
    promoter_feature_name=promoter_feature_name[~indx,0]
#    print('promoter_data.values[0,1]_shape', promoter_data.values[0,1].shape)    
#    print('promoter_feature_name_len', len(promoter_feature_name))

    # check all matched
    if (sum(deg_data.values[:,0] != mRNA_data.values[:,0]))>0:
        raise Exception("gene name does not math!")
    if (sum(deg_data.values[:,0] != promoter_data.values[:,0]))>0:
        raise Exception("gene name does not math!")

    return deg_data,mRNA_data,promoter_data, mRNA_feature_name, promoter_feature_name


def prep_ml_data_split3(deg_data_file,
                       mRNA_data_loc,
                       mRNA_annotation_data,
                       promoter_data_loc,
                       promoter_annotation_data,
                       test_genes,
                       outloc,
                       shuffle="None"):

    # deg_data_file: a full path to transcriptome data
    # mRNA_data_loc: a full path to a directory where mRNA (promoter in this study) annotation data is stored
    # mRNA_annotation_data: a file name of mRNA (promoter in this study) annotation data
    # promoter_data_loc:  a full path to a directory where promoter (enhancer in this study) annotation data is stored
    # promoter_annotation_data: a file name of promoter (enhancer in this study) annotation data
    # test_genes: a full path to a file contating gene ids for testing
    # outloc: a full path where scaling factors are saved
    # shuffle: which features are shuffled or not

    import pandas as pd
    import numpy as np

    deg_data, mRNA_data, promoter_data, mRNA_feature_name, promoter_feature_name=load_ml_data(deg_data_file,
                                                                                                         mRNA_data_loc,
                                                                                                         mRNA_annotation_data,
                                                                                                         promoter_data_loc,
                                                                                                         promoter_annotation_data)

    # Split data into training, validating, and testing subsets
    test = pd.read_csv(test_genes,sep="\t",header=None,compression='gzip').values[:,0]

    # Split mRNA feature
    X_mRNA_test=mRNA_data.query("Name in @test").copy()

    # Split promoter feature
    X_promoter_test=promoter_data.query("Name in @test").copy()

    # Split target data
    Y_test=deg_data.query("Name in @test").copy()

#    print('X_promoter_test2.shape', X_promoter_test.shape[0])
#    print('Y_test2.shape', Y_test.shape[0])
    # Scale fold changes
    std=Y_test.values[:,1:(Y_test.shape[1]-1)].std()
    if std==0:
        std=1
    for i in range(1,Y_test.shape[1]-1):
        Y_test.iloc[:,i]=Y_test.values[:,i]/std

    # Scale log2-TPM (assume it is located at the last column)
    std=Y_test.iloc[:,(Y_test.shape[1]-1)].std()
    Y_test.iloc[:,(Y_test.shape[1]-1)]=Y_test.values[:,(Y_test.shape[1]-1)]/std

    # Normalize mRNA features
    std = np.hstack(X_mRNA_test.values[:,1]).max(axis=1)

    # Special treatment for STD=0
    std[std==0]=1

    for i in range(len(X_mRNA_test.values[:,1])):
        X_mRNA_test.values[i,1]=((X_mRNA_test.values[i,1].transpose())/std).transpose()

    # Normalizing promoter features
    std = np.hstack(X_promoter_test.values[:,1]).max(axis=1)

    # Special treatment for STD=0
    std[std==0]=1

    # Scaling
    for i in range(len(X_promoter_test.values[:,1])):
        X_promoter_test.values[i,1]=((X_promoter_test.values[i,1].transpose())/std).transpose()

    # Shuffle features
    if shuffle=="all":
        np.random.seed(1234)
        shuffle_indices = np.random.permutation(np.arange(X_mRNA_test.shape[0]))
        X_mRNA_test.values[:,1] =X_mRNA_test.values[shuffle_indices,1]

        shuffle_indices = np.random.permutation(np.arange(X_promoter_test.shape[0]))
        X_promoter_test.values[:,1] =X_promoter_test.values[shuffle_indices,1]
    elif shuffle=="DNA":
        np.random.seed(1234)
        shuffle_indices = np.random.permutation(np.arange(X_promoter_test.shape[0]))
        X_promoter_test.values[:,1] =X_promoter_test.values[shuffle_indices,1]
    elif shuffle=="RNA":
        np.random.seed(1234)
        shuffle_indices = np.random.permutation(np.arange(X_mRNA_test.shape[0]))
        X_mRNA_test.values[:,1] =X_mRNA_test.values[shuffle_indices,1]

    return Y_test, X_mRNA_test, X_promoter_test

