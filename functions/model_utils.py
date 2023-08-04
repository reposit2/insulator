# This code was created based on codes in 'Tasaki, S., Gaiteri, C., Mostafavi, S. & Wang, Y. Deep learning decodes the principles of differential gene expression. Nature Machine Intelligence (2020)'

def test_prediction2(outloc,
                    best_model,
                    X_promoter_test,
                    Y_test):
    # Inputs
    # outloc: a full path to a directory of the best model
    # best_model: name of the best model
    # X_promoter_test: promoter annoatation data
    # Y_test: transcriptome data     
    
    # Outputs
    #.{outloc}/test_data/prediction.txt.gz: Predicted gene expression data
    #.{outloc}/test_data/actual.txt.gz: Actual gene expression data
    #.{outloc}/test_data/geneid.txt.gz: Genes in testing data.
    
    import sys
    import os
    sys.path.append('./functions/')
    import metrics
    from keras.models import load_model
    import data
    import numpy as np
    import os 
    
    if not os.path.exists(outloc+'test_data/'):
        os.makedirs(outloc+'test_data/')
    
    # Load best model
    model = load_model(outloc+best_model+'_model.h5',
                        custom_objects={'pcor': metrics.pcor})

    # Batching testing data
    batch_size=128
    test_steps, test_batches = data.batch_iter2(
                                               X_promoter_test.values[:,1],
                                               Y_test.values[:,1:],
                                               batch_size,
                                               shuffle=False)

    # Making prediction
    pred=[]
    actu=[]
    for i in range(test_steps):
        a=test_batches.__next__()
        b=model.predict(a[0])
        pred.append(b)
        actu.append(np.vstack(a[1]))

    pred=np.vstack(pred)
    actu=np.vstack(actu)
    
    # Save actual and predicted gene expression
    np.savetxt(outloc+'test_data/actual.txt',
               actu, delimiter='\t')
    np.savetxt(outloc+'test_data/prediction.txt',
               pred, delimiter='\t')
    X_promoter_test['Name'].to_csv(outloc+'test_data/geneid.txt', header=False, index=False, sep='\t')
    
    # gzip text files
    os.system("gzip "+outloc+'test_data/actual.txt')
    os.system("gzip "+outloc+'test_data/prediction.txt')    
    os.system("gzip "+outloc+'test_data/geneid.txt')    
    

def compute_DeepLIFT4(outloc,
             best_model,
             X_promoter_test,
             Y_test,
             X_promoter_test2,
             Y_test2):

    # Input
    # outloc: a full path to a directory of the best model
    # best_model: name of the best model
    # X_promoter_test: promoter and enhancer annoatation data
    # Y_test: transcriptome data    
    
    # Output
    #.{outloc}/{best_model}/DeepLIFT/DNA_{sample index}.txt.gz: DeepLIFT scores of DNA regulators for each sample. The sample index corresponds to the column index of gene expression data, which starts from 0. DeepLIFT scores of regulators were sparated with commas. Its order is identical to the one appeared in {outloc}/feature_norm_stats.txt   

    import sys
    import os
    sys.path.append('./functions/')
    import metrics
    from keras.models import load_model
    from keras.layers.core import Dense, Activation
    from keras.layers import Input, BatchNormalization, Concatenate, Conv1D
    from keras.models import Model, Sequential
    from keras.models import model_from_json
    import data
    import numpy as np
    import deeplift
    from deeplift.layers import NonlinearMxtsMode
    from deeplift.conversion import kerasapi_conversion as kc
    import tensorflow as tf

    def my_len(l):
        count = 0
        if isinstance(l, list):
            for v in l:
                count += my_len(v)
            return count
        else:
            return 1

    # Load model
    model = load_model(outloc+best_model+'_model.h5',
                        custom_objects={'pcor': metrics.pcor})

    # Get parameters for the last dens layer
    dens_parameter = model.layers[-1].get_weights()
    
    # Construct a single output model 
    # Adding new layers
    fc = Dense(1,activation='linear',name='out')(model.layers[-2].output)
    new_model = Model(inputs=model.input, outputs=fc)
    new_model = Sequential(layers=new_model.layers)
#    new_model2_file = outloc+'tmp_model.h5'

    # Paramter for background distribution
    med_promoter_len=int(np.median([x.shape[1] for x in X_promoter_test.values[:,1]]))
    med_promoter_len2=int(np.median([x.shape[1] for x in X_promoter_test2.values[:,1]]))
    gene_names_test = X_promoter_test.values[:,0]

    # DeepLIFT score for each sample
    for out_indx in range(0,1):
        # Set dens parameters
        new_model.layers[-1].set_weights([dens_parameter[0][:,out_indx:(out_indx+1)],dens_parameter[1][out_indx:(out_indx+1)]])
#        new_model.save(new_model2_file)
        if not os.path.exists(outloc+'DeepLIFT/'):
            os.makedirs(outloc+'DeepLIFT/')

        # Speficy output file names
        outfile_name_at2=outloc+'DeepLIFT/DNA_'+str(out_indx)+'.txt'

        # Batching testing data
#        batch_size=256*4
        batch_size=1000
        test_steps, test_batches = data.batch_iter_DeepLIFT2(
                                                            X_promoter_test.values[:,1],
                                                            Y_test.values[:,1:],
                                                            batch_size,
                                                            med_promoter_len,
                                                            shuffle=False)

        test_steps2, test_batches2 = data.batch_iter_DeepLIFT2(
                                                              X_promoter_test2.values[:,1],
                                                              Y_test2.values[:,1:],
                                                              batch_size,
                                                              med_promoter_len2,
                                                              shuffle=False)

        for i in range(test_steps):
            xs_test,ys_test=next(test_batches)
            xs_test2,ys_test2=next(test_batches2)
            xs_dna_test=xs_test

            # Reshape background
            xs_background=xs_test2
            xs_background=np.array(xs_background)
            xs_background=xs_background * 0.5
            ys_pred = new_model.predict(xs_test)

            # Compute DeepLIFT scores2
            h5_file=outloc+best_model+'_model.h5'
            revealcancel_model = kc.convert_model_from_saved_files(
                                 h5_file=outloc+best_model+'_model.h5',
                                 nonlinear_mxts_mode=NonlinearMxtsMode.DeepLIFT_GenomicsDefault)
            revealcancel_func = revealcancel_model.get_target_contribs_func(find_scores_layer_idx=0, target_layer_idx=-1)
            shap_values = np.array(revealcancel_func(
                              task_idx=0,
                              input_data_list=[xs_dna_test],
                              input_references_list=[xs_background],
#                              batch_size=500,
                              batch_size=250,
                              progress_update=None))
            shap_values = shap_values[:,30:230,:]

            with open(outfile_name_at2, 'a') as f_handle:
                    for j in range(shap_values.shape[0]):
                        seq_indx=xs_dna_test[j,30:230,0]>0
                        feature_vector=list(map(str, np.sum(shap_values[j,seq_indx,:],axis=0)))
                        out_txt=str(out_indx)+'\t'+gene_names_test[j+i*batch_size]+'\t'+','.join(feature_vector)+'\n'
                        f_handle.write(out_txt)
        # gzip text files
        os.system("gzip "+outfile_name_at2)

