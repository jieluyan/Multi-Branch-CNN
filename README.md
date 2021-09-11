# Ion-Prarllel-CNN

1. download the package
2. create a new conda environment and install the requirement packages follow:
        requirements.txt
3. activate the conda environment
4. if you want to test peptide sequences, do as the following codes:
    from tools_baseClasses import *
    from tools import *
    from FinalMdl import *
    Note: activity can be "sodium", "potassium", or "calcium"
    mdl_path = {"sodium": "~/path-to-ion-project/ion-parallel-cnn/CNNMdl/CNNMdl-calcium-CDHit.pkl",
         "potassium": "~/path-to-ion-project/ion-parallel-cnn/CNNMdl/CNNMdl-potassium-CDHit.pkl",
         "calcium": "~/path-to-ion-project/ion-parallel-cnn/CNNMdl/CNNMdl-calcium-CDHit.pkl"}
    print("Starting predict fasta sequences: ###########")
    Note: test path is a fasta file
    test_path = "your_test_peptide_sequences.fasta"
    pred_path = "test_result_path.csv" # path where you want to saved the prediction result
    pred = predictSequenceFromSaveKerasMdl(test_path, mdl_path["sodium"], pred_path, ft_list=best18fts_cdhit)

5. if you want to develop model by yourself, plese use the following codes to develop models:
    from tools_baseClasses import *
    from tools import *
    from FinalMdl import *
    Note: activity can be "sodium", "potassium", or "calcium"

    print("Starting develop and save model: ###########")
    cnn = FinalCNNMdl(best18fts_cdhit, activity="sodium", init=True)

6. for the data, please refer to the ./data/data_readme.txt
