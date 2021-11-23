from tools_baseClasses import *
from tools import *
from FinalMdl import *
mdl_path = {"sodium": "./CNNMdl/CNNMdl-sodium-CDHit.mdl", "potassium": "./CNNMdl/CNNMdl-potassium-CDHit.mdl", "calcium": "./CNNMdl/CNNMdl-calcium-CDHit.mdl"}
print("Starting predict fasta sequences: ###########")
test_path = "./input.fasta"
#test_path = "./test.fasta"
createFolder("./test_rs/")
def test_fasta(input_fasta_path="./input.fasta", output_csv_path="./test_rs/input_test_result_3channels.csv"):
    for channel in ["sodium", "potassium", "calcium"]:
        pred_path = "./test_rs/input_test_result_%s.csv" % channel
        x, y, t_df = predictSequenceFromSaveKerasMdl(input_fasta_path, mdl_path[channel], pred_path, ft_list=best18fts_cdhit)
        #t_df = pred[3]
        if channel == "sodium":
            df = t_df
        else:
            df = pd.merge(df, t_df, how ='inner', on =['name', 'Sequence'])
    df.to_csv(output_csv_path)
    print("prediction score of %s for sodium, potassium, and calcium channel were written to %s." % (input_fasta_path, output_csv_path))
    print("Note: score is a float between 0 ~ 1, >0.5 is positive for the corresponding ion channel.")
    return

if __name__ == "__main__":
    test_fasta(input_fasta_path="./input.fasta", output_csv_path="./test_rs/input_test_result_3channels.csv")