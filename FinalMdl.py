from tools_baseClasses import *
from tools import *
best18fts_cdhit = ["type3Braac9", "type7raac19", "type8raac16", "type7raac18", "type1raac18", "type10raac18",
                   "type8raac14", "CKSAAP", "DDE", "type7raac17", "type11raac10", "type12raac16", "type12raac17",
                   "type8raac13", "KSCTriad", "type7raac16", "type1raac19", "type10raac19"]
default_hyperParaDict={"layer_num":2, "dropOutRate":0.0, "filter_num":64, "whole_tune_epoch":20, "epoch":1000}
best2fts = ['CKSAAP', 'type8raac18']
#ori_def_lrParas = {"lr":1e-4 ,"decay_rate":0.9 ,"decay_steps":100, "patience":20, "monitor":"loss"}
def_lrParas = {"lr":1e-4 ,"decay_rate":0.9 ,"decay_steps":50, "patience":20, "monitor":"loss"}
class FinalMdlSet():
    def __init__(self, activity="sodium",  dataGeneMethod="Normal", rep=0, init=True):
        self.dataGeneMethod, self.rep = dataGeneMethod, rep
        self.activity = activity
        self.root_dir = os.path.dirname(os.path.realpath(__file__))
        self.data_dir = os.path.join(self.root_dir, "data")
        self.final_set = os.path.join(self.data_dir, "finalMdlSet")
        createFolder(self.final_set)
        self.po_path = os.path.join(self.final_set, "%s_finalMdlSet_po.fasta" % activity)
        self.ne_path = os.path.join(self.final_set, "%s_finalMdlSet_ne.fasta" % activity)
        if init:
            self.GeneAllData()

    def GeneAllData(self):
        s = dataset4GeneMethods(self.activity, self.dataGeneMethod, self.rep)
        if not (os.path.isfile(self.po_path)):
            cmd = 'cat %s %s %s %s > %s' % (s.train_po_path, s.test_po_path, s.novel_po_path, s.short_po_path, self.po_path)
            msg = os.popen(cmd).read()
        if not (os.path.isfile(self.ne_path)):
            cmd = 'cat %s %s %s %s > %s' % (s.train_ne_path, s.test_ne_path, s.novel_ne_path, s.short_ne_path, self.ne_path)
            msg = os.popen(cmd).read()
        return

def_lrParas = {"lr": 1e-4 ,"decay_rate": 0.9 ,"decay_steps": 100, "patience": 20, "monitor": "loss"}
cal_lrParas = {"lr": 1e-4 ,"decay_rate": 0.86 ,"decay_steps": 50, "patience": 20, "monitor": "loss"}
lrParasDict = {"sodium": def_lrParas, "potassium": def_lrParas, "calcium": cal_lrParas}
hyperParaDict = {"layer_num": 2, "dropOutRate": 0.0, "filter_num": 64, "whole_tune_epoch": 20, "epoch": 1000}
class FinalCNNMdl():
    def __init__(self, best_ft_list=best18fts_cdhit, activity="calcium", dataGeneMethod="CDHit", hyperParaDict=hyperParaDict,
                 lrParasDict=lrParasDict, test_path="", rep=0, init=False):
        self.hyperParas = hyperParaDict
        self.lrParas = lrParasDict[activity]
        self.dataGeneMethods, self.activity, self.test_path = activity, dataGeneMethod, test_path
        f = FinalMdlSet(activity, dataGeneMethod, rep, init=True)
        self.ft_list = best_ft_list
        self.po_path, self.ne_path = f.po_path, f.ne_path
        self.testFlag = os.path.isfile(test_path)
        self.GeneOutPath()
        if init:
            self.developMdl()
        if self.testFlag:
            self.initTestRsPath()

    def GeneOutPath(self):
        self.root_dir = os.path.dirname(os.path.realpath(__file__))
        self.mdl_dir = os.path.join(self.root_dir, "CNNMdl")
        createFolder(self.mdl_dir)
        out_name = "CNNMdl-%s-%s" % (self.dataGeneMethods, self.activity)
        mdl_path = os.path.join(self.mdl_dir, "%s.mdl" % out_name)
        mdl_path = getUnexistedName(mdl_path, file_type=".mdl")
        pkl_path = os.path.join(self.mdl_dir, "%s.pkl" % out_name)
        pkl_path = getUnexistedName(pkl_path, file_type=".pkl")
        self.pkl_path, self.mdl_path, self.out_name = pkl_path, mdl_path, out_name
        return

    def GetOutPath(self):
        self.root_dir = os.path.dirname(os.path.realpath(__file__))
        self.mdl_dir = os.path.join(self.root_dir, "CNNMdl")
        createFolder(self.mdl_dir)
        out_name = "CNNMdl-%s-%s" % (self.dataGeneMethods, self.activity)
        mdl_path = os.path.join(self.mdl_dir, "%s.mdl" % out_name)
        mdl_path = getExistedLargestName(mdl_path, file_type=".mdl")
        pkl_path = os.path.join(self.mdl_dir, "%s.pkl" % out_name)
        pkl_path = getExistedLargestName(pkl_path, file_type=".pkl")
        self.pkl_path, self.mdl_path, self.out_name = pkl_path, mdl_path, out_name
        return

    def initTestRsPath(self):
        self.rs_dir = os.path.join(self.root_dir, "TestRs")
        createFolder(self.rs_dir)
        test_name = os.path.basename(self.test_path)
        out_name = os.path.basename(self.mdl_path)
        testRs_path = os.path.join(self.rs_dir, "testRs_%s_of_%s.csv" % (test_name, out_name))
        self.testRs_path = getUnexistedName(testRs_path, file_type=".csv")
        return

    def developMdl(self):
        self.sets = CNNimportFtsDataSetsPoNe(self.po_path, self.ne_path, self.ft_list)
        self.mdl = DevelopfitAndSaveCNNMdl(self.sets, self.mdl_path, self.pkl_path, self.hyperParas, self.lrParas)
        if self.testFlag:
            self.initTestRsPath()
            self.Y, self.pred_Y = predictSequenceFromSaveKerasMdl(self.test_path, self.mdl_path, self.testRs_path, self.ft_list)
        return

if __name__ == "__main__":
    activities = ["sodium", "potassium", "calcium"]
    # print("Starting develop and save model: ###########")
    # cnn = FinalCNNMdl(best18fts_cdhit, activity="sodium", init=True)
    mdl_path = {"sodium": "~/path-to-ion-project/ion-parallel-cnn/CNNMdl/CNNMdl-calcium-CDHit.pkl",
         "potassium": "~/path-to-ion-project/ion-parallel-cnn/CNNMdl/CNNMdl-potassium-CDHit.pkl",
         "calcium": "~/path-to-ion-project/ion-parallel-cnn/CNNMdl/CNNMdl-calcium-CDHit.pkl"}
    print("Starting predict fasta sequences: ###########")
    # test path is a fasta file
    test_path = "your_test_peptide_sequences.fasta"
    pred_path = "test_result_path.csv"
    pred = predictSequenceFromSaveKerasMdl(test_path, mdl_path["sodium"], pred_path, ft_list=best18fts_cdhit)
