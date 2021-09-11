import os
import numpy as np
import pandas as pd
from tools import *
from tools_dataGeneMethods import *
import json
from pandas import ExcelWriter

class dataset4GeneMethods():
    def __init__(self, activity="sodium", dataGeneMethod="Normal", rep=0, init=True):
        # self.activity = "sodium" # "sodium", "potassium", "calcium"
        # self.dataGeneMethod = "Normal" # "Normal", "CDHit", "Shuffle", "Random"
        self.activity, self.dataGeneMethod, self.rep = activity, dataGeneMethod, rep
        self.ori_data_folder = "data"
        self.getDataPath()
        if init:
            self.GeneAllData()

    def getDataPath(self):
        self.cur_folder = os.path.dirname(os.path.realpath(__file__))
        self.root_folder = self.cur_folder
        self.ori_po_name = "%s.fasta" % self.activity
        self.short_po_name = "%s_len5ToLen10Seqs.fasta" % self.activity
        self.data_folder = os.path.join(self.root_folder, "data")
        self.short_data_folder = os.path.join(self.data_folder, "len5ToLen10Seqs")
        self.po_oriPath = os.path.join(self.data_folder, self.ori_po_name )
        self.po_shortPath = os.path.join(self.short_data_folder, self.short_po_name)
        self.short_po_path = self.po_shortPath
        self.ne_oriPath = os.path.join(self.data_folder, "negative.fasta")
        self.short_ne_name = "negative_len5ToLen10Seqs.fasta"
        self.ne_shortWholePath = os.path.join(self.short_data_folder, self.short_ne_name)
        self.split_folder = os.path.join(self.data_folder, "split_po")
        self.ne_folder = os.path.join(self.data_folder, self.dataGeneMethod)
        createFolder(self.ne_folder)
        createFolder(self.split_folder)
        self.novel_po_path = os.path.join(self.split_folder, "%s_novel.fasta" % self.activity)
        self.novel_ne_path = os.path.join(self.ne_folder, "%s_novel_ne.fasta" % self.activity)
        self.test_po_path = os.path.join(self.split_folder, "%s_test.fasta" % self.activity)
        self.test_ne_path = os.path.join(self.ne_folder, "%s_test_ne.fasta" % self.activity)
        self.train_po_path = os.path.join(self.split_folder, "%s_train.fasta" % self.activity)
        self.train_ne_path = os.path.join(self.ne_folder, "%s_train_ne%d.fasta" % (self.activity, self.rep))
        self.short_ne_folder = os.path.join(self.short_data_folder, self.dataGeneMethod)
        createFolder(self.short_ne_folder)
        self.short_ne_path = os.path.join(self.short_ne_folder, "%s_len5ToLen10Seqs_ne%d.fasta" % (self.activity, self.rep))
        return

    def GenePoData(self):
        if not(os.path.isfile(self.novel_po_path)):
            msg, out_folder, novel_path, test_path, train_path = executeCdhitAndSplitDataset(self.po_oriPath)
            cmd = 'cp %s %s' % (novel_path, self.novel_po_path)
            msg = os.popen(cmd).read()
            cmd = 'cp %s %s' % (test_path, self.test_po_path)
            msg = os.popen(cmd).read()
            cmd = 'cp %s %s' % (train_path, self.train_po_path)
            msg = os.popen(cmd).read()
        return

    def GeneNeData(self):
        if not (os.path.isfile(self.novel_ne_path)):
            if self.dataGeneMethod == "Normal":
               Normal(self.ne_oriPath, self.test_po_path, self.novel_po_path, self.test_ne_path, self.novel_ne_path)
            if self.dataGeneMethod == "CDHit":
                CDHit(self.ne_oriPath, self.test_po_path, self.novel_po_path, self.test_ne_path, self.novel_ne_path)
            if self.dataGeneMethod == "Shuffle":
                Shuffle(self.test_po_path, self.novel_po_path, self.test_ne_path, self.novel_ne_path)
            if self.dataGeneMethod == "Random":
                Random(self.test_po_path, self.novel_po_path, self.test_ne_path, self.novel_ne_path)
        if not (os.path.isfile(self.train_ne_path)):
            if self.dataGeneMethod == "Normal":
                NormalTrain(self.train_po_path, self.test_ne_path, self.novel_ne_path, self.train_ne_path,
                            self.ne_oriPath)
            if self.dataGeneMethod == "CDHit":
                CDHitTrain(self.train_po_path, self.train_ne_path, self.ne_oriPath)
            if self.dataGeneMethod == "Shuffle":
                ShuffleTrain(self.train_po_path, self.train_ne_path)
            if self.dataGeneMethod == "Random":
                RandomTrain(self.train_po_path, self.train_ne_path)
        return

    def GeneShortNeData(self):
        if not(os.path.isfile(self.short_ne_path)):
            if self.dataGeneMethod == "Normal":
                change2ndfastaDistributionTheSameWith1st(self.short_po_path, self.ne_shortWholePath, self.short_ne_path)
            if self.dataGeneMethod == "CDHit":
                change2ndfastaDistributionTheSameWith1st(self.short_po_path, self.ne_shortWholePath, self.short_ne_path)
            if self.dataGeneMethod == "Shuffle":
                ShuffleTrain(self.short_po_path, self.short_ne_path)
            if self.dataGeneMethod == "Random":
                RandomTrain(self.short_po_path, self.short_ne_path)
        return

    def GeneAllData(self):
        self.GenePoData()
        self.GeneNeData()
        self.GeneShortNeData()
        return


class cat3Acts():
    def __init__(self, dataGeneMethod="Normal", rep=0, init=True):
        self.dataGeneMethod, self.rep = dataGeneMethod, rep
        all = dataset4GeneMethods("positive", self.dataGeneMethod, self.rep, init=False)
        self.train_po_path, self.test_po_path, self.novel_po_path = all.train_po_path, all.test_po_path, all.novel_po_path
        self.train_ne_path, self.test_ne_path, self.novel_ne_path = all.train_ne_path, all.test_ne_path, all.novel_ne_path
        self.short_po_path, self.short_ne_path = all.short_po_path, all.short_ne_path
        if init:
            self.GeneAllData()

    def GeneAllData(self):
        acts = ["sodium", "potassium", "calcium"]
        s = dataset4GeneMethods("sodium", self.dataGeneMethod, self.rep)
        p = dataset4GeneMethods("potassium", self.dataGeneMethod, self.rep)
        c = dataset4GeneMethods("calcium", self.dataGeneMethod, self.rep)
        if not (os.path.isfile(self.train_po_path)):
            cmd = 'cat %s %s %s > %s' % (s.train_po_path, p.train_po_path, c.train_po_path, self.train_po_path)
            msg = os.popen(cmd).read()
        if not (os.path.isfile(self.test_po_path)):
            cmd = 'cat %s %s %s > %s' % (s.test_po_path, p.test_po_path, c.test_po_path, self.test_po_path)
            msg = os.popen(cmd).read()
        if not (os.path.isfile(self.novel_po_path)):
            cmd = 'cat %s %s %s > %s' % (s.novel_po_path, p.novel_po_path, c.novel_po_path, self.novel_po_path)
            msg = os.popen(cmd).read()
        if not (os.path.isfile(self.train_ne_path)):
            cmd = 'cat %s %s %s > %s' % (s.train_ne_path, p.train_ne_path, c.train_ne_path, self.train_ne_path)
            msg = os.popen(cmd).read()
        if not (os.path.isfile(self.test_ne_path)):
            cmd = 'cat %s %s %s > %s' % (s.test_ne_path, p.test_ne_path, c.test_ne_path, self.test_ne_path)
            msg = os.popen(cmd).read()
        if not (os.path.isfile(self.novel_ne_path)):
            cmd = 'cat %s %s %s > %s' % (s.novel_ne_path, p.novel_ne_path, c.novel_ne_path, self.novel_ne_path)
            msg = os.popen(cmd).read()
        if not (os.path.isfile(self.short_ne_path)):
            cmd = 'cat %s %s %s > %s' % (s.short_ne_path, p.short_ne_path, c.short_ne_path, self.short_ne_path)
            msg = os.popen(cmd).read()
        if not (os.path.isfile(self.short_po_path)):
            cmd = 'cat %s %s %s > %s' % (s.short_po_path, p.short_po_path, c.short_po_path, self.short_po_path)
            msg = os.popen(cmd).read()
        return

class cat2PreTrained():
    def __init__(self, activity="pretrain-sodium", dataGeneMethod="Normal", rep=0, init=True):
        self.dataGeneMethod, self.rep = dataGeneMethod, rep
        self.activity = activity.split("-")[-1]
        print(self.activity)
        self.activities = ["sodium", "potassium", "calcium"]
        all = dataset4GeneMethods("pretrain-%s" % self.activity, self.dataGeneMethod, self.rep, init=False)
        self.train_po_path, self.test_po_path, self.novel_po_path = all.train_po_path, all.test_po_path, all.novel_po_path
        self.train_ne_path, self.test_ne_path, self.novel_ne_path = all.train_ne_path, all.test_ne_path, all.novel_ne_path
        self.short_po_path, self.short_ne_path = all.short_po_path, all.short_ne_path
        if init:
            self.GeneAllData()

    def GeneAllData(self):
        acts = self.activities
        print(self.activity)
        acts.remove(self.activity)
        d = {}
        for act in acts:
            d[act] = dataset4GeneMethods(act, self.dataGeneMethod, self.rep)
        if not (os.path.isfile(self.train_po_path)):
            cmd = 'cat %s %s > %s' % (d[acts[0]].train_po_path, d[acts[1]].train_po_path, self.train_po_path)
            msg = os.popen(cmd).read()
        if not (os.path.isfile(self.test_po_path)):
            cmd = 'cat %s %s > %s' % (d[acts[0]].test_po_path, d[acts[1]].test_po_path, self.test_po_path)
            msg = os.popen(cmd).read()
        if not (os.path.isfile(self.novel_po_path)):
            cmd = 'cat %s %s > %s' % (d[acts[0]].novel_po_path, d[acts[1]].novel_po_path, self.novel_po_path)
            msg = os.popen(cmd).read()
        if not (os.path.isfile(self.train_ne_path)):
            cmd = 'cat %s %s > %s' % (d[acts[0]].train_ne_path, d[acts[1]].train_ne_path, self.train_ne_path)
            msg = os.popen(cmd).read()
        if not (os.path.isfile(self.test_ne_path)):
            cmd = 'cat %s %s > %s' % (d[acts[0]].test_ne_path, d[acts[1]].test_ne_path, self.test_ne_path)
            msg = os.popen(cmd).read()
        if not (os.path.isfile(self.novel_ne_path)):
            cmd = 'cat %s %s > %s' % (d[acts[0]].novel_ne_path, d[acts[1]].novel_ne_path, self.novel_ne_path)
            msg = os.popen(cmd).read()
        if not (os.path.isfile(self.short_ne_path)):
            cmd = 'cat %s %s > %s' % (d[acts[0]].short_ne_path, d[acts[1]].short_ne_path, self.short_ne_path)
            msg = os.popen(cmd).read()
        if not (os.path.isfile(self.short_po_path)):
            cmd = 'cat %s %s > %s' % (d[acts[0]].short_po_path, d[acts[1]].short_po_path, self.short_po_path)
            msg = os.popen(cmd).read()
        return



class CNNfeaturesExtraction():
    def __init__(self, best_ft_list, activity="sodium", dataGeneMethod="Normal", rep=0, init=True):
        self.best_ft_list = best_ft_list
        self.activity, self.dataGeneMethod, self.rep = activity, dataGeneMethod, rep
        self.ift_types = ["AAC", "CKSAAP", "DPC", "DDE", "TPC", "GAAC",
                         "CKSAAGP", "GDPC", "GTPC", "NMBroto",
                         "Moran", "Geary", "CTDC", "CTDT", "CTDD", "CTriad", "KSCTriad",
                         "SOCNumber", "QSOrder", "PAAC", "APAAC"]
        self.kraacs = ["type1", "type2", "type3A", "type3B", "type4", "type5", "type6A",
                 "type6B", "type6C", "type7", "type8", "type9", "type10", "type11",
                 "type12", "type13", "type14", "type15", "type16"]
        self.subtype = {'g-gap': 0, 'lambda-correlation': 4}
        self.ktuple = 2
        self.gap_lambda = 1
        self.nlag=4
        self.lambdaValue=4
        self.useAllFts = False
        if self.activity == "positive":
            self.d = cat3Acts(self.dataGeneMethod, self.rep)
        elif "pretrian" in self.activity:
            self.d = cat2PreTrained(self.activity, self.dataGeneMethod, self.rep)
        else:
            self.d = dataset4GeneMethods(self.activity, self.dataGeneMethod, self.rep)
        if init:
            self.initDatasets()

    def initBestRaacNames(self):
        names = []
        if self.useAllFts:
            for ft_type in self.kraacs:
                AAGroup = eval("%s.AAGroup" % ft_type)
                for raactype in AAGroup:
                    if raactype == 20 or raactype < 5:
                        continue
                    name = "%sraac%d" % (ft_type, raactype)
                    names.append(name)
            self.best_raac_names = names
        else:
            self.best_raac_names = []
            self.best_trd_names = []
            for ft_name in self.best_ft_list:
                if "raac" in ft_name:
                    self.best_raac_names.append(ft_name)
                else:
                    self.best_trd_names.append(ft_name)
        return

    def importDataByPsekraac(self, po_path, ne_path, ft_name="type1"):
        po_df = genePsekraac(po_path, ft_name, self.raactype, self.subtype, self.ktuple, self.gap_lambda, class_val=1)
        ne_df = genePsekraac(ne_path, ft_name, self.raactype, self.subtype, self.ktuple, self.gap_lambda, class_val=0)
        df = pd.concat((po_df,ne_df))
        return df

    def importDataByiFtCodes(self, po_path, ne_path, ft_name="AAC"):
        po_df = GeneIfeature(po_path, ft_name, self.nlag, self.lambdaValue, class_val=1)
        ne_df = GeneIfeature(ne_path, ft_name, self.nlag, self.lambdaValue, class_val=0)
        df = pd.concat((po_df,ne_df))
        return df

    def addZeroCols(self, df):
        row, col = df.shape
        df = df.copy()
        p = self.ktuple
        add_col = int(np.power(np.ceil(np.power(col, 1/p)), p) - col)
        for i in range(add_col):
            col_name = "zero%d" % i
            df[col_name] = [0] * row
        return df

    def import1FtData(self, ft_whole_name="type1raac10"):
        novel_po_path, novel_ne_path, test_po_path, test_ne_path, train_po_path, train_ne_path = \
            self.d.novel_po_path, self.d.novel_ne_path, self.d.test_po_path, self.d.test_ne_path, self.d.train_po_path, self.d.train_ne_path
        short_po_path, short_ne_path = self.d.short_po_path, self.d.short_ne_path
        if "type" in ft_whole_name:
            self.raactype = int(ft_whole_name.split("raac")[-1])
            if "Ktuple" in ft_whole_name:
                ft_name = ft_whole_name.split("Ktuple")[0]
            else:
                ft_name = ft_whole_name.split("raac")[0]
            train_df = self.importDataByPsekraac(train_po_path, train_ne_path, ft_name)
            test_df = self.importDataByPsekraac(test_po_path, test_ne_path, ft_name)
            novel_df = self.importDataByPsekraac(novel_po_path, novel_ne_path, ft_name)
            novel_short_df = self.importDataByPsekraac(short_po_path, short_ne_path, ft_name)
        else:
            ft_name = ft_whole_name
            train_df = self.importDataByiFtCodes(train_po_path, train_ne_path, ft_name)
            test_df = self.importDataByiFtCodes(test_po_path, test_ne_path, ft_name)
            novel_df = self.importDataByiFtCodes(novel_po_path, novel_ne_path, ft_name)
            novel_short_df = self.importDataByiFtCodes(short_po_path, short_ne_path, ft_name)
        novel_df = pd.concat([novel_df, novel_short_df])
        return novel_df, test_df, train_df

    def ftset2ftImgAndLabel(self, df):
        label = df["default"].to_list()
        new_df = df.drop(["default"], 1)
        # add zero columns to a square or 3 power
        new_df = self.addZeroCols(new_df)
        new_df = new_df.to_numpy()
        row, col = new_df.shape
        label = np.array(label).reshape((row,1))
        label = label.astype(dtype=np.uint8)
        img_width = int(np.power(col,1/self.ktuple))
        if self.ktuple == 2:
            new_df = np.reshape(new_df, (row, img_width, img_width, 1))
        if self.ktuple == 3:
            new_df = np.reshape(new_df, (row, img_width, img_width, img_width))
        new_df = new_df.astype(dtype=np.float32)
        dict = {"img": new_df,
                "cls": label}
        return dict

    def geneAllRaacFeatureSets(self):
        novel_sets, test_sets, train_sets = {}, {}, {}
        ft_list = self.best_ft_list
        in_shapes = []
        for i in range(len(ft_list)):
            # for ft_whole_name in self.best_ft_list:
            # print("Processing file:\n\t": % novel_ne_path)
            ft_whole_name = ft_list[i]
            print("Processing feature: %s" % ft_whole_name)
            novel_df, test_df, train_df = self.import1FtData(ft_whole_name)
            novel_dict = self.ftset2ftImgAndLabel(novel_df)
            test_dict = self.ftset2ftImgAndLabel(test_df)
            train_dict = self.ftset2ftImgAndLabel(train_df)
            in_shapes.append([train_dict['img'].shape, test_dict['img'].shape, novel_dict['img'].shape])
            novel_sets[ft_whole_name] = novel_dict
            test_sets[ft_whole_name] = test_dict
            train_sets[ft_whole_name] = train_dict
        self.x_shapes = in_shapes
        self.train_po_num, self.train_ne_num = sum(train_df["default"]), len(train_df) - sum(train_df["default"])
        self.test_po_num, self.test_ne_num = sum(test_df["default"]), len(test_df) - sum(test_df["default"])
        self.novel_po_num, self.novel_ne_num = sum(novel_df["default"]), len(novel_df) - sum(novel_df["default"])
        print("Processing %s of %s:", (self.activity, self.dataGeneMethod))
        print("\tTrain positive No.: ", self.train_po_num)
        print("\tTrain negative No.: ", self.train_ne_num)
        print("\tTest positive No.: ", self.test_po_num)
        print("\tTest negative No.: ", self.test_ne_num)
        print("\tNovel positive No.: ", self.novel_po_num)
        print("\tNovel negative No.: ", self.novel_ne_num)
        return novel_sets, test_sets, train_sets

    def standardInputOutput(self, sets):
        X = []
        for ft_type in sets:
            x = sets[ft_type]['img']
            X.append(x)
            Y = sets[ft_type]['cls']
        return X, Y

    def initDatasets(self):
        self.novel_sets, self.test_sets, self.train_sets = \
                self.geneAllRaacFeatureSets()
        self.train_X, self.train_Y = self.standardInputOutput(self.train_sets)
        self.test_X, self.test_Y = self.standardInputOutput(self.test_sets)
        self.novel_X, self.novel_Y = self.standardInputOutput(self.novel_sets)
        return

class MLfeaturesExtraction():
    def __init__(self, best_ft_list=["AAC"], activity="sodium", dataGeneMethod="Normal", rep=0, init=True):
        self.best_ft_list = best_ft_list
        self.activity, self.dataGeneMethod, self.rep = activity, dataGeneMethod, rep
        self.ift_types = ["AAC", "CKSAAP", "DPC", "DDE", "TPC", "GAAC",
                         "CKSAAGP", "GDPC", "GTPC", "NMBroto",
                         "Moran", "Geary", "CTDC", "CTDT", "CTDD", "CTriad", "KSCTriad",
                         "SOCNumber", "QSOrder", "PAAC", "APAAC"]
        self.kraacs = ["type1", "type2", "type3A", "type3B", "type4", "type5", "type6A",
                 "type6B", "type6C", "type7", "type8", "type9", "type10", "type11",
                 "type12", "type13", "type14", "type15", "type16"]
        self.initAllFtNames()
        self.subtype = {'g-gap': 0, 'lambda-correlation': 4}
        self.ktuple = 2
        self.gap_lambda = 1
        self.nlag=4
        self.lambdaValue=4
        self.useAllFts = False
        if self.activity == "positive":
            self.d = cat3Acts(self.dataGeneMethod, self.rep)
        elif "pretrian" in self.activity:
            self.d = cat2PreTrained(self.activity, self.dataGeneMethod, self.rep)
        else:
            self.d = dataset4GeneMethods(self.activity, self.dataGeneMethod, self.rep)
        if init:
            self.initDatasets()

    def initAllFtNames(self):
        raac_list = []
        raac_name = []
        ft_names = []
        for ft_type in self.kraacs:
            AAGroup = eval("%s.AAGroup" % ft_type)
            for raactype in AAGroup:
                if raactype < 20:
                    raac_name.append("%sraac%d" % (ft_type, raactype))
                    if raactype > 5:
                        raac_list.append([ft_type, raactype])
        self.all_kraac_list = raac_list
        self.all_kraac_typeNames = raac_name
        ft_names.extend(raac_name)
        ft_names.extend(self.ift_types)
        self.all_ft_typeNames = ft_names
        return

    def importDataByPsekraac(self, po_path, ne_path, ft_name="type1"):
        po_df = genePsekraac(po_path, ft_name, self.raactype, self.subtype, self.ktuple, self.gap_lambda, class_val=1)
        ne_df = genePsekraac(ne_path, ft_name, self.raactype, self.subtype, self.ktuple, self.gap_lambda, class_val=0)
        df = pd.concat((po_df,ne_df))
        return df

    def importDataByiFtCodes(self, po_path, ne_path, ft_name="AAC"):
        po_df = GeneIfeature(po_path, ft_name, self.nlag, self.lambdaValue, class_val=1)
        ne_df = GeneIfeature(ne_path, ft_name, self.nlag, self.lambdaValue, class_val=0)
        df = pd.concat((po_df,ne_df))
        return df

    def import1FtData(self, ft_whole_name="type1raac10"):
        novel_po_path, novel_ne_path, test_po_path, test_ne_path, self.train_po_path, self.train_ne_path = \
            self.d.novel_po_path, self.d.novel_ne_path, self.d.test_po_path, self.d.test_ne_path, self.d.train_po_path, self.d.train_ne_path
        short_po_path, short_ne_path = self.d.short_po_path, self.d.short_ne_path
        if "type" in ft_whole_name:
            self.raactype = int(ft_whole_name.split("raac")[-1])
            if "Ktuple" in ft_whole_name:
                ft_name = ft_whole_name.split("Ktuple")[0]
            else:
                ft_name = ft_whole_name.split("raac")[0]
            train_df = self.importDataByPsekraac(self.train_po_path, self.train_ne_path, ft_name)
            test_df = self.importDataByPsekraac(test_po_path, test_ne_path, ft_name)
            novel_df = self.importDataByPsekraac(novel_po_path, novel_ne_path, ft_name)
            novel_short_df = self.importDataByPsekraac(short_po_path, short_ne_path, ft_name)
        else:
            ft_name = ft_whole_name
            train_df = self.importDataByiFtCodes(self.train_po_path, self.train_ne_path, ft_name)
            test_df = self.importDataByiFtCodes(test_po_path, test_ne_path, ft_name)
            novel_df = self.importDataByiFtCodes(novel_po_path, novel_ne_path, ft_name)
            novel_short_df = self.importDataByiFtCodes(short_po_path, short_ne_path, ft_name)
        novel_df = pd.concat([novel_df, novel_short_df])
        return novel_df, test_df, train_df

    def ftset2ftImgAndLabel(self, df):
        label = df["default"].to_list()
        new_df = df.drop(["default"], 1)
        dict = {"img": new_df.copy(),
                "cls": label}
        return dict

    def standardInputOutput(self, sets):
        c = -1
        for ft_type in sets:
            c += 1
            if c == 0:
                ft = sets[ft_type]['img']
            else:
                temp = sets[ft_type]['img']
                ft = pd.concat([ft, temp], axis=1)
        row, col = ft.shape
        col_names = ["ft%d" % i for i in range(col)]
        ft.columns = col_names
        ft['default'] = sets[ft_type]['cls']
        return ft

    def geneAllRaacFeatureSets(self):
        novel_sets, test_sets, train_sets = {}, {}, {}
        ft_list = self.best_ft_list
        in_shapes = []
        for i in range(len(ft_list)):
            # for ft_whole_name in self.best_ft_list:
            # print("Processing file:\n\t": % novel_ne_path)
            ft_whole_name = ft_list[i]
            print("Processing feature: %s" % ft_whole_name)
            novel_df, test_df, train_df = self.import1FtData(ft_whole_name)
            novel_dict = self.ftset2ftImgAndLabel(novel_df)
            test_dict = self.ftset2ftImgAndLabel(test_df)
            train_dict = self.ftset2ftImgAndLabel(train_df)
            novel_sets[ft_whole_name] = novel_dict
            test_sets[ft_whole_name] = test_dict
            train_sets[ft_whole_name] = train_dict
        self.train_po_num, self.train_ne_num = sum(train_df["default"]), len(train_df) - sum(train_df["default"])
        self.test_po_num, self.test_ne_num = sum(test_df["default"]), len(test_df) - sum(test_df["default"])
        self.novel_po_num, self.novel_ne_num = sum(novel_df["default"]), len(novel_df) - sum(novel_df["default"])
        print("Processing %s of %s:", (self.activity, self.dataGeneMethod))
        print("\tTrain positive No.: ", self.train_po_num)
        print("\tTrain negative No.: ", self.train_ne_num)
        print("\tTest positive No.: ", self.test_po_num)
        print("\tTest negative No.: ", self.test_ne_num)
        print("\tNovel positive No.: ", self.novel_po_num)
        print("\tNovel negative No.: ", self.novel_ne_num)
        return novel_sets, test_sets, train_sets

    def initDatasets(self):
        self.novel_sets, self.test_sets, self.train_sets = \
                self.geneAllRaacFeatureSets()
        self.train = self.standardInputOutput(self.train_sets)
        self.test = self.standardInputOutput(self.test_sets)
        self.novel = self.standardInputOutput(self.novel_sets)
        self.x_shape = [self.train.shape, self.test.shape, self.novel.shape]
        return

class Metrics():
    def __init__(self, ori, sco, threshold=0.5):
        self.ori, self.sco, self.threshold = ori, sco, threshold
        self.calMetrics()

    def F1(self):
        ori, sco, threshold = self.ori, self.sco, self.threshold
        recall = metrics.Recall(thresholds=threshold)
        recall.update_state(ori, sco)
        Recall = tf.cast(recall.result(), tf.float32)
        precision = metrics.Precision(thresholds=threshold)
        precision.update_state(ori, sco)
        Precision = tf.cast(precision.result(), tf.float32)
        one = tf.constant(1, dtype=tf.float32)
        two = tf.constant(2, dtype=tf.float32)
        f1 = tf.divide(two, tf.add(tf.divide(one, Recall), tf.divide(one, Precision)))
        return f1

    def Kappa(self):
        ori, sco, threshold = self.ori, self.sco, self.threshold
        pre = tf.math.round(sco)
        new_ori = tf.squeeze(ori)
        new_pre = tf.squeeze(pre)
        [tn, fp], [fn, tp] = tf.cast(tf.math.confusion_matrix(new_ori, new_pre), tf.float32)
        denominator = tf.add(tf.multiply(tf.add(tp, fp), tf.add(fp, tn)), tf.multiply(tf.add(tp, fn), tf.add(fn, tn)))
        two = tf.constant(2, dtype = tf.float32)
        numerator = tf.multiply(two, tf.add(tf.multiply(tp, tn), -tf.multiply(fn, fp)))
        kappa = tf.divide(numerator, denominator)
        return kappa

    def MCC(self):
        ori, sco, threshold = self.ori, self.sco, self.threshold
        pre = tf.math.round(sco)
        new_ori = tf.squeeze(ori)
        new_pre = tf.squeeze(pre)
        [tn, fp], [fn, tp] = tf.cast(tf.math.confusion_matrix(new_ori, new_pre), tf.float32)
        numerator = tf.add(tf.multiply(tp, tn), -tf.multiply(fp, fn))
        denominator = tf.multiply(tf.multiply(tf.add(tp, fp), tf.add(tp, fn)),
                                  tf.multiply(tf.add(tn, fn), tf.add(tn, fp)))
        mcc = tf.divide(numerator, tf.sqrt(denominator))
        return mcc

    def TNR(self):
        ori, sco, threshold = self.ori, self.sco, self.threshold
        pre = tf.math.round(sco)
        new_ori = tf.squeeze(ori)
        new_pre = tf.squeeze(pre)
        # print("ori: ", ori)
        # print("sco: ", sco)
        # print(tf.math.confusion_matrix(new_ori, new_pre))
        [tn, fp], [fn, tp] = tf.cast(tf.math.confusion_matrix(new_ori, new_pre), tf.float32)
        tnr = tf.divide(tn, tf.add(tn, fp))
        return tnr

    def calMetrics(self):
        ori, sco, threshold = self.ori, self.sco, self.threshold
        #sco = tf.keras.activations.sigmoid(tf.convert_to_tensor(sco))
        sco = tf.convert_to_tensor(sco)
        sco = sco[:,1]
        ori = tf.convert_to_tensor(ori)
        pre = tf.math.round(sco)
        new_ori = tf.squeeze(ori)
        new_pre = tf.squeeze(pre)
        self.ori, self.sco, self.pre = ori, sco, pre
        # print(tf.cast(tf.math.confusion_matrix(new_ori, new_pre), tf.float32))
        [tn, fp], [fn, tp] = tf.cast(tf.math.confusion_matrix(new_ori, new_pre), tf.float32)
        self.tn, self.fp, self.fn, self.tp = tn, fp, fn, tp
        all = tf.add_n([tn, fp, fn, tp])
        p = tf.add(tp,tn)
        Accuracy = tf.math.divide(p,all)
        # acc = metrics.BinaryAccuracy(threshold=threshold); acc.update_state(ori, sco); Accuracy = acc.result().numpy()
        auc = metrics.AUC()
        auc.update_state(ori,sco)
        AUC = auc.result()
        self.auc, self.acc = AUC, Accuracy
        recall = metrics.Recall(thresholds=threshold)
        recall.update_state(ori,sco)
        Recall = tf.cast(recall.result(), tf.float32)
        precision = metrics.Precision(thresholds=threshold)
        precision.update_state(ori, sco)
        Precision = tf.cast(precision.result(), tf.float32)
        self.recall, self.precision = Recall, Precision
        one = tf.constant(1, dtype = tf.float32)
        two = tf.constant(2, dtype = tf.float32)
        F1 = tf.divide(two, tf.add(tf.divide(one, Recall), tf.divide(one, Precision)))
        denominator = tf.add(tf.multiply(tf.add(tp, fp), tf.add(fp, tn)), tf.multiply(tf.add(tp, fn), tf.add(fn, tn)))
        numerator = tf.multiply(two, tf.add(tf.multiply(tp, tn), -tf.multiply(fn, fp)))
        Kappa = tf.divide(numerator, denominator)
        numerator = tf.add(tf.multiply(tp, tn), -tf.multiply(fp, fn))
        denominator = tf.multiply(tf.multiply(tf.add(tp, fp), tf.add(tp, fn)),
                                  tf.multiply(tf.add(tn, fn), tf.add(tn, fp)))
        MCC = tf.divide(numerator, tf.sqrt(denominator))
        # MAE, MSE, RMSE, R2, RMSLE, MAPE.
        TNR = tf.divide(tn, tf.add(tn, fp))
        names = ["Accuracy", "AUC", "Sp/TNR", "Recall/Sn/TPR", "Prec.", "F1", "Kappa", "MCC"]
        self.kappa, self.f1, self.mcc, self.tnr = Kappa, F1, MCC, TNR
        perf = tf.stack([Accuracy, AUC, TNR, Recall, Precision, F1, Kappa, MCC],axis=0)
        self.metrics = perf
        return

if __name__ == "__main__":
    # test class of dataset4GeneMethods()
    # dataGeneMethods = ["Normal", "CDHit", "Shuffle", "Random"]
    # acts = ["sodium", "potassium", "calcium"]
    # for gene_method in dataGeneMethods:
    #     for act in acts:
    #         for rep in range(10):
    #             d = dataset4GeneMethods(act, gene_method, rep)


    # test class of cat3Acts()
    dataGeneMethods = ["Normal", "CDHit", "Shuffle", "Random"]
    for gene_method in dataGeneMethods:
        for rep in range(10):
            c = cat3Acts(gene_method, rep)


    # test class of CNNfeaturesExtraction()
    # dataGeneMethods = ["Normal", "CDHit", "Shuffle", "Random"]
    # acts = ["sodium", "potassium", "calcium"]
    #
    # for gene_method in dataGeneMethods:
    #     for act in acts:
    #         for rep in range(10):
    #             f = CNNfeaturesExtraction(best10fts, act, gene_method, rep)

    # test class of MLfeatureExtraction()
    # dataGeneMethods = ["Normal", "CDHit", "Shuffle", "Random"]
    # acts = ["sodium", "potassium", "calcium"]
    # best10fts = ['type7raac18', 'type7raac19', 'CKSAAP', 'type8raac18',
    #              'type8raac15', 'type7raac16', 'type9raac20', 'type14raac18',
    #              'type14raac17', 'DPC']
    #
    # for gene_method in dataGeneMethods:
    #     for act in acts:
    #         for rep in range(10):
    #             f = MLfeatureExtraction(best10fts, act, gene_method, rep)

    # test class of cat2PreTrained()
    # dataGeneMethods = ["Normal", "CDHit", "Shuffle", "Random"]
    # acts = ["sodium", "potassium", "calcium"]
    # for gene_method in dataGeneMethods:
    #     for act in acts:
    #         tar_act = "pretrain-%s" % act
    #         for rep in range(10):
    #             d = cat2PreTrained(tar_act, gene_method, rep)