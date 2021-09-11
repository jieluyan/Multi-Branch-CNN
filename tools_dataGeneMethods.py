from tools import *
def Normal(ne_path, test_po_path, novel_po_path, test_ne_path, novel_ne_path, ratio=1):
    createFolder("./temp/")
    change2ndfastaDistributionTheSameWith1st(test_po_path, ne_path, test_ne_path, ratio)
    executeCdhit2d(test_ne_path, ne_path, "./temp/ne_whole_delTest.fasta", c=1, n=2)
    change2ndfastaDistributionTheSameWith1st(novel_po_path, "./temp/ne_whole_delTest.fasta", novel_ne_path, ratio)
    return

def NormalTrain(train_po_path, test_ne_path, novel_ne_path, train_ne_path, whole_ne_path, ratio=1):
    createFolder("./temp/")
    cmd = 'cat %s %s > ./temp/test_add_novel.fasta' % (test_ne_path, novel_ne_path)
    msg = os.popen(cmd).read()
    executeCdhit2d('./temp/test_add_novel.fasta', whole_ne_path, "./temp/ne_whole_delTestNovel.fasta", c=1, n=2)
    change2ndfastaDistributionTheSameWith1st(train_po_path, "./temp/ne_whole_delTestNovel.fasta", train_ne_path, ratio)
    return

def CDHit(ne_path, test_po_path, novel_po_path, test_ne_path, novel_ne_path, ratio=1):
    msg, out_folder, novel_path, test_path, train_path = findPathOfCdhitAndSplitDataset(ne_path)
    if not(os.path.isfile(train_path)):
        msg, out_folder, novel_path, test_path, train_path = executeCdhitAndSplitDataset(ne_path)
    change2ndfastaDistributionTheSameWith1st(test_po_path, test_path, test_ne_path, ratio)
    change2ndfastaDistributionTheSameWith1st(novel_po_path, novel_path, novel_ne_path, ratio)
    return train_path

def CDHitTrain(train_po_path, train_ne_path, ne_path, ratio=1):
    msg, out_folder, novel_path, test_path, train_path = findPathOfCdhitAndSplitDataset(ne_path)
    if not(os.path.isfile(train_path)):
        msg, out_folder, novel_path, test_path, train_path = executeCdhitAndSplitDataset(ne_path)
    change2ndfastaDistributionTheSameWith1st(train_po_path, train_path, train_ne_path, ratio)
    return

def Shuffle(test_po_path, novel_po_path, test_ne_path, novel_ne_path):
    in_records, out_records = shuffleGeneNegaFromPosFasta(test_po_path, test_ne_path)
    in_records, out_records = shuffleGeneNegaFromPosFasta(novel_po_path, novel_ne_path)
    return

def ShuffleTrain(train_po_path, train_ne_path):
    in_records, out_records = shuffleGeneNegaFromPosFasta(train_po_path, train_ne_path)
    return

def Random(test_po_path, novel_po_path, test_ne_path, novel_ne_path):
    in_records, out_records = randomGeneNegaFromPosFasta(test_po_path, test_ne_path)
    in_records, out_records = randomGeneNegaFromPosFasta(novel_po_path, novel_ne_path)
    return

def RandomTrain(train_po_path, train_ne_path):
    in_records, out_records = randomGeneNegaFromPosFasta(train_po_path, train_ne_path)
    return