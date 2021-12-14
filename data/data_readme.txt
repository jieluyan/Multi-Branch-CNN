1. calcium.fasta, potassium.fasta, sodium.fasta are fasta files that contain sequences >= 11 amino acid (a.a.)
2. ./data/len5ToLen10Seqs folders are fasta files that contain sequences in the length range [5 a.a. ~ 10 a.a.]
    the reason of different length in different files is that:
        we use CD-Hit, a bioinformatic tool, to handle original data, but CD-Hit delete all the sequences whose length < 11 a.a.
3. CDHit, Normal, Random, Shuffle folders are negative data generation methods
    each folder contain:
        a. 1 test negative dataset
        b. 1 novel-test negative dataset
        c. 10 train negative dataset
        (I always use train set with 0 order like "calcium_train_ne0.fasta",
        you can select one you like or you can also use multiples to do further investigating
        (develop different models by different train sets to calculate std or...))
4. For the novel-test, test and train sets used in comparison with TML13, TML13-Stack, and Single-Branch-CNN (in the Result section of our paper):
    a. For the positive sets:
        novel-test dataset positive part includes all the sequences of two files which are (use calcium channel as example):
            1) ./data/len5ToLen10Seqs/calcium_len5ToLen10Seqs.fasta (short sequences whose length range is [5,10], cause CD-HIT only processes the sequences whose length >= 11)
            2) ./data/split_po/calcium_novel.fasta (contain the sequnces from all the clusters that only contain 1 sequence by CD-HIT -n 2 -c 0.4)
        test dataset positive part includes all the sequences of three files which are (use calcium channel as example):
            1) ./data/len5ToLen10Seqs/calcium_len5ToLen10Seqs.fasta
            2) ./data/split_po/calcium_novel.fasta
            3) ./data/split_po/calcium_test.fasta
        train dataset positive part includes all the sequences of 1 file which are (use calcium channel as example):
            1) ./data/split_po/calcium_train.fasta
    b. for CDHit negative data generation method (use calcium channel as example):
        novel-test dataset negative part includes all the sequences of two files which are (use calcium channel as example):
            1) ./data/len5ToLen10Seqs/CDHit/calcium_len5ToLen10Seqs_ne0.fasta 
            2) ./data/CDHit/calcium_novel_ne.fasta
        test dataset negative part includes all the sequences of three files which are (use calcium channel as example):
            1) ./data/len5ToLen10Seqs/CDHit/calcium_len5ToLen10Seqs_ne0.fasta
            2) ./data/CDHit/calcium_novel_ne.fasta
            3) ./data/CDHit/calcium_test_ne.fasta
        train dataset positive part includes all the sequences of 1 file which are (use calcium channel as example):
            1) ./data/CDHit/calcium_train_ne0.fasta
4.  ./data/finalMdlSet is the sets to develop the final dataset each of them conatins:
    a. all the corresponding train sequences
    a. all the corresponding test sequences
    a. all the corresponding novel-test sequences
