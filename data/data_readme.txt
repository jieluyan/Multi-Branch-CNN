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
4.  ./data/finalMdlSet is the sets to develop the final dataset each of them conatins:
    a. all the corresponding train sequences
    a. all the corresponding test sequences
    a. all the corresponding novel-test sequences