from codes.draw_length_distribution_tools import *

def readFastaAndLeftSpecificLengthRange(in_fasta, out_fasta):
    records = []
    for record in SeqIO.parse(in_fasta, "fasta"):
        # get sequence length
        length = len(record.seq)
        if length in list(range(5, 11)):
        # get sequence name record.id = "P61542|ori"
            name = record.id.split("|")[0]
            record.id = name
            records.append(record)
    SeqIO.write(records, out_fasta, "fasta")
    return records

if __name__ == "__main__":
    in_fasta = "negative_remove_dupliate_seqs.fasta"
    out_fasta = "negative_len5ToLen10Seqs.fasta"
    records = readFastaAndLeftSpecificLengthRange(in_fasta, out_fasta)