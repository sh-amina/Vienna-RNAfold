import RNA

# Parse the t-RNA Sequence
file = open('t-RNA/sequence (2).FASTA', mode = 'r', encoding = 'utf-8-sig')
lines = file.readlines()
seq = ''
for line in lines[1:]:
	seq += line.strip()

# compute minimum free energy (MFE) and corresponding structure
(ss, mfe) = RNA.fold(seq)
# print output
print("{}\n\n{} [ {:6.2f} ]".format(seq, ss, mfe))