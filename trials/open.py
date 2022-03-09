with open("t-RNA/cov.fasta") as fh:
    sequence = fh.readlines()
coal = ""
for line in sequence[1:]:
	line = line[:-1]
	coal += line
print(coal)
# print(sequence[0])