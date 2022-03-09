with open("t-RNA/Mg_1M7_1_29903.map") as fh:
    lines = fh.readlines()
out = open('csv-2.csv','w')
for line in lines:
	line = line.split()
	out.write(line[0])
	out.write(',')
	if line[1] == 0:
		out.write('nan')
	else:
		out.write(line[1])
	out.write("\n")
	print(line[0],line[1],line[2])

