import RNA
import random

def getShapeFromCSV(csv_file):
    retVec = []
    retVec.append(-999.0);  # data list is 1-based, so we add smth. at pos 0
    with open(csv_file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        line = line.split(',')
        if line[1] == 'nan\n':
            retVec.append(-999.0)
        else:
            retVec.append(float(line[1][:-1]))
    return retVec


def initialize_fold(sequence_file, output_file, size, shape_file = None, m=None, b=None):
    with open(sequence_file) as fh:
        string = fh.readlines()
    sequence = ""
    for line in string[1:]:
        if line[-1] == "\n":
            line = line[:-1]
        sequence += line
## preferable, more flexible method: local folding via fold_compound interface
    md = RNA.md()
    md.window_size = size
    md.max_bp_span = size
    fc_A = RNA.fold_compound(sequence,md,RNA.OPTION_WINDOW)
	# Check for shape data
		
    
    fc_B = fc_A
    
    shape_A = getShapeFromCSV(shape_file[0])
    shape_B = getShapeFromCSV(shape_file[1])
    fc_A.sc_add_SHAPE_deigan(shape_A,m,b)     # A with A SHAPE data   fc_A
    fc_A_B = fc_A
    fc_A_B.sc_add_SHAPE_deigan(shape_B,m,b)   # A with B SHAPE data   fc_A_B 
    fc_B.sc_add_SHAPE_deigan(shape_B,m,b)   # B with B SHAPE data   fc_B
    fc_B_A = fc_B
    fc_B_A.sc_add_SHAPE_deigan(shape_A,m,b) # B with A SHAPE data   fc_B_A

    file_A, file_AB, file_B, file_BA = 'A.fasta', 'AB.fasta', 'B.fasta', 'BA.fasta'
    with open(file_A,"w") as fh:
        mfe = fc_A.mfe_window(fh)
    with open(file_AB,"w") as fh:
        mfe = fc_A_B.mfe_window(fh)
    with open(file_B,"w") as fh:
        mfe = fc_B.mfe_window(fh)
    with open(file_BA,"w") as fh:
        mfe = fc_B_A.mfe_window(fh)

    return(file_A, file_AB, file_B, file_BA,sequence)
	# shape = getShapeFromCSV(shape_file)
	# # print(shape[:100])
	# fc.sc_add_SHAPE_deigan(shape,m,b)
        
    # with open(output_file,"w") as fh:
    #     mfe = fc.mfe_window(fh)
    # # file = open(output_file, mode = 'r', encoding = 'utf-8-sig')
    # # f = open("myfile2.txt", "x")
    # # f.write(file.read())
    # # f.close()
    # return(output_file,sequence)

def parse_and_merge(file_A, file_AB):
	file_A = open(file_A, mode = 'r', encoding = 'utf-8-sig')
	file_AB = open(file_AB, mode = 'r', encoding = 'utf-8-sig')
	lines = file_A.readlines()
	seq_A = []
	for line in lines:
		window = line.split()
        #window:  ['.(((((.........))))).', '(', '-0.10)', '50']
        # To remove the '(' member of the list
		if '(' in window:
			window.remove('(')
            # To remove the '-0.10)' parenthese from one side
			window[1] = window[1][:-1]
			seq_A.append(window)
		elif window != [] :
            # To remove the '(-10.3)' parentheses at both ends
			window[1] = window[1][1:-1]
			seq_A.append(window)
		
	lines = file_AB.readlines()
	seq_AB = []
	for line in lines:
		window = line.split()
        #window:  ['.(((((.........))))).', '(', '-0.10)', '50']
        # To remove the '(' member of the list
		if '(' in window:
			window.remove('(')
            # To remove the '-0.10)' parenthese from one side
			window[1] = window[1][:-1]
			seq_AB.append(window)
		elif window != [] :
            # To remove the '(-10.3)' parentheses at both ends
			window[1] = window[1][1:-1]
			seq_AB.append(window)
	count = 0 
	for i in range(len(seq_A)):
		if seq_A[i][2] == seq_AB[i][2]:
			count += 1
			print('done: ', count)
			seq_A[i][1] = float(seq_A[i][1]) + float(seq_AB[i][1])

	return seq_A

# res = initialize_fold('t-RNA/seq_isolat_PS.fasta', 'METHOD.lfold.chain',200,('full1.csv','full2.csv'), 1.8, -0.6)
res = parse_and_merge('AB.fasta', 'A.fasta')
res2 = parse_and_merge('B.fasta', 'BA.fasta')
print(res[:10],'\n')
print(res2[:10])