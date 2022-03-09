import RNA
import random

#convert .map file to CSV
# with open("t-RNA/Mg_1M7_1_29903.map") as fh:
#     lines = fh.readlines()
# out = open('csv-2.csv','w')
# for line in lines:
# 	line = line.split()
# 	out.write(line[0])
# 	out.write(',')
# 	if line[1] == 0:
# 		out.write('nan')
# 	else:
# 		out.write(line[1])
# 	out.write("\n")
# 	print(line[0],line[1],line[2])


def getShapeDataFromFile(shape_file):
    retVec = []
    retVec.append(-999.0);  # data list is 1-based, so we add smth. at pos 0
    count=1
    with open(shape_file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        line = line.split()
        pos = int(line[0])
        value = float(line[1])

        if(pos==count):
            retVec.append(value)
        else:
            for i in range(pos-count):
                retVec.append(-999.0)
            retVec.append(value)
            count=pos
        count+=1
    return retVec


# res = getShapeDataFromFile('t-RNA/Mg_1M7_1_29903.map')
# print(res)
# A positive slope m penalizes high reactivities in paired regions, 
# while a negative intercept b results in a confirmatory `‘bonus’' 
# free energy for correctly predicted base pairs

#  m= 1.9
#  b = -0.7

def initialize_fold(sequence_file, output_file, size, shape_file = None, m=None, b=None, hc= None, sc = None):
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
    fc = RNA.fold_compound(sequence,md,RNA.OPTION_WINDOW)
	# Check for shape data
    if shape_file is not None:
        shape = getShapeDataFromFile(shape_file)
        fc.sc_add_SHAPE_deigan(shape,m,b)
    # Check for constraints
    if hc is not None:
        # constrain some bases as unpaired and fold again
        a,b = hc
        for i in range(a,b): 
            fc.hc_add_up(i) 
    if sc is not None:
        a,b = sc
        # add a base pair constraint
        fc.hc_add_bp(a,b) 
    with open(output_file,"w") as fh:
        mfe = fc.mfe_window(fh)
    file = open(output_file, mode = 'r', encoding = 'utf-8-sig')
    return(file,sequence)

# x = initialize_fold('t-RNA/shape-seq.fasta', 't-RNA/output-shape.fasta', 50, 't-RNA/shape-react.fasta', 1, -1)
# print(x)

def window_fold(input_file, output_file, size, shape_file = None, m=None, b=None, hc= None, sc = None):
    # Create and open file
    file, sequence = initialize_fold(input_file, output_file, size, shape_file, m, b)#the function returns (file,sequence)
    # Parse the t-RNA Sequence
    lines = file.readlines()
    seq = []
    for line in lines:
        window = line.split()
        #window:  ['.(((((.........))))).', '(', '-0.10)', '50']
        # To remove the '(' member of the list
        if '(' in window:
            window.remove('(')
            # To remove the '-0.10)' parenthese from one side
            window[1] = window[1][:-1]
            seq.append(window)
        elif window != [] :
            # To remove the '(-10.3)' parentheses at both ends
            window[1] = window[1][1:-1]
            seq.append(window)
    # Reformat data into useful lists
    intervals = []
    for window in seq:
        (left,right,energy) = (int(window[2]),'L'), (int(window[2])+len(window[0]),'R'), float(window[1])
        intervals.append((left,right,energy))   
    intervals.reverse()
    I = []
    for i in range(len(intervals)):
        # a c = location start/end
        # b d = Left or Right 
        # i = index of tuples from 'Interval'
        a, b = intervals[i][0]
        c, d = intervals[i][1]
        I.append((a,b,i))
        I.append((c,d,i))
    I.sort()
    #Algorithm
    minimum = 0
    temp_min = -1
    V = [ (j[2],[j]) for j in intervals]
    v = [ j[2] for j in intervals]
    for i in range(len(I)):
        # Left of interval j, V[j] = v[j] + minimum
        if I[i][1] == 'L':
            j = I[i][2]
            V_0 = minimum + v[j]
            if temp_min != -1:
                V_1 = V[temp_min][1]+[intervals[j]]
            else:
                V_1 = [intervals[j]]
            V[j] = (V_0,V_1)
        # Right of interval j, minimum = min(minimum, V[j])
        if I[i][1] == 'R':
            j = I[i][2]
            if V[j][0] <= minimum:
                minimum = V[j][0]
                temp_min = j  
    result = []
    for p in V[temp_min][1]:
        left, right, energy = p[0][0], p[1][0], p[2]
        result.append(([left,right], energy))
    return(result,sequence) 


def local_fold_list(input_file, output_file, size, shape_file = None, m=None, b=None, hc= None, sc = None):
    # Create and open file
    file = initialize_fold(input_file, output_file, size, shape_file, m, b)[0] #the function returns (file,sequence)
    # Parse the t-RNA Sequence
    lines = file.readlines()
    out = open('t-RNA/output_output.fasta', "w")
    n = out.write(file)
    out.close()
    d = {}
    for line in lines:
        window = line.split()
        #window:  ['.(((((.........))))).', '(', '-0.10)', '50']
        # dictionary with {50 : '.(((((.........))))).' , etc.}
        if window != []:
            d[int(window[-1])] =  window[0]
    return d


def concatinate_local_fold(input_file, output_file, size, shape_file = None, m=None, b=None, hc= None, sc = None):
    result, sequence = window_fold(input_file, output_file, size, shape_file, m, b) #the function returns (file,sequence)
    # results = [([1, 37], -16.8), ([46, 74], -3.3), ([76, 122], -6.4), ([123, 164], -9.1), ([167, 199], -11.7)]
    seq_len = len(sequence)
    # print(result)
    D = local_fold_list(input_file, output_file, size, shape_file, m, b)
    # D = {13: '.((((.((((((((((.(.......).))))))........))))..)))).',
    #       2: '.((((..(((((((.....))))))))))).',
    #       1: '(((((.(((((((((.....))))))).))))))).' , etc.}
    sequence += '\n'
    for i in range(len(result)):
        start = result[i][0][0]
        if i == 0:
            for j in range(start-1):
                sequence += '.'
        end = start+len(D[start])
        if i != len(result)-1:
            next = result[i+1][0][0]
        sequence += D[start]
        for j in range(next-end):
            sequence += '.'
        
        if i == len(result)-1:
            for j in range(end,seq_len+1):
                sequence += '.'
    output = ">"+ input_file + "\n" #set the header
    output += sequence + "\n"
    out = open(output_file, "w")
    n = out.write(output)
    out.close()
    return(output)

# res = concatinate_local_fold('cov.fasta',200,'output')
# res = concatinate_local_fold('t-RNA/shape-seq.fasta', 't-RNA/output-shape.fasta', 50, 't-RNA/shape-react.fasta', 1, -1)
#  m= 1.9
#  b = -0.7
# res = concatinate_local_fold('t-RNA/sequence.fasta', 'browser files/sequence.lfold.chain', 20)
res = concatinate_local_fold('t-RNA/seq_isolat_PS.fasta', 'browser files/COVID-300.lfold.chain',300,'t-RNA/Mg_1M7_1_29903.map', 1.9,-0.7)
# res = concatinate_local_fold('t-RNA/seq_isolat_PS.fasta', 'browser files/COVID-300-noshape.lfold.chain',300)
# res = concatinate_local_fold('t-RNA/sequence.fasta', 't-RNA/sequence-out.lfold.chain',20)
# res = window_fold('cov.fasta',200,'output')
# print(res)
 

# print('Heya')