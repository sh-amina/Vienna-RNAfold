import RNA
import random
import pandas as pd

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
# res = getShapeFromCSV('full1.csv')

def get_SHAPE_average(file1,file2):
    retVec = []
    retVec.append(-999.0);  # data list is 1-based, so we add smth. at pos 0
    with open(file1, 'r') as f:
        lines1 = f.readlines()
    with open(file2, 'r') as f:
        lines2 = f.readlines()
    for i in range(len(lines1)):
        line1 = lines1[i].split(',')
        if line1[1] == 'nan\n':
            a = -999
        else:
            a = float(line1[1][:-1])

        line2 = lines2[i].split(',')
        if line2[1] == 'nan\n':
            b = -999
        else:
            b = float(line2[1][:-1])
        avr = (a+b)*0.5
        retVec.append(avr)
    # print(retVec)
    return retVec

# res = get_SHAPE_average('full1.csv','full2.csv')

# A positive slope m penalizes high reactivities in paired regions, 
# while a negative intercept b results in a confirmatory `‘bonus’' 
# free energy for correctly predicted base pairs
#  m = 1.9
#  b = -0.7
def initialize_fold(sequence_file, output_file, size, shape_file = None, shape_flag = None, m=None, b=None, hc= None, sc = None):
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
        if len(shape_file) == 2:
            if shape_flag:
                file1, file2 = shape_file
                shape = get_SHAPE_average(file1,file2)
                fc.sc_add_SHAPE_deigan(shape,m,b)
            else:
                fc2 = fc
                shape1 = getShapeFromCSV(shape_file[0])
                shape2 = getShapeFromCSV(shape_file[1])
                fc.sc_add_SHAPE_deigan(shape1,m,b)
                fc2.sc_add_SHAPE_deigan(shape2,m,b)
                file1, file2 = 'file1.fasta', 'file2.fasta'
                with open(file1,"w") as fh:
                    mfe = fc.mfe_window(fh)
                with open(file2,"w") as fh:
                    mfe = fc2.mfe_window(fh)
                return(file1,file2,sequence)
        else:
            shape = getShapeFromCSV(shape_file)
            # print(shape[:100])
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
    # file = open(output_file, mode = 'r', encoding = 'utf-8-sig')
    # f = open("myfile2.txt", "x")
    # f.write(file.read())
    # f.close()
    return(output_file,sequence)


def window_fold(input_file, output_file, size, shape_file = None, shape_flag = None, m=None, b=None, hc= None, sc = None):
    # Create and open file
    if  not shape_flag and len(shape_file) == 2:
        # file2, sequence = initialize_fold(input_file, output_file, size, shape_file[0], m, b)#the function returns (file,sequence)
        # file1 = initialize_fold(input_file, output_file, size, shape_file[1], m, b)[0]
        file1, file2, sequence = initialize_fold(input_file, output_file, size, shape_file, shape_flag, m, b)#the function returns (file1,file2,sequence)
        file2 = open(file2, "r")
        data2 = file2.read()
        file2.close()
        file1 = open(file1, "r")
        data1 = file1.read()
        file1.close()
        file = open('mix.fasta', "a")
        file.write(data2)
        file.write(data1)
        file.close()
        file = open('mix.fasta', mode = 'r', encoding = 'utf-8-sig')
    else:
        file, sequence = initialize_fold(input_file, output_file, size, shape_file, shape_flag, m, b)#the function returns (file,sequence)
        file = open(file, mode = 'r', encoding = 'utf-8-sig')
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
        (left,right,energy,string) = (int(window[2]),'L'), (int(window[2])+len(window[0]),'R'), float(window[1]), window[0]
        intervals.append((left,right,energy,string))   
    intervals.reverse()
    I = []
    for i in range(len(intervals)):
        # a,c = location start/end     b,d = Left or Right   i = index of tuples from 'Interval'
        a, b = intervals[i][0]
        c, d = intervals[i][1]
        I.append((a,b,i))
        I.append((c,d,i))
    I.sort()
    #Algorithm
    minimum = 0
    temp_min = -1
    V = [ (j[2],[j]) for j in intervals]  # [(energy, [left,right,energy,string]),....]
    v = [ j[2] for j in intervals]        # [energy,....]

    for i in range(len(I)):
        # Left of interval j, V[j] = v[j] + minimum
        # intervals (left endpoint, right endpoint, v(j))
        # I : (position, left/right, j pointer to intervals)
    
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
    # print('V: \n',V,'\n')
    result = []
    d = {}
    for p in V[temp_min][1]:
        # print('p :',p,'\n')
        left, right, energy = p[0][0], p[1][0], p[2]
        d[left] = p[3]
        result.append(([left,right], energy))
    # print('res: \n',result, '\n')
    # print('d: ', d, '\n')
    return(result,sequence,d) 

# res = window_fold('cov700.fasta', 'merge.lfold.chain', 200, 'cov1.csv' , 1.9, -0.7)
# res = window_fold('cov700.fasta', 'merge.lfold.chain', 200, ('cov1.csv','cov2.csv') , 1.9, -0.7)

def concatinate_local_fold(input_file, output_file, size, shape_file = None, shape_flag = None, m=None, b=None, hc= None, sc = None):
    result, sequence, D = window_fold(input_file, output_file, size, shape_file, shape_flag, m, b) #the function returns (file,sequence,d)
    # results = [([1, 37], -16.8), ([46, 74], -3.3), ([76, 122], -6.4), ([123, 164], -9.1), ([167, 199], -11.7)]
    seq_len = len(sequence)
    # D = local_fold_list(input_file, output_file, size, shape_file, m, b)
    # D = {13: '.((((.((((((((((.(.......).))))))........))))..)))).',
    #       2: '.((((..(((((((.....))))))))))).',
    #       1: '(((((.(((((((((.....))))))).))))))).' , etc.}
    sequence += '\n'
    next = 1
    energy = 0
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
        energy += result[i][1]
    output = ">"+ input_file + "\n" #set the header
    output += sequence + " ("
    output +=  str(energy) + ") \n"
    out = open(output_file, "w")
    n = out.write(output)
    out.close()
    return(output)

def similarity(input_1, input_2):
    with open(input_1, 'r') as fp:
        seq_1 = fp.readlines()[2]
    with open(input_2, 'r') as fp:
        seq_2 = fp.readlines()[2]
    yes = 0
    count = 0
    for i, c in enumerate(seq_1):
        if c == ' ':
            break
        if c == seq_2[i]:
            yes += 1
        count += 1
    percentage = f"{100*(yes/count):.2f}%"
    return percentage

def SHAPE_similarity(input_chain, input_csv, threshold_num):
    df = pd.read_csv (input_csv, header=None, names=["index", "value"])
    df = df.dropna()
    all = df.size
    df = df[df["value"] > threshold_num]
    threshold = df["index"].tolist()
    with open(input_chain, 'r') as fp:
        seq = fp.readlines()[2]
    count = len(threshold)
    yes = 0
    for index in threshold:
        if seq[index-1] == '.':
            yes += 1
    percentage = f"{100*(yes/count):.2f}%"
    return percentage

# res = similarity('COVID-200.lfold.chain', 'average-200.lfold.chain')
# res = SHAPE_similarity('full2/COVID-400.lfold.chain', 'full2.csv', 1)
# print(res)
for i in ['100','150','200','250','300','350','400']:
    res = similarity('full1/COVID-'+i+'.lfold.chain', 'full2/COVID-'+i+'.lfold.chain')
    print(i,res)
      

# res = concatinate_local_fold('cov700.fasta', 'merge.lfold.chain', 100, ('cov1.csv','cov2.csv') , 1.9, -0.7)
# import time
# start_time = time.time()
# res = concatinate_local_fold('t-RNA/seq_isolat_PS.fasta', 'full2/COVID-300.lfold.chain',300,'full2.csv',False, 1.8, -0.6)
# print("--- %s seconds ---" % (time.time() - start_time))

# res = concatinate_local_fold('cov-162.fasta', 'res-cov-81.lfold.chain',300,'full2.csv',False, 1.8, -0.6)


# res = concatinate_local_fold('fullcov.fasta', 'average-250.lfold.chain', 250, ('full1.csv','full2.csv') ,True, 1.9, -0.7)
# res = concatinate_local_fold('fullcov.fasta', 'double-320.lfold.chain', 320, ('full1.csv','full2.csv') , 1.9, -0.7)
# res = concatinate_local_fold('cov.fasta',200,'output')
# res = concatinate_local_fold('t-RNA/shape-seq.fasta', 't-RNA/output-shape.fasta', 50, 't-RNA/shape-react.fasta', 1, -1)
# res = concatinate_local_fold('t-RNA/sequence.fasta', 'browser files/sequence.lfold.chain', 20)
# res = concatinate_local_fold('t-RNA/seq_isolat_PS.fasta', 'browser files/COVID-300.lfold.chain',300,'t-RNA/Mg_1M7_1_29903.map', 1.9,-0.7)
# res = concatinate_local_fold('t-RNA/seq_isolat_PS.fasta', 'browser files/COVID-300-noshape.lfold.chain',300)
# res = concatinate_local_fold('t-RNA/sequence.fasta', 't-RNA/sequence-out.lfold.chain',20)
# res = window_fold('cov.fasta',200,'output')
