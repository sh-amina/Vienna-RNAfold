import RNA
import random

def initialize_fold(input_file, size, output_file, hc = None, sc = None ):
    ## preferable, more flexible method: local folding via fold_compound interface
    md = RNA.md()
    md.window_size = size
    md.max_bp_span = size
    with open(input_file) as fh:
        sequence = fh.read()
    fc = RNA.fold_compound(sequence,md,RNA.OPTION_WINDOW)
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
    return file


def window_fold(input_file, size, output_file):
    # Create and open file
    file = initialize_fold(input_file, size, output_file)
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
    return result   


def local_fold_list(input_file, size, output_file):
    # Create and open file
    file = initialize_fold(input_file, size, output_file)
    # Parse the t-RNA Sequence
    lines = file.readlines()
    d = {}
    for line in lines:
        window = line.split()
        #window:  ['.(((((.........))))).', '(', '-0.10)', '50']
        # dictionary with {50 : '.(((((.........))))).' , etc.}
        if window != []:
            d[int(window[-1])] =  window[0]
    return d


def concatinate_local_fold(input_file,size, output_file):
    result = window_fold(input_file, size, output_file)
    # results = [([1, 37], -16.8), ([46, 74], -3.3), ([76, 122], -6.4), ([123, 164], -9.1), ([167, 199], -11.7)]
    D = local_fold_list(input_file, size, output_file)
    # D = {13: '.((((.((((((((((.(.......).))))))........))))..)))).',
    #       2: '.((((..(((((((.....))))))))))).',
    #       1: '(((((.(((((((((.....))))))).))))))).' , etc.}
    seq = ''
    for i in range(len(result)):
        start = result[i][0][0]
        length = len(seq)
        for i in range(start-length-2):
            seq += '.'
        seq += D[start]
    return(seq)

res = concatinate_local_fold('cov.fasta',200,'output')
# res = window_fold('cov.fasta',200,'output')
print(res)
 

# print('Heya')