import RNA
import random

def getShapeDataFromFile(filepath):
    retVec = []
    retVec.append(-999.0);  # data list is 1-based, so we add smth. at pos 0
    count=1
    with open(filepath, 'r') as f:
        lines = f.readlines()
    for line in lines:
        print(line)
        pos = int(line.split()[0])
        value = float(line.split()[1])

        if(pos==count):
            retVec.append(value)
        else:
            for i in range(pos-count):
                retVec.append(-999.0)
            retVec.append(value)
            count=pos
        count+=1
    return retVec

# A positive slope m penalizes high reactivities in paired regions, 
# while a negative intercept b results in a confirmatory `‘bonus’' 
# free energy for correctly predicted base pairs
def predict_window_shape(seq_file,shape_file,output_file,size,m,b):
    shape = getShapeDataFromFile(shape_file)
    print(shape)
    with open(seq_file) as fh:
        sequence = fh.read()
    #
    sequence = 'AUUAAAGUUAUACCUUCCCAGGUAACAAACCAACCAACUUUCGAUCUCUUGUAGAUCUGUUCUCUAAACGAACUUUAAAAUCUGUGUGGCUGUCACUCGGCUGCAUGCUUAGUGCACUCACGCAGUAUAAUUAAUAACUAAUUACUGUCGUUGACAGGACACGAGUAACUCGUCUAUCUUCUGCAGGCUGCUUACGGU'
    #
    md = RNA.md()
    md.window_size = size
    md.max_bp_span = size
    fc = RNA.fold_compound(sequence,md,RNA.OPTION_WINDOW)
    fc.sc_add_SHAPE_deigan(shape,m,b)  
    with open(output_file,"w") as fh:
        mfe = fc.mfe_window(fh)
    file = open(output_file, mode = 'r', encoding = 'utf-8-sig')
    return file

z = predict_window_shape('t-RNA/shape-seq.fasta','t-RNA/shape-react.fasta','t-RNA/output-shape.fasta',50,1,-1)
print(z)


def initialize_fold(seq_file, output_file, size):
    with open(seq_file) as fh:
        sequence = fh.read()
    #
    sequence = 'AUUAAAGUUAUACCUUCCCAGGUAACAAACCAACCAACUUUCGAUCUCUUGUAGAUCUGUUCUCUAAACGAACUUUAAAAUCUGUGUGGCUGUCACUCGGCUGCAUGCUUAGUGCACUCACGCAGUAUAAUUAAUAACUAAUUACUGUCGUUGACAGGACACGAGUAACUCGUCUAUCUUCUGCAGGCUGCUUACGGU'
    #
    ## preferable, more flexible method: local folding via fold_compound interface
    md = RNA.md()
    md.window_size = size
    md.max_bp_span = size
    fc = RNA.fold_compound(sequence,md,RNA.OPTION_WINDOW)

    with open(output_file,"w") as fh:
        mfe = fc.mfe_window(fh)
    file = open(output_file, mode = 'r', encoding = 'utf-8-sig')
    return file

x = initialize_fold('t-RNA/shape-seq.fasta','t-RNA/output-no-shape.fasta',50)
print(x)