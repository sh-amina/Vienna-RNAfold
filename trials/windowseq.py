import RNA
import os
import subprocess

#list_dir = subprocess.Popen(["ls", "-l"])
#list_dir.wait()


#out = subprocess.run(["RNALfold", "-L25"],text=True , input= "GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACGGUUACGUAGGCAAUGGCCAAUUUGCCAA",stdout=subprocess.PIPE)
#out = os.system("RNALfold -L25 < sequence.fasta")
#print("********")
#print(subprocess.PIPE)

seq = "GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACGGUUACGUAGGCAAUGGCCAAUUUGCCAA"
l=[]

# create a fold_compound object for the current sequence
#fc = RNA.fold_compound(seq)
#fc.mfe_window(25)
fc = RNA.Lfold(seq,25)
# compute the MFE and corresponding structure
(mfe_struct, mfe) = fc.mfe()
#print("********")
#print(type(fc))
#print(l)