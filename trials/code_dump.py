
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

#get SHAPE data from a .map file
# def getShapeDataFromFile(shape_file):
#     retVec = []
#     retVec.append(-999.0);  # data list is 1-based, so we add smth. at pos 0
#     count = 1
#     with open(shape_file, 'r') as f:
#         lines = f.readlines()
#     for line in lines:
#         line = line.split()
#         pos = int(line[0])
#         value = float(line[1])
#         if(pos==count):
#             retVec.append(value)
#         else:
#             for i in range(pos-count):
#                 retVec.append(-999.0)
#             retVec.append(value)
#             count=pos
#         count += 1
#     return retVec


# def local_fold_list(input_file, output_file, size, shape_file = None, m=None, b=None, hc= None, sc = None):
#     # Create and open file
#     if len(shape_file) == 2:
#         file1, file2, sequence = initialize_fold(input_file, output_file, size, shape_file, m, b) #the function returns (file1,file2,sequence)
#         file2 = open(file2, "r")
#         data2 = file2.read()
#         file2.close()
#         file1 = open(file1, "r")
#         data1 = file1.read()
#         file1.close()
#         file = open('mix.fasta', "a")
#         file.write(data2)
#         file.write(data1)
#         file.close()
#         file = open('mix.fasta', mode = 'r', encoding = 'utf-8-sig')
#     else:
#         file = initialize_fold(input_file, output_file, size, shape_file, m, b)[0] #the function returns (file,sequence)
#         file = open(file, mode = 'r', encoding = 'utf-8-sig')
        
#     # Parse the t-RNA Sequence
#     lines = file.readlines()
#     d = {}
#     for line in lines:
#         window = line.split()
#         #window:  ['.(((((.........))))).', '(', '-0.10)', '50']
#         # dictionary with {50 : '.(((((.........))))).' , etc.}
#         if window != []:
#             key = int(window[-1])
#             d[key] =  window[0]
#     # for key in sorted(d.keys()) :
#     #     print(key , " :: " , d[key])
#     return d