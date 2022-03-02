# Parse the t-RNA Sequence
file = open('covid.FASTA', mode = 'r', encoding = 'utf-8-sig')
lines = file.readlines()
seq = []
for line in lines[1:-2]:
    window = line.split()
    print(window)
    #window:  ['.(((((.........))))).', '(', '-0.10)', '50']
    # To remove the '-0.10)' parentheses at the end
    window[2] = window[2][:-1]
    # To remove the '(' member of the list
    window.remove('(')
    seq.append(window)
print("sequence result: \n",seq,"\n")

#intervals = []
#for window in seq:
#    (a,b,c) = 'L'+window[2], 'R'+str(int(window[2])+len(window[0])), float(window[1])
#    intervals.append((a,b,c))
#intervals.reverse()
#for i in range(len(intervals)):
#    I.append(str(i)+intervals[i][0])
#    I.append(str(i)+intervals[i][1])
#############################################################
#
# Makes I as a list of tuples: 
# I = [(index on the original sequence, start or end of interval (L or R), interval no. )]
#
#############################################################
intervals = []
for window in seq:
    (a,b,c) = (int(window[2]),'L'), (int(window[2])+len(window[0]),'R'), float(window[1])
    # a  = left end (start)
    # b = right end (finish)
    # c = energy (value)
    intervals.append((a,b,c))
    
intervals.reverse()
I = []
for i in range(len(intervals)):
    a, b = intervals[i][0]
    c, d = intervals[i][1]
    # a and c are the location start/end
    # b and d are the Left or Right 
    # i is the index of the tuples from 'Interval'
    I.append((a,b,i))
    I.append((c,d,i))

I.sort()
print("intervals: \n",intervals,"\n")
print("I: \n",I, "\n")

def pop_until_j(l, j):
    if j not in l: return
    while (len(l) > 0):
        current = l.pop()
        if current == j:
            
            l.append(current)
            break
#
#minimum = 0
#V = [ j[2] for j in intervals]
#v = V.copy()
#tracker = []
#t = []
#for i in range(len(I)):
#    # If Left position of interval j, set V[j] = v[j] + minimum
#    if I[i][1] == 'L':
#        j = I[i][2]
#        V[j] = minimum + v[j]
#       
#    # If Right position of interval j, set minimum = min(minimum, V[j])
#    if I[i][1] == 'R':
#        j = I[i][2]
#        if V[j] <= minimum:
#            minimum = V[j]
#            tracker.append(j)
#        
#    print("step:", i, V,"min:",minimum )
#
##print("V = \n", V,"\n")      
minimum = 0
min_int_temp = -1
V = [ (j[2],[j]) for j in intervals]
v = [ j[2] for j in intervals]

for i in range(len(I)):
    # If Left position of interval j, set V[j] = v[j] + minimum
    if I[i][1] == 'L':
        j = I[i][2]
        V_0 = minimum + v[j]
        if min_int_temp != -1:
            V_1 = V[min_int_temp][1]+[intervals[j]]
        else:
            V_1 = [intervals[j]]
        V[j] = (V_0,V_1)
       
    # If Right position of interval j, set minimum = min(minimum, V[j])
    if I[i][1] == 'R':
        j = I[i][2]
        if V[j][0] <= minimum:
            minimum = V[j][0]
            min_int_temp = j
        
#    print("step:", i, V,"min:",minimum )

#print("V = \n", V,"\n")      
print("intervals for minimum: ", V[min_int_temp][1])  
print("minimum:",minimum, "\n")  
