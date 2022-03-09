# with open("t-RNA/cov.fasta") as fh:
#     sequence = fh.readlines()
# string = ""
# for line in sequence[1:]:
# 	line = line[:-1]
# 	string += line
# print(coal)
# # print(sequence[0])

# with open("cleandis.fasta") as fh:
#     sequence = fh.readlines()
# string = ""
# for line in sequence[1:]:
# 	line = line[:-1]
# 	string += line
# out = open('cleandis.fasta','w')
# out.write(string)
# print(string)

# def check_parentheses(st):
#     flag = True
#     s = []
#     for i in st:
#         # push if opening bracket
#         if i == "(":
#             s.append(i)
#         else:
#             if len(s) > 0:
#                 # check if top of s is pair of current
#                 # element
#                 temp = s[-1]
#                 s.pop()
#                 if i == "(" and temp != ")":
#                     flag = False
#                     break
#             # if stack is empty, not balanced
#             else:
#                 flag = False
#                 break
#     # If stack is not empty after traversal
#     # then not balanced
#     if len(s) > 0:
#         flag = False

#     if flag:
#         print("balanced")
#     else:
#         print("Not balanced")

# # Driver code
# s = "()..(..)"

open_list = ["[","{","("]
close_list = ["]","}",")"]
def check(myStr):
    stack = []
    for i in myStr:
        if i in open_list:
            stack.append(i)
        elif i in close_list:
            pos = close_list.index(i)
            if ((len(stack) > 0) and
                (open_list[pos] == stack[len(stack)-1])):
                stack.pop()
            else:
                return "Unbalanced"
    if len(stack) == 0:
        return "Balanced"
    else:
        return "Unbalanced"
  
  
# Driver code
# with open('browser files/example files for the SARS structure genome browser/1M7_Mg-L300-D-full_genome.lfold.chain') as fh:
with open('t-RNA/COVID-200.lfold.chain') as fh:
    sequence = fh.readlines()
string = sequence[1]
print("len seq", len(sequence[1]))

print("len paran", len(sequence[2]))
print(check(string))