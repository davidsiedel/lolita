import symfem
import re

p_qua = symfem.create_element("quadrilateral", "S", 3)
for fun in p_qua.map_to_cell(((-1, -1), (1, -1), (-1, 1), (1, 1))):
    print("---")
    # print(fun)
    count = 0
    terms = []
    flag = False
    cc = str(fun)
    term = ""
    for i, c in enumerate(cc):
        term += c
        if cc[i] == '(':
            flag = True
        if cc[i] == ')':
            flag = False
        if (cc[i] == ' ' and cc[i - 1] in ['+', '-'] and cc[i - 2] == ' ' and not flag) or i == len(cc) - 1:
            terms.append(term)
            term = ""
    for t in terms:
        print(t)
t = p_qua.map_to_cell(((-1, -1), (1, -1), (-1, 1), (1, 1)))[0]
# for c in t:
# print(t.grad(1))

# p_tri = symfem.create_element("triangle", "P", 1)
# for fun in p_tri.map_to_cell(((-1, -1), (1, -1), (-1, 1))):
#     print(fun)