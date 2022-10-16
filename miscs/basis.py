import symfem
import re

p_qua = symfem.create_element("quadrilateral", "Lagrange", 4)
# p_qua.plot_basis_function(0, "/home/dsiedel/projetcs/lolita/miscs/basis.svg")
p_qua.plot_dof_diagram("/home/dsiedel/projetcs/lolita/miscs/dofs.svg")
for fun in p_qua.map_to_cell(((-1, -1), (1, -1), (-1, 1), (1, 1))):
    print("---")
    # print(fun)
    count = 0
    terms = []
    flag = False
    cc = str(fun)
    term = ""
    s_term = "+ "
    for i, c in enumerate(cc):
        term += c
        if cc[i] == '(':
            flag = True
        if cc[i] == ')':
            flag = False
        if (cc[i] == ' ' and cc[i - 1] in ['+', '-'] and cc[i - 2] == ' ' and not flag) or i == len(cc) - 1:
            terms.append(term)
            print(term)
            # print(re.sub(r'\((.*?)\)\*\*(.)', r'std::pow(\1, \2)', term))
            # print(re.findall(r"\d+", term))
            c = re.sub(r'\((.*?)\)\*\*(/[1-9][0-9]*(?:\/[1-9][0-9])*/g)', r'std::pow(\1, \2)', term)
            c = re.sub(r"(.\d+)", r'\1.0', c)
            c = re.sub(r'x', r'basis_vector(0)', c)
            c = re.sub(r'y', r'basis_vector(1)', c)
            c = re.sub(r'\*', r' * ', c)
            c = re.sub(r'\/', r' / ', c)
            # print(s_term + c)
            if i == len(cc) - 1:
                print(s_term + c)
            else:
                print(s_term + c[:-2])
            s_term = "{} ".format(cc[i - 1])
            term = ""
    # for t in terms:
    #     print(t)
t = p_qua.map_to_cell(((-1, -1), (1, -1), (-1, 1), (1, 1)))[0]
# for c in t:
# print(t.grad(1))

# p_tri = symfem.create_element("triangle", "P", 1)
# for fun in p_tri.map_to_cell(((-1, -1), (1, -1), (-1, 1))):
#     print(fun)