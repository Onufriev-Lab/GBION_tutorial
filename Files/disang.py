import math
namedisang = "disang_NaCl.txt"
Number_of_atoms_for_com = 10
number_of_ions=10316
starting_index_of_ion=24942
sum_x = 0
sum_y = 0
sum_z = 0
dist = []
xval = []
yval = []
zval = []
atom_index = []
with open('nucleosome.pdb', 'r') as fin:
    for line in fin:
        if line.startswith('ATOM'):
            line_vals = line.strip().split()
            for i in range(2, len(line_vals)):
                if line_vals[i] == '1':
                    index_in_pdb = i
            xval.append(float(line_vals[index_in_pdb + 1]))
            yval.append(float(line_vals[index_in_pdb + 2]))
            zval.append(float(line_vals[index_in_pdb + 3]))
            atom_index.append(int(line_vals[1]))
if len(xval) < Number_of_atoms_for_com:
    Number_of_atoms_for_com = len(xval)
for i in range(len(xval)):
    sum_x += xval[i]
    sum_y += yval[i]
    sum_z += zval[i]
com_x = sum_x/len(xval)
com_y = sum_y/len(yval)
com_z = sum_z/len(zval)
for i in range(len(xval)):
    dist.append( math.sqrt( ( xval[i] - com_x ) ** 2 + ( yval[i] - com_y ) ** 2 + ( zval[i] - com_z ) ** 2 ) )
min_dist = dist[0]
atoms_for_com = [1] * Number_of_atoms_for_com
for i in range(len(dist)):
    if dist[i] < min_dist:
        for j in range(0, len(atoms_for_com) - 1):
            atoms_for_com[len(atoms_for_com) - j - 1] = atoms_for_com[len(atoms_for_com) - j - 2]
        atoms_for_com[0] = atom_index[i]
        min_dist = dist[i]
with open(namedisang,'w') as disang:
    for i in range(number_of_ions):
        id_ion=starting_index_of_ion + i
        disang.write("&rst\n")
        disang.write("iresid=0,\n")
        disang.write("  iat=-1,-1,r1=0.0,r2=0.0,r3=240,r4=250,rk2=0.0, rk3=20.0,\n")
        disang.write("igr1=")
        for j in range(len(atoms_for_com)):
            disang.write(str(atoms_for_com[j]))
            disang.write(",")
        disang.write("igr2={}\n".format(id_ion))
        disang.write(" /\n")
