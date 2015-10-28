import shutil

basename = ''
z = open('moleculename.csv')
basename = z.readline()
basename = basename.strip()

z = open('basisset.csv')
basisset = z.readline()
basisset = basisset.strip()

z = ''
job = basename + '.inp'
x = open('{0}'.format(job))
print('Constructing input files for {0} from base file'.format(basename))
f = x.readline()
swit = 0
zpt = ''
zpo = ''
zo = ''
zt = ''
zs = ''
za = ''
atoms = []
while f:
    if '$molecule' in f:
        swit = 3;
    g = f.split()
    if swit==3 and len(g) > 2:
        atoms.append(int(g[0])) 
    if '$end' in f:
        swit = 0;
    if '$rem' in f:
        swit = 1;   
        zpt = zpt + f
        zpo = zpo + f
        zo = zo +f
        za = za + f + 'pathn = 20 \n DIRECT_SCF = true \n max_SCF_CYCLES = 50 \n THRESH = 10 \n'
        zs = zs + f + 'pathn = 20 \n DIRECT_SCF = true \n max_SCF_CYCLES = 50 \n THRESH = 10 \n'
        zt = zt + f
        zo = zo + '   symmetry = off\n   iprint = 200 \n   pathn = 10\n'
        zt = zt + '   symmetry = off\n   aoints_debug_print = 4 \n   pathn = 10 \n '
        zpo = zpo + '   symmetry = off\n   iprint = 200 \n   pathn = 10 \n'  
        zpt = zpt + '   symmetry = off\n   aoints_debug_print = 4 \n   pathn = 10 \n'  
        f = x.readline()
    if 'BASIS' in f or 'basis' in f: 
        f = '   BASIS  =  general \n'
    	zpt = zpt + f
    	zpo = zpo + f
    	zo = zo +f
    	zt = zt + f
    	zs 	= zs + f + 'SCF_CONVERGENCE = 2 \n'
    	if 'pG' not in basisset:
             za = za + 'BASIS = 6-31G \n SCF_CONVERGENCE = 6 \n SCF_FINAL_PRINT = 2 \n'
        else:
             za = za + 'BASIS = 6-31+G \n SCF_CONVERGENCE = 6 \n SCF_FINAL_PRINT = 2 \n'
    else:
	zpt = zpt + f
   	zpo = zpo + f
    	zo  = zo + f
    	zt = zt + f
    	zs = zs + f
    	za = za + f
    f = x.readline()
x.close()

# Add the basis function section

atoms = list(set(atoms))
output2 = ''
output2=output2+"\n \n$basis \n"
for atom in atoms:
    y3=open('{1}/Basis_pseudo/{0}.txt'.format(atom,basisset))
    f2=y3.readline()
    while f2:
        output2=output2+f2
        f2=y3.readline()
    y3.close()
    output2 = output2 
zpo = zpo + output2 + '$end' 
zpt = zpt + output2 + '$end'

output2 = ''
output2=output2+"\n \n$basis \n"
for atom in atoms:
    y3=open('{1}/Basis_reduced/{0}.txt'.format(atom,basisset))
    f2=y3.readline()
    while f2:
        output2=output2+f2
        f2=y3.readline()
    y3.close()
    output2 = output2 
zt = zt + output2+ '$end'
zo = zo + output2+ '$end'
zs = zs + output2 + '$end'

jobout = 'allgaus.inp'
y = open('{0}'.format(jobout),'w')
y.write(za)
y.close()


