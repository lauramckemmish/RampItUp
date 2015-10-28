import shutil

x = open('moleculename.csv')
basename = x.readline()
basename = basename.strip('\n')
basename = basename.strip(' ')
x.close()

job = 'allgaus.inp'
x = open('{0}'.format(job))
f=x.readline()
swit=0
zz = ''
atid = 0;
while f:
    if '$molecule' in f:
        swit = 1
        f = x.readline()
        g = f.split()
        charge = g[0]
        noelec = -int(charge);
        multiplicity = g[1]
        f = x.readline()
    if '$end' in f:
        swit = 0
    if swit == 1:
        atid = atid+1
        h = f.split()
        noelec = noelec+int(h[0])
        j='\n'.join(h)
        zz = zz  + str(atid) + '\n'+j + '\n'
    f=x.readline()
#print(k)
x.close()
y = open('geom.csv'.format(basename),'w')
z2 = str(atid) + '\n' +str(noelec) + '\n' +multiplicity + '\n' +zz
y.write(z2)
y.close()

job ='allgaus' + '.out'
x = open('{0}'.format(job))
print('6-31G All-Gaussian calculation')
f = x.readline()
swit = 0
while f:
    if 'Current DIIS' in f:
        swit = 1
    if (swit == 1) and 'Total number of integrals generated' in f:
        print(f)
        f = x.readline()
        print(f)
    f = x.readline()
x.close()

job ='allgaus' + '.out'
x = open('{0}'.format(job))
#print('In All-Gaussian calculation')
f = x.readline()
swit = 0
while f:
    if 'Convergence criterion met' in f:
        swit = 1;
   # if 'One-Electron    Energy' in f and swit ==1:
    #    print(f)
    #if 'Nuclear Repu.   Energy ' in f and swit ==1:
    #    print(f)
    if 'Total energy in the final basis set =' in f and swit ==1:
     #   print(f)
        g = f.split()
        gausenergy = g[len(g)-1]
    f = x.readline()
x.close()

y=open('{0}qchemenergies.out'.format(basename),"w")
y.write(gausenergy)
y.close()


job = 'allgaus.inp'
x = open('{0}'.format(job))
f=x.readline()
swit=0
zz = ''
atid = 0;
while f:
    if '$molecule' in f:
        swit = 1
        f = x.readline()
        g = f.split()
        charge = g[0]
        noelec = -int(charge);
        multiplicity = g[1]
        f = x.readline()
    if '$end' in f:
        swit = 0
    if swit == 1:
        atid = atid+1
        h = f.split()
        noelec = noelec+int(h[0])
        j='\n'.join(h)
        zz = zz  + str(atid) + '\n'+j + '\n'
    f=x.readline()
#print(k)
x.close()
y = open('geom.csv'.format(basename),'w')
z2 = str(atid) + '\n' +str(noelec) + '\n' +multiplicity + '\n' +zz
y.write(z2)
y.close()


