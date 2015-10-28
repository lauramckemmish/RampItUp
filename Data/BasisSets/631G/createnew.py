import shutil

for i in range(1,2):
    out = '';
    nobf = 0;
    x = open('Basis_reduced/{0}.txt'.format(i))
    y = x.readline()
    while y:
        y = x.readline()
        z = y.split()
        
        if (len(z) > 0 and z[0] == 'S'):
            nobf = nobf+1;
            noprim = z[1]
            out = out + '{0} \n 1 \n\n'.format(noprim)
            for j in range(int(noprim)):
                y = x.readline()
                z=y.split()
                out = out + '{0}d0 \n{1}d0\n'.format(z[0],z[1])
            out = out + '\n'
        elif (len(z) > 1 and z[0] == 'SP'):
            nobf = nobf+4
            noprim = z[1]
            outs = '';
            outp = '';
            for j in range(int(noprim)):
                y=x.readline()
                z=y.split()
                outs = outs + '{0}d0 \n {1}d0\n'.format(z[0],z[1])
                outp = outp + '{0}d0 \n {1}d0\n'.format(z[0],z[2])
            out = out + '{0} \n 1 \n\n'.format(noprim)
            out = out + outs + '\n'
            for j in range(2,5):
                out = out + '{0} \n {1} \n\n'.format(noprim,j)+outp+'\n'

    h = open('Basis_fortran/{0}.bas'.format(i),'w');
    h.write('{0} \n \n'.format(nobf))
    h.write(out)
    h.close()
    
for i in range(3,11):
    out = '';
    nobf = 0;
    x = open('Basis_reduced/{0}.txt'.format(i))
    y = x.readline()
    while y:
        y = x.readline()
        z = y.split()
        
        if (len(z) > 0 and z[0] == 'S'):
            nobf = nobf+1;
            noprim = z[1]
            out = out + '{0} \n 1 \n\n'.format(noprim)
            for j in range(int(noprim)):
                y = x.readline()
                z=y.split()
                out = out + '{0}d0 \n{1}d0\n'.format(z[0],z[1])
            out = out + '\n'
        elif (len(z) > 1 and z[0] == 'SP'):
            nobf = nobf+4
            noprim = z[1]
            outs = '';
            outp = '';
            for j in range(int(noprim)):
                y=x.readline()
                z=y.split()
                outs = outs + '{0}d0 \n {1}d0\n'.format(z[0],z[1])
                outp = outp + '{0}d0 \n {1}d0\n'.format(z[0],z[2])
            out = out + '{0} \n 1 \n\n'.format(noprim)
            out = out + outs + '\n'
            for j in range(2,5):
                out = out + '{0} \n {1} \n\n'.format(noprim,j)+outp+'\n'

    h = open('Basis_fortran/{0}.bas'.format(i),'w');
    h.write('{0} \n \n'.format(nobf))
    h.write(out)
    h.close()
    
        
