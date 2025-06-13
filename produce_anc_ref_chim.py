import gzip,os,sys

prefix = sys.argv[1]
if sys.argv[2] == '0':
    suffix = ''
else:
    suffix = sys.argv[2]
length = int(sys.argv[3])
rep = sys.argv[4]

na = prefix+'_'+str(rep)+suffix
f = gzip.open(na+'_chr1.vcf.gz','rt')
while True:
    a = f.readline().split()
    if a[0][0] == '#' and a[0][1] != '#':
        break

k = 1
fout1 = open(na+'_homo_sapiens_ancestor_1.fa','w')
fout1.write('>ANCESTOR_for_chromosome:GRCh37:1:1:'+str(length)+':1\n')
fout2 = open(na+'_refgenome_chr1.fa','w')
fout2.write('>chr1\n')
m1=0
m2=0

a = f.readline().split()
while k <= length:
    while k<int(a[1]):
        fout1.write('A')
        fout2.write('A')
        if k%100 == 0:
            fout1.write('\n')
        if k%50 == 0:
            fout2.write('\n')
        k += 1
    if k>length:
        break
    fout1.write(a[3])
    fout2.write(a[3])
    if k%100 == 0:
        fout1.write('\n')
    if k%50 == 0:
        fout2.write('\n')
    a = f.readline().split()
    if not a:
        a = ['1',str(length+1)]
    k += 1
f.close()
fout1.close()
fout2.close()
os.system('bgzip -k '+na+'_homo_sapiens_ancestor_1.fa')

na = prefix+'_'+str(rep)+suffix
f = gzip.open('Chim_'+na+'_chr1.vcf.gz','rt')
while True:
    a = f.readline().split()
    if a[0][0] == '#' and a[0][1] != '#':
        break

k = 1
fout3 = open('chr1.'+na+'.chimp.fa','w')
fout3.write('>chr1\n')
m1=0
m2=0

a = f.readline().split()
while k <= length:
    while k<int(a[1]):
        fout3.write('A')
        if k%50 == 0:
            fout3.write('\n')
        k += 1
    if k>length:
        break
    if a[9][0] == '0':
        fout3.write(a[3])
    elif a[9][0] == '1':
        fout3.write(a[4])
    if k%50 == 0:
        fout3.write('\n')
    a = f.readline().split()
    if not a:
        a = ['1',str(length+1)]
    k += 1
f.close()
fout3.close()
os.system('bgzip chr1.'+na+'.chimp.fa')


