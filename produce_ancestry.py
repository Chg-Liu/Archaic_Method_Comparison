import sys

def sampleID(x):
    return 'tsk_'+str(int(x/2))+'_'+str(x%2+1)

if __name__ == '__main__':

    Anc = {'2':'Chim','3':'Archaic','5':'AFR','6':'EUR_EAS'}

    f = open(sys.argv[1],'r')
    A = [i.split() for i in f.readlines()[1:]]
    f.close()

    fout = open(sys.argv[2],'w')
    fout.write('\t'.join(['sample','left','right','ancestry'])+'\n')

    for i in A:
        fout.write('\t'.join([sampleID(int(i[0])), str(int(float(i[1]))), str(int(float(i[2]))), Anc[i[-1]]])+'\n')

    fout.close()
    

