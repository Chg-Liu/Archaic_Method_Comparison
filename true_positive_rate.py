def readresult(x,m):
    A = {}
    total = 0
    dlen = 0
    tlen = 0
    seg = {}
    f = open(x,'r')
    while True:
        a = f.readline().split()
        if not a:
            break
        if '/' in a[:3]:
            continue
        if tuple(a[:3]) not in seg.keys():
            seg[tuple(a[:3])] = 0
            total += 1
            for i in range(0,110,10):
                if i not in A.keys():
                    A[i] = 0
                if i == 0:
                    if float(a[5]) > i/100:
                        A[i] += 1
                elif float(a[5])>=i/100:
                    A[i]+=1
            dlen += int(a[4])
            tlen += int(a[3])

    f.close()
    ratelist = []
    for i in A.keys():
        if i == 0:
            label = '>0%'
        elif i==100:
            label = '=100%'
        else:
            label = '>='+str(i)+'%'
        ratelist.append([str(round(A[i]/total,4)),label,m])
    ratelist.append([str(round(dlen/tlen,4)),'length',m])
    return ratelist

if __name__ == '__main__':
    fout = open('basic_true_positive_rate.tsv','w')
    fout.write('\t'.join(['TPRate','Category','Method'])+'\n')
    for i in range(1,21):
        k = readresult('basic_'+str(i)+'.groundtruth_sprime.tsv','SPrime')
        fout.writelines(['\t'.join(j)+'\n' for j in k])
        k = readresult('basic_'+str(i)+'.groundtruth_hmmix.tsv','HMMix')
        fout.writelines(['\t'.join(j)+'\n' for j in k])
    fout.close()

        



            
            


