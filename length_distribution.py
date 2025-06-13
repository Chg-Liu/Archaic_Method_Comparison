def readseg(x,rep):
    model = prefix+'_'+str(rep)+suffix
    seggt = {}
    seginfer = {}
    gt = []
    infer = []
    f = open(model+'.groundtruth_'+x+'.tsv','r')
    while True:
        a = f.readline().split() 
        if not a:
            break
        if tuple(a[:3]) not in seggt.keys() and a[1] !='/':
            gt.append(a[3])
            seggt[tuple(a[:3])] = 0
        if tuple([a[0],a[6],a[7]]) not in seginfer.keys() and a[6]!= '/':
            infer.append(a[8])
            seginfer[tuple([a[0],a[6],a[7]])] = 0
    f.close()
    return gt,infer

if __name__ == '__main__':
    prefix = 'basic'
    suffix = ''
    md = {'sprime':'SPrime','hmmix':'HMMix'}
    fout = open(prefix+'_length_distribution.tsv','w')
    fout.write('\t'.join(['length','Method','Replicate'])+'\n')
    for m in ['sprime','hmmix']:
        for i in range(1,21):
            gt,infer = readseg(m,i)
            if m == 'sprime':
                fout.writelines(['\t'.join([k,'GroundTruth','GroundTruth_'+str(i)])+'\n' for k in gt])
            fout.writelines(['\t'.join([k,md[m],md[m]+'_'+str(i)])+'\n' for k in infer])
    fout.close()



