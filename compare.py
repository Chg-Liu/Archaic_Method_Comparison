import sys

def readgt():
    gt = {'tsk_'+str(i)+'_'+str(j):[] for i in range(100,200) for j in range(1,3)}
    f = open('../data/'+model+'_'+rep+suffix+'_chr1.ancestry','r')
    start = '0'
    end = '0'
    hap = '0'
    while True:
        a = f.readline().split()
        if not a:
            break
        if a[0] not in gt.keys():
            continue
        if a[0] != hap:
            if hap != '0':
                if end != '0':
                    gt[hap].append((start,end))
            start = '0'
            end = '0'
            hap = a[0]
        if a[-1] == 'Archaic':
            if int(a[1]) > int(end):
                if end != '0':
                    gt[hap].append((start, end))
                start = a[1]
                end = a[2]
            else:
                end = a[2]
    if end != '0': 
        gt[hap].append((start,end))
    f.close()

    return gt

def readhmmix():
    hmmix = {}
    for i in range(100,200):
        for j in range(1,3):
            hmmix['tsk_'+str(i)+'_'+str(j)] = []
            f = open('../'+model+'_'+rep+suffix+'/HMMix/hmmix/'+model+'_'+rep+suffix+'_tsk_'+str(i)+'.decoded-admixpop.hap'+str(j)+'.txt','r')
            f.readline()
            while True:
                a = f.readline().split()
                if not a:
                    break
                if a[4] == 'Archaic':
                    hmmix['tsk_'+str(i)+'_'+str(j)].append((a[1],a[2]))
            f.close()
    return hmmix

def readsprime():
    sprime = {}
    for i in range(100,200):
        for j in range(1,3):
            sprime['tsk_'+str(i)+'_'+str(j)] = []
    f = open('../'+model+'_'+rep+suffix+'/Sprime/Sprime/EUR_EAS_Sprime_phased_tracts_perSample_maxgap0.bed','r')
    while True:
        a = f.readline().split()
        if not a:
            break
        ind = a[3].split(';')[1].split('=')[1]
        if len(sprime[ind])>0:
            if int(sprime[ind][-1][1])>int(a[1]):
                sprime[ind] = sprime[ind][:-1] +[(str(min(int(sprime[ind][-1][0]),int(a[1]))),str(max(int(sprime[ind][-1][1]),int(a[2]))))]
            else:
                sprime[ind].append((a[1],a[2]))
        else:
            sprime[ind].append((a[1],a[2]))
    f.close()
    return sprime

def readas2():
    as2 = {}
    for i in range(100,200):
        for j in range(1,3):
            as2['tsk_'+str(i)+'_'+str(j)] = []
    f = open('../'+model+'_'+rep+suffix+'/ArchaicSeeker2/ArchaicSeeker/AS2_EUR_EAS.seg','r')
    f.readline()
    while True:
        a = f.readline().split()
        if not a:
            break
        as2[a[0]].append((a[2],a[3]))
    f.close()
    return as2

def compare_list(x,y):
    ol_pairs = []
    if x == []:
        for j in y:
            ol_pairs.append(['/' for j in range(5)]+[j[0], j[1], str(int(j[1])-int(j[0])), '0','0']+['0','0','0'])
        return ol_pairs
    elif y == []:
        for i in x:
            ol_pairs.append([i[0], i[1], str(int(i[1])-int(i[0])), '0','0']+['/' for i in range(5)]+['0','0','0'])
        return ol_pairs
    else:
        sx = {}
        sy = {}
        i = 0
        j = 0
        si = {'len':int(x[i][1])-int(x[i][0]), 'n':0, 'overlap':{}}
        sj = {'len':int(y[j][1])-int(y[j][0]), 'n':0, 'overlap':{}}
        while True:
            if i == len(x) and j == len(y):
                break
            elif i == len(x):
                sy[tuple(y[j])] = sj
                j += 1
                if j != len(y):
                    sj = {'len':int(y[j][1])-int(y[j][0]), 'n':0, 'overlap':{}}
                continue
            elif j == len(y):
                sx[tuple(x[i])] = si
                i += 1
                if i != len(x):
                    si = {'len':int(x[i][1])-int(x[i][0]), 'n':0, 'overlap':{}}
                continue

            if int(x[i][1]) <= int(y[j][0]):
                sx[tuple(x[i])]=si
                i += 1
                if i != len(x):
                    si = {'len':int(x[i][1])-int(x[i][0]), 'n':0, 'overlap':{}}
                continue
            elif int(x[i][1]) < int(y[j][1]):
                ollen = min(int(x[i][1]),int(y[j][1]))-max(int(x[i][0]), int(y[j][0]))
                si['overlap'][tuple(y[j])] = ollen
                si['n'] += ollen
                sj['overlap'][tuple(x[i])] = ollen
                sj['n'] += ollen
                sx[tuple(x[i])]=si
                i += 1
                if i!= len(x):
                    si = {'len':int(x[i][1])-int(x[i][0]), 'n':0, 'overlap':{}}
                continue
            elif int(x[i][1]) == int(y[j][1]):
                ollen = min(int(x[i][1]),int(y[j][1]))-max(int(x[i][0]), int(y[j][0]))
                si['overlap'][tuple(y[j])] = ollen
                si['n'] += ollen
                sj['overlap'][tuple(x[i])] = ollen
                sj['n'] += ollen
                sx[tuple(x[i])]=si
                sy[tuple(y[j])]=sj
                i += 1
                j += 1
                if i != len(x):
                    si = {'len':int(x[i][1])-int(x[i][0]), 'n':0, 'overlap':{}}
                if j != len(y):
                    sj = {'len':int(y[j][1])-int(y[j][0]), 'n':0, 'overlap':{}}
                continue
            elif int(x[i][0]) >= int(y[j][1]):
                sy[tuple(y[j])]=sj
                j += 1
                if j != len(y):
                    sj = {'len':int(y[j][1])-int(y[j][0]), 'n':0, 'overlap':{}}
            else:
                ollen = min(int(x[i][1]),int(y[j][1]))-max(int(x[i][0]), int(y[j][0]))
                si['overlap'][tuple(y[j])] = ollen
                si['n'] += ollen
                sj['overlap'][tuple(x[i])] = ollen
                sj['n'] += ollen
                sy[tuple(y[j])]=sj
                j += 1
                if j != len(y):
                    sj = {'len':int(y[j][1])-int(y[j][0]), 'n':0, 'overlap':{}}
        for segx in sx:
#            print(segx,sx[segx])
            sx[segx]['percentage'] = round(sx[segx]['n']/sx[segx]['len'],4)
        for segy in sy:
#            print(segy,sy[segy])
            sy[segy]['percentage'] = round(sy[segy]['n']/sy[segy]['len'],4)
        j = 0
        for seg in x:
            if j<len(y):
                while int(y[j][1]) <= int(seg[0]) and sy[y[j]]['overlap'] == {}:
                    ol_pairs.append(['/' for i in range(5)]+[y[j][0], y[j][1], str(sy[y[j]]['len']), '0','0']+['0','0','0'])
                    j += 1
                    if j == len(y):
                        break
                
            if sx[seg]['overlap'] == {}:
                ol_pairs.append([seg[0], seg[1], str(sx[seg]['len']), '0','0']+['/' for i in range(5)]+['0','0','0'])

            else:
                for k,m in sx[seg]['overlap'].items():
                    ol_pairs.append([seg[0], seg[1], str(sx[seg]['len']), str(sx[seg]['n']), str(sx[seg]['percentage'])] + [k[0], k[1], str(sy[k]['len']), str(sy[k]['n']), str(sy[k]['percentage'])] + [str(m), str(round(m/sx[seg]['len'],4)), str(round(m/sy[k]['len'],4))])
#                print(list(sx[seg]['overlap'].keys()))
                if list(sy[list(sx[seg]['overlap'].keys())[-1]]['overlap'].keys())[-1] == seg:
                    j += 1
        return ol_pairs

def compare_to_gt(method, result):
    f = open(model+'_'+rep+suffix+'.groundtruth_'+method+'.tsv','w')
    for hap in gt.keys():
        print(hap,'analyzing')
        print(gt[hap])
        print(result[hap])
        ol = compare_list(gt[hap], result[hap])
        f.writelines(['\t'.join([hap]+line)+'\n' for line in ol])
        print(hap,'written')
    f.close()
    

if __name__ == '__main__':
    model = sys.argv[1]
    rep = sys.argv[2]
    suffix = sys.argv[3]
    if suffix == '0':
        suffix = ''

    gt = readgt()
    sprime = readsprime()
#    hmmix = readhmmix()
#    as2 = readas2()
    compare_to_gt('sprime',sprime)
#    compare_to_gt('hmmix',hmmix)
#    compare_to_gt('as2',as2)

