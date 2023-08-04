#!/home/user/.pyenv/shims/python3

import sys
args = sys.argv
print(args[0])
print(args[1])
print(args[2])
print(args[3])
print(args[4])

path = args[4] + args[1]
outfile = args[4] + args[1] + '/stat_out1_' + args[1]
outfile2 = args[4] + args[1] + '/stat_out2_' + args[1]
outfileb = args[4] + args[1] + '/stat_out1_' + args[1]
outfileb2 = args[4] + args[1] + '/stat_out2_' + args[1]
outfile3 = args[4] + args[1] + '/stat_out3_' + args[1]
outfile4 = args[4] + args[1] + '/stat_outlog_' + args[1] + '.txt'

import pandas as pd
import numpy as np
from scipy import stats
import re
import statistics
import os
import random
import math
import shutil

pth0 = 0.05
pth1 = 0.001
pth2 = 0.0001
scth0 = 0
scth1 = 0
scth2 = 0
total = 0
pcou1 = 0
pcou2 = 0
pcou3 = 0
pcou4 = 0
pcou5 = 0
pcou6 = 0
pcou1b = 0
pcou2b = 0
pcou3b = 0
pcou4b = 0
pcou5b = 0
pcou6b = 0
corlist1 = []
coufrrf = 0

frx1 = args[2]
frx2 = args[3]

pthb = str(pth1)
scthb = str(scth1)
pth2b = str(pth2)
scth2b = str(scth2)
pth0b = str(pth0)
scth0b = str(scth0)
print('pvalue_threshold_1\t',pthb)
print('score_threshold_1\t',scthb)
print('pvalue_threshold_2\t',pth2b)
print('score_threshold_2\t',scth2b)
print('pvalue_threshold_3\t',pth0b)
print('score_threshold_3\t',scth0b)
pthc = re.sub('\W', '', pthb)
scthc = re.sub('\W', '', scthb)
pth2c = re.sub('\W', '', pth2b)
scth2c = re.sub('\W', '', scth2b)
pth0c = re.sub('\W', '', pth0b)
scth0c = re.sub('\W', '', scth0b)
#outfile = outfile + '_pth' + pthc + '_scth' + scthc + '.txt'
#outfile2 = outfile2 + '_pth' + pthc + '_scth' + scthc + '.txt'
#outfileb = outfileb + '_pth' + pth2c + '_scth' + scth2c + '.txt'
#outfileb2 = outfileb2 + '_pth' + pth2c + '_scth' + scth2c + '.txt'
#outfile3 = outfile3 + '_pth' + pth0c + '_scth' + scth0c + '.txt'
outfile = outfile + '_pth' + pthc + '.txt'
outfile2 = outfile2 + '_pth' + pthc + '.txt'
outfileb = outfileb + '_pth' + pth2c + '.txt'
outfileb2 = outfileb2 + '_pth' + pth2c + '.txt'
outfile3 = outfile3 + '_pth' + pth0c + '.txt'
print('outfile',outfile)
print('outfile2',outfile2)
print('outfileb',outfileb)
print('outfileb2',outfileb2)
print('outfile3',outfile3)
print('outfile4',outfile4)

if os.path.isfile(outfile):
    os.remove(outfile)
if os.path.isfile(outfile2):
    os.remove(outfile2)
if os.path.isfile(outfileb2):
    os.remove(outfileb2)
if os.path.isfile(outfile3):
    os.remove(outfile3)
if os.path.isfile(outfile4):
    os.remove(outfile4)

def robust_z(x):
    from sklearn.preprocessing import robust_scale
    from scipy.stats import norm

    coefficient = norm.ppf(0.75)-norm.ppf(0.25)
    robust_z_score = robust_scale(x)*coefficient

    return robust_z_score

data5 = {}
data5b = ''
data5n = 0
fn = 0

files = os.listdir(path)
for infiletest in files:
    if 'calc_deepliftscore_out' in infiletest:
        print(infiletest)
        tfid = re.sub('calc_deepliftscore_out_', '', infiletest)
        tfid = re.sub('.txt', '', tfid)
        infiletest = path + '/' + infiletest
        with open(infiletest) as f:
            lines = f.read().splitlines()

        k = 0
        infile = []
        infile2 = []
        if len(lines) != 4 and len(lines) != 2:
            continue
        for i in range(0, len(lines)):
            if 'promoter' not in lines[i] and 'enhancer' not in lines[i]:
                continue
            filename = re.sub('DeepLIFT\/.+', "DeepLIFT/", lines[i])
            filename2 = re.sub('DeepLIFT\/.+', "test_data/cor_tbl.txt", lines[i])
            filename = re.sub('.+\/Btrainout.+?\/(.+)', r'\1', filename)
            filename = path + '/' + filename
            filename2 = re.sub('.+\/Btrainout.+?\/(.+)', r'\1', filename2)
            filename2 = path + '/' + filename2
            if os.path.exists(filename2) == 0 or os.path.getsize(filename2) == 0:
                continue
            infile.insert(k,filename)
            infile2.insert(k,filename2)
            k += 1
        if len(infile) < 2 or len(infile2) < 2:
            continue

        for i in range(0,1):
            if i == 0:
                if frx1 in infile[0] and frx2 in infile[1]:
                    score1 = infile[0]
                    score2 = infile[1]
                    cor1 = infile2[0]
                    cor2 = infile2[1]
                else:
                    continue

            with open(outfile, mode='a') as f:
                print(score1,file=f)
                print(score2,file=f)

            data1h = {"A":0}
            data2h = {"A":0}
            data1b = []
            data2b = []
            data12c = []
            for j in range(0,1): # DNA_0.txt.gz
                score1b = score1 + 'DNA_' + str(j) + '.txt.gz'
                score2b = score2 + 'DNA_' + str(j) + '.txt.gz'
                data1 = pd.read_csv(score1b,header=0,sep='\t',compression='gzip')
                data2 = pd.read_csv(score2b,header=0,sep='\t',compression='gzip')

                data1 = data1.values
                data2 = data2.values

                data12 = np.intersect1d(data1[:,1], data2[:,1])
                data12b = data12.tolist()

                for j2 in range(data1.shape[0]):
                        data1c = np.array([float(s) for s in data1[j2,2].split(',')])
                        data1c_indx=np.array(data1c)>0
                        data1c = data1c[data1c_indx] / sum(abs(data1c))
                        data1b.extend(data1c)
                for j2 in range(data2.shape[0]):
                        data2c = np.array([float(s) for s in data2[j2,2].split(',')])
                        data2c_indx=np.array(data2c)>0
                        data2c = data2c[data2c_indx] / sum(abs(data2c))
                        data2b.extend(data2c)

                data1h.clear
                data2h.clear
            if (len(data1b) == 0 or len(data2b) == 0):
                 continue

            corname1 = ''
            corname2 = ''
            corname1b = []
            corname2b = []
            corsc1 = ''
            corsc2 = ''
            corsc1b = []
            corsc2b = []
            with open(cor1) as f:
                 lines = f.read().splitlines()
                 for i in range(1, len(lines)):
                      cor = lines[i].split('\t')
                      corname1 = corname1 + '\t\t' + cor[0]
                      corname1b.append(cor[0])
                      corsc1 = corsc1 + '\t' + cor[1]
                      corsc1b.append(cor[1])
            with open(cor2) as f:
                 lines = f.read().splitlines()
                 for i in range(1, len(lines)):
                      cor = lines[i].split('\t')
                      corname2 = corname2 + '\t\t' + cor[0]
                      corname2b.append(cor[0])
                      corsc2 = corsc2 + '\t' + cor[1]
                      corsc2b.append(cor[1])

            data3 = data1b + data2b
            data3b = random.sample(data3,len(data3))
            lendata1b = len(data1b)
            lendata3 = len(data3)
            data3 = data3b[0:lendata1b]
            data4 = data3b[lendata1b:lendata3]
            mean1 = statistics.mean(data1b)
            median1 = statistics.median(data1b)
            mean2 = statistics.mean(data2b)
            median2 = statistics.median(data2b)
            mean3 = statistics.mean(data3)
            median3 = statistics.median(data3)
            mean4 = statistics.mean(data4)
            median4 = statistics.median(data4)
            mean = '\t' + str(mean1) + '\t' + str(mean2)
            median = '\t' + str(median1) + '\t' + str(median2)
            with open(outfile, mode='a') as f:
                print('\t\tNo. of scores\t\tmean\t\t\tmedian\t',corname1,file=f)
                print('epa1\t\t',len(data1b),'(',len(data12b),')\t',str(mean1),'\t',str(median1),'\t',corsc1,file=f)
                print('epa2\t\t',len(data2b),'(',len(data12b),')\t',str(mean2),'\t',str(median2),'\t',corsc2,file=f)
#                diffmean = abs(mean1 - mean2)
                diffmean = mean1 - mean2
                print('shuffle1\t',len(data3),'\t',str(mean3),'\t',str(median3),file=f)
                print('shuffle2\t',len(data4),'\t',str(mean4),'\t',str(median4),file=f)
                print('',file=f)
#                res = stats.ttest_rel(data1b, data2b)
#                print('result\t',res,file=f)
#                print('pvalue\t',res.pvalue,file=f)
                res = stats.shapiro(data1b)
                print('result data1\t',res,file=f)
                print('pvalue data1\t',res.pvalue,file=f)
                res = stats.shapiro(data2b)
                print('result data2\t',res,file=f)
                print('pvalue data2\t',res.pvalue,file=f)
                res = stats.kstest(data1b, 'norm')
                print('result data1\t',res,file=f)
                print('pvalue data1\t',res.pvalue,file=f)
                res = stats.kstest(data2b, 'norm')
                print('result data2\t',res,file=f)
                print('pvalue data2\t',res.pvalue,file=f)
#                if 'FRRF' in score1:
                if frx1 in score1:
                        total += 1
#                if 'FRRF' in score1 and res.pvalue < pth:
#                        pcou1 += 1
                res = stats.ttest_ind(data1b, data2b, equal_var=True)
                print('result\t',res,file=f)
                print('pvalue\t',res.pvalue,file=f)
#                if 'FRRF' in score1 and res.pvalue < pth:
                if frx1 in score1 and res.pvalue < pth1:
                        pcou2 += 1
                if frx1 in score1 and res.pvalue < pth2:
                        pcou2b += 1
                res = stats.ttest_ind(data1b, data2b, equal_var=False)
                print('result\t',res,file=f)
                print('pvalue\t',res.pvalue,file=f)
#                if 'FRRF' in score1 and res.pvalue < pth:
                if frx1 in score1 and res.pvalue < pth1:
                        pcou3 += 1
                if frx1 in score1 and res.pvalue < pth2:
                        pcou3b += 1
                res = stats.mannwhitneyu(data1b, data2b, alternative='two-sided')
                u1 = res.statistic
                pval = res.pvalue
                e1 = (len(data1b)*len(data2b))/2
                v1 = math.sqrt(len(data1b)*len(data2b)*(len(data1b)+len(data2b)+1)/12)
                z1 = (u1-e1)/v1
                print('result\t',res,'\t','zscore=',z1,file=f)
                print('pvalue\t',pval,'\t',z1,file=f)
#                r1 = math.sqrt(z1**2/(z1**2+len(data1b)+len(data2b)-1))
                u2 = (len(data1b)*len(data2b)) - u1
                if u1 < u2:
                    r1 = 1-(2*u1)/(len(data1b)*len(data2b))
                else:
                    r1 = 1-(2*u2)/(len(data1b)*len(data2b))
                if frx1 in score1 and pval < pth1:
                        pcou5 += 1
                        with open(outfile2, mode='a') as f2:
                            tfids = tfid.split('_',1)
#                            print(tfids[0], '\t', tfids[1], '\t', pval, '\t', z1, '\t', diffmean, '\t', r1, '\t', u1, '\t', len(data1b), '\t', len(data2b), file=f2)
                            print(tfids[0], '\t', pval, '\t', z1, '\t', diffmean, '\t', r1, '\t', u1, '\t', len(data1b), '\t', len(data2b), file=f2)
                if frx1 in score1 and pval < pth2:
                        pcou5b += 1
                        with open(outfileb2, mode='a') as f2:
                            tfids = tfid.split('_',1)
#                            print(tfids[0], '\t', tfids[1], '\t', pval, '\t', z1, '\t', diffmean, '\t', r1, '\t', u1, '\t', len(data1b), '\t', len(data2b), file=f2)
                            print(tfids[0], '\t', pval, '\t', z1, '\t', diffmean, '\t', r1, '\t', u1, '\t', len(data1b), '\t', len(data2b), file=f2)
                if frx1 in score1:
                        tfids = tfid.split('_',1)
#                        data5b = tfids[0] + '\t' + tfids[1] + '\t' + str(diffmean) + '\t' + str(pval) + '\t' + str(z1) + '\t' + str(r1) + '\t' + str(len(data1b)) + '\t' + str(len(data2b))
                        data5b = tfids[0] + '\t' + str(diffmean) + '\t' + str(pval) + '\t' + str(z1) + '\t' + str(r1) + '\t' + str(len(data1b)) + '\t' + str(len(data2b))
                        if pval in data5:
                            data5[float(pval)] = data5[float(pval)] + '@' + data5b
                            data5n += 1
                        else:
                            data5[float(pval)] = data5b
                            data5n += 1
                res = stats.brunnermunzel(data1b, data2b, nan_policy='omit')
                print('result\t',res,file=f)
                print('pvalue\t',res.pvalue,file=f)
                if frx1 in score1 and res.pvalue < pth1:
                        pcou6 += 1
                if frx1 in score1 and res.pvalue < pth2:
                        pcou6b += 1

                print('---shuffle of a pair of scores---',file=f)
#                res = stats.ttest_rel(data3, data4)
#                print('result\t',res,file=f)
#                print('pvalue\t',res.pvalue,file=f)
                res = stats.ttest_ind(data3, data4, equal_var=True)
                print('result\t',res,file=f)
                print('pvalue\t',res.pvalue,file=f)
                res = stats.ttest_ind(data3, data4, equal_var=False)
                print('result\t',res,file=f)
                print('pvalue\t',res.pvalue,file=f)
#                res = stats.wilcoxon(data3, data4)
#                print('result\t',res,file=f)
#                print('pvalue\t',res.pvalue,file=f)
                res = stats.mannwhitneyu(data3, data4, alternative='two-sided')
                print('result\t',res,file=f)
                print('pvalue\t',res.pvalue,file=f)
                res = stats.brunnermunzel(data3, data4, nan_policy='omit')
                print('result\t',res,file=f)
                print('pvalue\t',res.pvalue,file=f)
                print('--------------------------------',file=f)

            data1b.clear
            data2b.clear
            corname1 = ''
            corname2 = ''
            corname1b = ''
            corname2b = ''
            corsc1 = ''
            corsc2 = ''
            corsc1b = ''
            corsc2b = ''

if total != 0:
    data5c = 0
    data5q = 0
    data5s = []
    pcou5c = 0
    data5_sorted = sorted(data5.items(), key=lambda x:x[0])
    print('TF_name\tMean_difference\tpvalue\tzscore\trvalue\t#_score_1\t#_score_2')
    for key1, value1 in data5_sorted:
        data5s = value1.split('@')
        for value2 in data5s:
            data5c += 1
            data5q = key1 * data5n / data5c
            if data5q < pth0:
#                print(value2, '\t' , data5c, '\t', data5q)
                print(value2)
                with open(outfile3, mode='a') as f2:
                    print(value2, '\t' , data5c, '\t', data5q, file=f2)
                pcou5c += 1
        
    shutil.copy2(outfile, outfileb)
#    ratio1 = pcou1 / total
    ratio2 = pcou2 / total
    ratio3 = pcou3 / total
#    ratio4 = pcou4 / total
    ratio5 = pcou5 / total
    ratio6 = pcou6 / total
    ratio2b = pcou2b / total
    ratio3b = pcou3b / total
    ratio5b = pcou5b / total
    ratio6b = pcou6b / total
    ratio5c = pcou5c / total

#    print('outfile4',outfile4)
    print('')
    print('pvalue_threshold_1\t',pthb)
    print('score_threshold_1\t',scthb)
    print('outfile',outfile)
    print('outfile2',outfile2)
    print('\tTotal\tSignificant\tRatio')
#    print(total,'\t',pcou1,'\t',ratio1)
    print('Test1\t',total,'\t',pcou2,'\t',ratio2)
    print('Test2\t',total,'\t',pcou3,'\t',ratio3)
#    print(total,'\t',pcou4,'\t',ratio4)
    print('Test3\t',total,'\t',pcou5,'\t',ratio5)
    print('Test4\t',total,'\t',pcou6,'\t',ratio6)
    print('')

    print('pvalue_threshold_2\t',pth2b)
    print('score_threshold_2\t',scth2b)
    print('outfileb',outfileb)
    print('outfileb2',outfileb2)
    print('\tTotal\tSignificant\tRatio')
#    print(total,'\t',pcou1b,'\t',ratio1b)
    print('Test1\t',total,'\t',pcou2b,'\t',ratio2b)
    print('Test2\t',total,'\t',pcou3b,'\t',ratio3b)
#    print(total,'\t',pcou4b,'\t',ratio4b)
    print('Test3\t',total,'\t',pcou5b,'\t',ratio5b)
    print('Test4\t',total,'\t',pcou6b,'\t',ratio6b)
    print('')

    print('pvalue_threshold_3\t',pth0b)
    print('score_threshold_3\t',scth0b)
    print('outfile3',outfile3)
    print('Total\tSignificant\tRatio')
    print(total,'\t',pcou5c,'\t',ratio5c)
    print('')

    with open(outfile4, mode='a') as f:
        print('pth1',pthb,file=f)
        print('scth1',scthb,file=f)
        print('outfile',outfile,file=f)
        print('outfile2',outfile2,file=f)
        print('Total\tSignificant\tRatio',file=f)
        print(total,'\t',pcou2,'\t',ratio2,file=f)
        print(total,'\t',pcou3,'\t',ratio3,file=f)
        print(total,'\t',pcou5,'\t',ratio5,file=f)
        print(total,'\t',pcou6,'\t',ratio6,file=f)
        print('',file=f)
        print('pth2',pth2b,file=f)
        print('scth2',scth2b,file=f)
        print('outfile',outfile,file=f)
        print('outfileb2',outfileb2,file=f)
        print('Total\tSignificant\tRatio',file=f)
        print(total,'\t',pcou2b,'\t',ratio2b,file=f)
        print(total,'\t',pcou3b,'\t',ratio3b,file=f)
        print(total,'\t',pcou5b,'\t',ratio5b,file=f)
        print(total,'\t',pcou6b,'\t',ratio6b,file=f)
        print('',file=f)
        print('pth0',pth0b,file=f)
        print('scth0',scth0b,file=f)
        print('outfile3',outfile3,file=f)
#        print('Total\tSignificant\tRatio',file=f)
#        print(total,'\t',pcou5c,'\t',ratio5c,file=f)

