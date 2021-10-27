'''
Isomer_finder
'''

import sys

results = []
results2 = []


'''
datas = ['/home/joe/Projects/clusters/Au6Co1/gms_s/_data.txt',
         '/home/joe/Projects/clusters/Au6Co1/gms_t/_data.txt',
         '/home/joe/Projects/clusters/Au6Co1/gms_5/_data.txt',
         '/home/joe/Projects/clusters/Au6Co1/gms_7/_data.txt',]

datas = ['/home/joe/Projects/clusters/Au6Ni1/gms_2/_data.txt',
         '/home/joe/Projects/clusters/Au6Ni1/gms_4/_data.txt',
         '/home/joe/Projects/clusters/Au6Ni1/gms_6/_data.txt',]
'''
datas = ['/media/joe/Seagate Portable Drive/Research/clusters/CuAu6/gms_a1/_data.txt',
         '/media/joe/Seagate Portable Drive/Research/clusters/CuAu6/gms_a3/_data.txt',
         '/media/joe/Seagate Portable Drive/Research/clusters/CuAu6/gms_a5/_data.txt']

for d in datas:
    with open(d) as data:
        for line in data:
            if line.find('DONE') < 0:
                results.append(line.split(';'))

results.sort(key=lambda e: e[1])
for r in results:
    if float(r[1].split('=')[1]) == 0:
        results.remove(r)
    else:
        results2.append([r[0], float(r[1].split('=')[1])])


l=len(results2)
print(results2[l-1][0] + '    0')
E=results2[l-1][1]
for k in range(1, l):
    d = results2[l-k-1][1]-results2[l-k][1]
    if d > 0.0001:
        print(format(results2[l-k-1][0] + '   ' + str(27.2*(results2[l-k-1][1] - E))))