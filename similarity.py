
# coding: utf-8

# In[1]:


from util import *
from objects import *
from collections import defaultdict

import numpy as np
import seaborn as sb
import matplotlib.pylab as plt

drugs,proteins = parse(
    '/home/stanislaw/Downloads/git_drugbank/download/full_database.xml')

def write_smiles():
    f = open('drugs.smi','w')
    for d in drugs.values():
        if d.smiles is not None:
            f.write('\t'.join([d.smiles,d.dbid])+'\n')
    f.close()


f = open('similarity-slim.tsv','r')
d = defaultdict(list)
sims = []
f.readline()
for line in f:
    l = line.split()
    sims.append(float(l[2]))
    d2[l[0]].append((l[1],float(l[2])))

f.close()

sim_arr = np.array(sims)
plt.figure(figsize=(10,8))
sb.distplot(sim_arr)
plt.title('Drug similarity distribution',fontsize=20)
plt.show()


common_targets = {}
for c in d:
    li = []
    for b,a in d[c]:
        try:
            li.append((b,len(set.intersection(set(drugs[c].edges),set(drugs[b].edges)))))
        except KeyError:
            pass
    common_targets[c] = li


v = {}
for c in common_targets:
    res = []
    for did,count in common_targets[c]:
        try:
            res.append((count/len(drugs[c].edges),did,count/len(drugs[did].edges)))
        except ZeroDivisionError:
            res.append((0,did,0))
        
    v[c]=res


a=0
b=0
for c in v.values():
    a+=len(c)
    for el in c:
        b+=el[0]+el[2]

print('Average:',b/a)

