
# coding: utf-8

# In[1]:


from util import *
from objects import *
from collections import defaultdict, Counter

import seaborn as sb
import matplotlib.pylab as plt



# In[46]:


drugs,proteins = parse(
    '/home/stanislaw/Downloads/git_drugbank/download/full_database.xml')


def get_multitarget_drugs():
    result = []
    dlist = list(drugs.values())
    for d in dlist:
        if len(d.edges) > 1:
            result.append(d)
    return result

multi = get_multitarget_drugs()


s = []
for a in multi:
    for b in a.edges:
        try:
            s+= list(b[1].bio_processes)
        except AttributeError:
            pass


cnt = Counter()
for el in s:
    cnt[el]+=1


t = sorted(list(cnt.items()),key = lambda f: f[1])

t2 = t[-20:]
a = [i[1] for i in t2]
b = [i[0] for i in t2]

plt.figure(figsize=(10,8))
sb.barplot(x=a,y=b,color='b')
plt.title('Processes influenced by multitarget drugs',fontsize=20)

plt.show()


cnt3 = Counter()
for dru in multi:
    if dru.patent_date is not None:
        cnt3[int(dru.patent_date.split('-')[0])]+=1

cnt4 = Counter()
for dru in drugs.values():
    if dru.patent_date is not None and dru not in m2:
        cnt4[int(dru.patent_date.split('-')[0])]+=1


for key in cnt4:
    if key not in cnt3:
        cnt3[key]=0
        
for key in cnt3:
    if key not in cnt4:
        cnt4[key]=0



t = sorted(list(cnt3.items()),key=lambda u: u[0])[10:]
t2 = sorted(list(cnt4.items()),key=lambda u: u[0])[10:]
x = [u[0] for u in t]
y1 = [u[1] for u in t]
y2 = [u[1] for u in t2]
su = [y1[i]+y2[i] for i in range(32)]
y1n= [y1[i]/su[i] for i in range(32)]
y2n= [y2[i]/su[i] for i in range(32)]
plt.figure(figsize=(10,8))
ps = [plt.bar(x,y1n),plt.bar(x,y2n,bottom=y1n)]
plt.title('Ratio of multitarget drugs vs. other drugs by year',fontsize=20)
plt.xlim((1982,2021))
plt.ylim((0,1.2))
plt.legend(ps,['multitarget','other'])
plt.show()


