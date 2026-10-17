from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[2]
r=root/'results';o=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.bbox':'tight'})
def rows(name):return list(csv.DictReader((r/name).open()))
cs=['#126782','#b84b28','#7161a8','#52854c']
a=rows('one_erasure.csv');f,ax=plt.subplots(figsize=(4.6,2.9))
for k,l,c,st in [('complete','Complete grid',cs[0],'-'),('one_missing','One missing: optimal',cs[1],'-.'),('asymptotic','Asymptotic',cs[2],'--'),('polynomial_lower_bound','Polynomial lower bound',cs[3],':')]:
 ax.semilogy([float(x['N']) for x in a],[float(x[k]) for x in a],st,label=l,color=c)
ax.set(xlabel='Grid size N',ylabel='Sharp unit-ball error');ax.legend(fontsize=8);f.tight_layout();f.savefig(o/'one_erasure.pdf');plt.close(f)
a=rows('geometry.csv');f,axs=plt.subplots(2,1,figsize=(4.8,5.5))
for p,c,st in zip(['block','dispersed','random'],cs,['o-','s--','^-.']):
 xs=[1,2,4,8,12,16];ys=[np.median([float(x['error_norm']) for x in a if x['pattern']==p and int(x['erased'])==e and x['method']=='recovery']) for e in xs]
 axs[0].semilogy(xs,ys,st,label=p,color=c,markersize=4)
axs[0].set(xlabel='Missing evaluations',ylabel='Sharp unit-ball error',title='(a) Geometry: N = 64, q = 0.8');axs[0].legend(fontsize=8)
names=['recovery','equal','nonuniform_trapezoid','fourier8','fourier_max']
for p,st,c in [('block','o-',cs[0]),('dispersed','s--',cs[1])]:
 vals=[float(next(x['error_norm'] for x in a if x['pattern']==p and x['erased']=='8' and x['method']==name)) for name in names]
 axs[1].semilogy(range(5),vals,st,label=p,color=c,markersize=4)
axs[1].set_xticks(range(5),['Recovery','Equal','Trap.','Fourier 8','Fourier max'],rotation=15)
axs[1].set(ylabel='Sharp unit-ball error',title='(b) Same eight erasures');axs[1].legend(fontsize=8);f.tight_layout(h_pad=2);f.savefig(o/'geometry.pdf');plt.close(f)
a=json.loads((r/'retries.json').read_text());f,ax=plt.subplots(figsize=(4.6,3.0));xs=range(13)
ax.fill_between(xs,a['random_q10'],a['random_q90'],alpha=.18,color=cs[1],label='Random: 10–90%')
ax.semilogy(xs,a['random_median'],'--',color=cs[1],label='Random: median')
ax.semilogy(xs,a['greedy_errors'],'o-',color=cs[0],markersize=3,label='Exact one-step gain')
ax.semilogy(range(5),a['optimal_first_five'],'x',color='black',label='Exhaustive optimum')
ax.set(xlabel='Restored evaluations',ylabel='Sharp unit-ball error');ax.legend(fontsize=8);f.tight_layout();f.savefig(o/'retries.pdf');plt.close(f)
a=rows('timing.csv');f,ax=plt.subplots(figsize=(4.4,2.9))
for name,st,c in [('fft','o-',cs[0]),('dense','s--',cs[1])]:
 t=[x for x in a if x['method']==name];ax.loglog([int(x['N']) for x in t],[1000*float(x['median_seconds']) for x in t],st,label=name,color=c,markersize=4)
ax.set(xlabel='Grid size N',ylabel='Median setup time (ms)');ax.legend(fontsize=8);f.tight_layout();f.savefig(o/'timing.pdf');plt.close(f)
