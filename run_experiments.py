"""Reproduce all tables, figures and recorded numerical evidence."""
import csv,itertools,json,platform,time,sys
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import scipy, scipy.linalg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.contour_recovery import recover,dense_rule,error_norm,kernel,greedy_retries

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results';FIGS=ROOT/'paper'/'figures'
RESULTS.mkdir(exist_ok=True);FIGS.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'savefig.bbox':'tight','pdf.fonttype':42})
COLORS=['#126782','#d05d36','#7161a8','#52854c']

def savecsv(name,rows):
    with (RESULTS/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def frac_solve(a,b):
    n=len(b);a=[list(x)+[y] for x,y in zip(a,b)]
    for j in range(n):
        p=next(i for i in range(j,n) if a[i][j]);a[j],a[p]=a[p],a[j]
        t=a[j][j];a[j]=[v/t for v in a[j]]
        for i in range(n):
            if i!=j:
                t=a[i][j];a[i]=[x-t*y for x,y in zip(a[i],a[j])]
    return [x[-1] for x in a]

def exact_counterexample():
    q=F(1,4);cos=[F(1),F(0),F(-1),F(0)]
    k=[(1-q*q)/(1-2*q*c+q*q) for c in cos]
    def v(a):
        w=frac_solve([[k[(i-j)%4] for j in a] for i in a],[F(1)]*len(a))
        return 1-sum(w)
    a=[0,2];b=[0,2,3];j=1
    d0=v(a)-v(a+[j]);d1=v(b)-v(b+[j])
    assert d1>d0
    out={'q':'1/4','N':4,'A':a,'B':b,'j':j,'gain_A':str(d0),'gain_B':str(d1),
         'difference':str(d1-d0),'gain_A_float':float(d0),'gain_B_float':float(d1)}
    (RESULTS/'exact_counterexample.json').write_text(json.dumps(out,indent=2))
    return out

def baselines(n,missing):
    a=np.array([j for j in range(n) if j not in set(missing)]);theta=2*np.pi*a/n;m=len(a)
    uniform=np.zeros(n);uniform[a]=1/m
    gaps=(np.roll(theta,-1)-theta)%(2*np.pi)
    periodic=np.zeros(n);periodic[a]=(gaps+np.roll(gaps,1))/(4*np.pi)
    out={'equal':uniform,'nonuniform_trapezoid':periodic}
    for label,d in [('fourier8',min(8,(m-1)//2)),('fourier_max',(m-1)//2)]:
        modes=np.arange(-d,d+1)
        mat=np.exp(1j*theta[:,None]*modes)
        w=np.zeros(n);w[a]=np.linalg.pinv(mat,rcond=1e-12)[d,:].real
        out[label]=w
    return out

def geometry():
    rows=[];n=64;q=.8;rng=np.random.default_rng(101)
    for e in [1,2,4,8,12,16]:
        for pattern in ['block','dispersed','random']:
            masks=[list(range(e))] if pattern=='block' else [list((np.arange(e)*n//e).astype(int))]
            if pattern=='random':masks=[sorted(rng.choice(n,e,replace=False).tolist()) for _ in range(100)]
            for rep,missing in enumerate(masks):
                r=recover(n,q,missing)
                ws={'recovery':r.weights,**baselines(n,missing)}
                for name,w in ws.items():
                    rows.append(dict(N=n,q=q,erased=e,pattern=pattern,rep=rep,method=name,
                        error_norm=error_norm(w,q),weight_l1=float(abs(w).sum())))
    savecsv('geometry.csv',rows)
    fig,ax=plt.subplots(1,2,figsize=(10.5,3.6))
    for p,c in zip(['block','dispersed','random'],COLORS):
        xs=[1,2,4,8,12,16];ys=[]
        for e in xs:ys.append(np.median([r['error_norm'] for r in rows if r['pattern']==p and r['erased']==e and r['method']=='recovery']))
        ax[0].semilogy(xs,ys,'o-',label=p,color=c)
    ax[0].set(xlabel='Missing evaluations',ylabel='Sharp unit-ball error',title='(a) Geometry at N = 64, q = 0.8');ax[0].legend()
    names=['recovery','equal','nonuniform_trapezoid','fourier8','fourier_max']
    for pattern,style in [('block','-'),('dispersed','--')]:
        vals=[next(r['error_norm'] for r in rows if r['pattern']==pattern and r['erased']==8 and r['method']==name) for name in names]
        ax[1].semilogy(range(5),vals,'o'+style,label=pattern)
    ax[1].set_xticks(range(5),['Recovery','Equal','Trap.','Fourier 8','Fourier max'],rotation=22)
    ax[1].set(ylabel='Sharp unit-ball error',title='(b) Same eight erasures');ax[1].legend()
    fig.tight_layout();fig.savefig(FIGS/'geometry.pdf');plt.close(fig)
    return rows

def rational_benchmark():
    rng=np.random.default_rng(421);n=64;q=.8;z=np.exp(2j*np.pi*np.arange(n)/n);rows=[]
    for case in range(40):
        pi=rng.uniform(.25,.72,3)*np.exp(2j*np.pi*rng.random(3))
        po=rng.uniform(1.4,2.2,3)*np.exp(2j*np.pi*rng.random(3))
        ri=rng.normal(size=3)+1j*rng.normal(size=3)
        ro=rng.normal(size=3)+1j*rng.normal(size=3)
        c0=rng.normal()+1j*rng.normal()
        a0=c0-np.sum(ro/po)
        ii=pi[:,None]*pi.conj()[None,:];oo=po[:,None]*po.conj()[None,:]
        norm2=abs(a0)**2+np.sum(ri[:,None]*ri.conj()[None,:]/(q-ii)).real
        norm2+=np.sum(ro[:,None]*ro.conj()[None,:]/(oo*(q*oo-1))).real
        norm=np.sqrt(norm2)
        y=(c0+np.sum(ri/(z[:,None]-pi),axis=1)+np.sum(ro/(z[:,None]-po),axis=1))/norm
        exact=a0/norm
        for pattern in ['block','random','dispersed']:
            if pattern=='block':e=list((np.arange(8)+case)%n)
            elif pattern=='dispersed':e=list((np.arange(8)*8+case)%n)
            else:e=sorted(rng.choice(n,8,replace=False).tolist())
            r=recover(n,q,e);ws={'recovery':r.weights,**baselines(n,e)}
            for name,w in ws.items():
                err=abs(w@y-exact);bound=error_norm(w,q)
                assert err<=bound+1e-12
                rows.append(dict(case=case,pattern=pattern,method=name,error=float(err),bound=bound,ratio=float(err/bound)))
    savecsv('rational.csv',rows)
    return rows

def matrix_benchmark():
    rng=np.random.default_rng(872);d=48;n=64;q=.8
    u,_=np.linalg.qr(rng.normal(size=(d,d)));v,_=np.linalg.qr(rng.normal(size=(d,d)))
    V=u@np.diag(np.geomspace(1,10,d))@v.T
    eigs=np.r_[.65*np.exp(2j*np.pi*np.arange(16)/16),1.6*np.exp(2j*np.pi*np.arange(32)/32)]
    A=V@np.diag(eigs)@np.linalg.inv(V);b=rng.normal(size=d);c=np.linalg.solve(V,b)
    reference=V[:,:16]@c[:16]
    M=np.linalg.norm(V,2)*np.sqrt(np.sum(abs(c[:16])**2/(1-abs(eigs[:16])**2/q))+
                np.sum(abs(c[16:])**2/(q*abs(eigs[16:])**2-1)))
    z=np.exp(2j*np.pi*np.arange(n)/n)
    y=np.array([zz*np.linalg.solve(zz*np.eye(d)-A,b) for zz in z])
    rows=[]
    for pattern,e in [('block',list(range(8))),('dispersed',list(range(0,64,8))),
                      ('random',sorted(rng.choice(n,8,replace=False).tolist()))]:
        r=recover(n,q,e)
        for name,w in {'recovery':r.weights,**baselines(n,e)}.items():
            error=float(np.linalg.norm(w@y-reference));bound=M*error_norm(w,q)
            assert error<=bound+1e-10
            rows.append(dict(pattern=pattern,method=name,error=error,relative_error=error/np.linalg.norm(reference),bound=bound,norm_bound=float(M),dimension=d))
    savecsv('projector.csv',rows)
    return rows

def timings():
    rng=np.random.default_rng(983);rows=[]
    for n in [64,128,256,512,1024]:
        q=float(np.exp(-12/n));e=sorted(rng.choice(n,8,replace=False).tolist())
        records={}
        for name,fun in [('fft',lambda:recover(n,q,e).weights),('dense',lambda:dense_rule(n,q,e))]:
            fun();times=[]
            for _ in range(7):
                t=time.perf_counter();w=fun();times.append(time.perf_counter()-t)
            records[name]=(w,times)
        delta=float(np.max(abs(records['fft'][0]-records['dense'][0])))
        for name,(_,ts) in records.items():
            rows.append(dict(N=n,q=q,erased=8,method=name,median_seconds=float(np.median(ts)),min_seconds=min(ts),max_seconds=max(ts),weight_difference=delta))
    savecsv('timing.csv',rows)
    fig,ax=plt.subplots(figsize=(5.7,3.5))
    for name,c in zip(['fft','dense'],COLORS):
        rs=[r for r in rows if r['method']==name]
        ax.loglog([r['N'] for r in rs],[1000*r['median_seconds'] for r in rs],'o-',label=name,color=c)
    ax.set(xlabel='Grid size N',ylabel='Median setup time (ms)',title='Eight erasures; q = exp(−12/N)');ax.legend();fig.tight_layout();fig.savefig(FIGS/'timing.pdf');plt.close(fig)
    return rows

def retries():
    n=64;q=.8;e=[0,1,2,3,4,5,6,7,16,24,40,56]
    order,errs=greedy_retries(n,q,e,len(e));rng=np.random.default_rng(784)
    random_errors=[]
    for _ in range(100):
        p=list(rng.permutation(e));random_errors.append([recover(n,q,[j for j in e if j not in p[:k]]).error for k in range(len(e)+1)])
    random_errors=np.array(random_errors)
    opts=[]
    for k in range(5):
        opts.append(min(recover(n,q,set(e)-set(s)).error for s in itertools.combinations(e,k)))
    out={'missing':e,'greedy_order':order,'greedy_errors':errs,'optimal_first_five':opts,
         'random_median':np.median(random_errors,axis=0).tolist(),
         'random_q10':np.quantile(random_errors,.1,axis=0).tolist(),'random_q90':np.quantile(random_errors,.9,axis=0).tolist()}
    (RESULTS/'retries.json').write_text(json.dumps(out,indent=2))
    fig,ax=plt.subplots(figsize=(6.0,3.7));ks=np.arange(len(e)+1)
    ax.fill_between(ks,out['random_q10'],out['random_q90'],alpha=.15,color=COLORS[1],label='Random: 10–90%')
    ax.semilogy(ks,out['random_median'],'--',color=COLORS[1],label='Random: median')
    ax.semilogy(ks,errs,'o-',color=COLORS[0],label='Exact marginal gain')
    ax.semilogy(range(5),opts,'x',markersize=8,color='black',label='Exhaustive optimum')
    ax.set(xlabel='Restored evaluations',ylabel='Sharp unit-ball error',title='Repairing a mixed outage');ax.legend();fig.tight_layout();fig.savefig(FIGS/'retries.pdf');plt.close(fig)
    return out

def one_erasure():
    q=.8;rows=[]
    s_inf=.5+2*sum(1/(q**k+q**(-k)) for k in range(1,1000))
    for n in range(8,193,8):
        m=n//2;t=q**n
        s=sum(1/(q**k+q**(-k)) for k in range(-m,m))
        closed=2*t/(1+t)+(1-t)*q**m/((1+t)**2*s)
        r=recover(n,q,[0])
        assert abs(r.error**2/closed-1)<1e-7
        poly_norm=1+2*sum(q**(-k) for k in range(1,m))+q**(-m)
        rows.append(dict(N=n,q=q,complete=np.sqrt(2*t/(1+t)),one_missing=r.error,
                         closed_form=np.sqrt(closed),asymptotic=q**(n/4)/np.sqrt(s_inf),
                         polynomial_lower_bound=1/np.sqrt(poly_norm)))
    savecsv('one_erasure.csv',rows)
    fig,ax=plt.subplots(figsize=(6.2,3.8))
    for key,label,c,sty in [('complete','Complete grid',COLORS[0],'-'),('one_missing','One missing: optimal recovery',COLORS[1],'-'),
                           ('asymptotic','Exact asymptotic',COLORS[2],'--'),('polynomial_lower_bound','Vanishing polynomial lower bound',COLORS[3],':')]:
        ax.semilogy([r['N'] for r in rows],[r[key] for r in rows],sty,label=label,color=c)
    ax.set(xlabel='Grid size N',ylabel='Worst-case absolute error (unit ball)',title='One erasure halves the asymptotic exponent (q = 0.8)')
    ax.legend(fontsize=8.5);fig.tight_layout();fig.savefig(FIGS/'one_erasure.pdf');plt.close(fig)
    return rows

def main():
    meta={'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,
          'matplotlib':matplotlib.__version__,'platform':platform.platform(),
          'seeds':[101,421,872,983,784], 'timing_repeats':7}
    (RESULTS/'environment.json').write_text(json.dumps(meta,indent=2))
    out={'exact':exact_counterexample()}
    for name,fun in [('one_erasure',one_erasure),('geometry',geometry),('rational',rational_benchmark),('projector',matrix_benchmark),('timing',timings),('retries',retries)]:
        print('Running',name,flush=True);data=fun();out[name]=len(data);print('Done',name,flush=True)
    (RESULTS/'run_summary.json').write_text(json.dumps(out,indent=2))

if __name__=='__main__':main()
