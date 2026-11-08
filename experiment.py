"""Reproducible Newman1 PCPSP pilot. MineLib: continuous destination fractions.
Run: python experiment.py --stage prepare|baseline|search|summarize
No reference solution is used to initialize the proposed method.
"""
from pathlib import Path
import argparse, json, time, heapq, platform, os, hashlib
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import highspy as hs

ROOT=Path(__file__).resolve().parent; OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
N=1060; T=6; M=N*T; NV=2*M
raw=(ROOT/'upload/newman1(4).pcpsp').read_text().splitlines()
v=np.array([list(map(float,l.split()[1:])) for l in raw[22:1082]])
blocks=pd.read_csv(ROOT/'upload/newman1(3).blocks',sep=r'\s+',header=None)
ton=blocks[6].to_numpy(); z=blocks[3].to_numpy()
pred=[[] for _ in range(N)]; succ=[[] for _ in range(N)]
for l in (ROOT/'upload/newman1(4).prec').read_text().splitlines():
    b,n,*pp=map(int,l.split()); assert n==len(pp);pred[b]=pp
    for p in pp: succ[p].append(b)
edges=[(b,p) for b in range(N) for p in pred[b]]
disc=1.08**(-np.arange(T))
c=np.r_[(v[:,0,None]*disc).ravel(),((v[:,1]-v[:,0])[:,None]*disc).ravel()]
rr=[];cc=[];aa=[];lo=[];up=[]
def row(ind,val,l=-np.inf,u=0):
    k=len(lo);rr.extend([k]*len(ind));cc.extend(ind);aa.extend(val);lo.append(l);up.append(u)
for b in range(N): row(list(range(b*T,(b+1)*T)),[1]*T,u=1)
for b,p in edges:
    for t in range(T): row([b*T+j for j in range(t+1)]+[p*T+j for j in range(t+1)],[1]*(t+1)+[-1]*(t+1))
for j in range(M): row([M+j,j],[1,-1])
capacity_start=len(lo)
for r,cap in [(0,2e6),(1,1.1e6)]:
    for t in range(T): row([r*M+b*T+t for b in range(N)],list(ton),u=cap)
A=coo_matrix((aa,(rr,cc)),shape=(len(lo),NV)).tocsr();lo=np.array(lo);up=np.array(up)
def model(integer=False,limit=60,seed=42):
    lp=hs.HighsLp();lp.num_col_=NV;lp.num_row_=A.shape[0];lp.col_cost_=c
    lp.col_lower_=np.zeros(NV);lp.col_upper_=np.ones(NV);lp.row_lower_=lo;lp.row_upper_=up
    lp.sense_=hs.ObjSense.kMaximize
    lp.a_matrix_.format_=hs.MatrixFormat.kRowwise;lp.a_matrix_.start_=A.indptr;lp.a_matrix_.index_=A.indices;lp.a_matrix_.value_=A.data
    h=hs.Highs();h.setOptionValue('output_flag',False);h.setOptionValue('threads',1)
    h.setOptionValue('random_seed',int(seed));h.setOptionValue('time_limit',float(limit));h.setOptionValue('mip_rel_gap',1e-6)
    h.passModel(lp)
    if integer: h.changeColsIntegrality(M,np.arange(M,dtype=np.int32),np.ones(M,dtype=np.uint8))
    return h
def check(s):
    x=s[:M].reshape(N,T);y=s[M:].reshape(N,T)
    mine=(ton[:,None]*x).sum(0);proc=(ton[:,None]*y).sum(0)
    precedence=max([max(np.cumsum(x[b])-np.cumsum(x[p])) for b,p in edges])
    ax=A@s
    return dict(objective=float(c@s),max_row_violation=float(max(np.max(ax-up),np.max(lo-ax),0)),
       integrality_error=float(np.max(abs(x-np.rint(x)))),bound_violation=float(max(0,-s.min(),s.max()-1)),
       precedence_violation=float(max(0,precedence)),mined_blocks=int(np.rint(x).sum()),
       mine_tonnes=mine.tolist(),processed_tonnes=proc.tolist(),
       fractional_routing_blocks=int(((y.sum(1)>1e-7)&(y.sum(1)<1-1e-7)).sum()))
def valid(s):
    q=check(s);return q['max_row_violation']<.05 and q['integrality_error']<1e-5 and q['bound_violation']<1e-5
def periods(s):
    x=s[:M].reshape(N,T);return np.where(x.sum(1)>.5,x.argmax(1),T)
def route(ts):
    x=np.zeros((N,T));y=np.zeros((N,T))
    for b,t in enumerate(ts):
        if t<T:x[b,t]=1
    for t in range(T):
        cap=1.1e6
        for b in sorted(np.where(ts==t)[0],key=lambda b:-(v[b,1]-v[b,0])/ton[b]):
            if v[b,1]<=v[b,0]:continue
            a=min(1,cap/ton[b]);y[b,t]=a;cap-=a*ton[b]
            if cap<1e-8:break
    return np.r_[x.ravel(),y.ravel()]
def prepare():
    t=time.perf_counter();h=model();h.run();elapsed=time.perf_counter()-t
    assert h.getModelStatus()==hs.HighsModelStatus.kOptimal
    ss=h.getSolution();s=np.array(ss.col_value);dual=np.array(ss.row_dual);ub=c@s
    x=s[:M].reshape(N,T);tau=x@np.arange(T)+T*(1-x.sum(1))
    # Positive mass set plus all ancestors. Dummy period T encodes not mined.
    selected=set(np.where(x.sum(1)>1e-7)[0]);stack=list(selected)
    while stack:
        b=stack.pop()
        for p in pred[b]:
            if p not in selected:selected.add(p);stack.append(p)
    # Topological priority sweep: LP expected time, then LP marginal economic value.
    ts=np.full(N,T,dtype=int);done=set()
    for t in range(T):
        cap=2e6;queue=[]
        def push(b):heapq.heappush(queue,(tau[b],-v[b].max()/ton[b],b))
        for b in selected-done:
            if all(p in done for p in pred[b]):push(b)
        while queue:
            _,_,b=heapq.heappop(queue)
            if b in done or ton[b]>cap+1e-7:continue
            ts[b]=t;done.add(b);cap-=ton[b]
            for q in succ[b]:
                if q in selected and q not in done and all(p in done for p in pred[q]):push(q)
    initial=route(ts);assert valid(initial)
    ref=np.zeros(NV)
    for l in (ROOT/'upload/newman1_pcpsp_gmunoz120723.sol').read_text().splitlines():
        b,d,t,a=l.split();b=int(b);d=int(d);t=int(t);a=float(a)
        ref[b*T+t]+=a
        if d==1:ref[M+b*T+t]+=a
    assert valid(ref)
    supplied=np.loadtxt(ROOT/'upload/newman1.LPcpit');file_lp=np.zeros(NV)
    for b,d,t,a in supplied:
        b=int(b);d=int(d);t=int(t);file_lp[b*T+t]+=a
        if d==1:file_lp[M+b*T+t]+=a
    np.savez(OUT/'prepared.npz',lp=s,dual=dual,tau=tau,initial=initial,reference=ref,ub=ub)
    report=dict(lp_objective=ub,lp_seconds=elapsed,lp_status=str(h.getModelStatus()),variables=NV,binary_variables=M,
      rows=A.shape[0],precedence_arcs=len(edges),initial=check(initial),reference=check(ref),
      supplied_cp_expected_objective=float(sum(v[int(b)].max()*a/1.08**t for b,d,t,a in supplied)),
      supplied_lp_note='CPIT file: destination 0 denotes its sole fixed-destination variable, not the PCPSP waste destination. Not used as PCPSP LP.',
      total_tonnes=float(ton.sum()),selected_lp_mass=float(x.sum()),unmined_lp_blocks=int((x.sum(1)<1e-7).sum()),
      capacity_duals=dual[capacity_start:].reshape(2,T).tolist(),
      versions=dict(highs=h.version(),python=platform.python_version(),platform=platform.platform(),threads=1,cpu_count=os.cpu_count()),
      sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'upload').iterdir()})
    (OUT/'prepare.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
def baseline(seconds,start_kind='initial'):
    d=np.load(OUT/'prepared.npz');initial=d[start_kind];h=model(True,seconds)
    stem='baseline' if start_kind=='initial' else 'baseline_reference'
    h.setOptionValue('log_file',str(OUT/f'{stem}_solver.log'));h.setOptionValue('log_to_console',False);h.setOptionValue('output_flag',True)
    h.setSolution(NV,np.arange(NV,dtype=np.int32),initial)
    t=time.perf_counter();h.run();elapsed=time.perf_counter()-t;s=np.array(h.getSolution().col_value)
    if not valid(s):raise RuntimeError(check(s))
    inf=h.getInfo();q=dict(method='full_mip' if start_kind=='initial' else 'full_mip_reference',seconds=elapsed,status=str(h.getModelStatus()),dual_bound=inf.mip_dual_bound,nodes=inf.mip_node_count,**check(s))
    (OUT/f'{stem}.json').write_text(json.dumps(q,indent=2));np.save(OUT/f'{stem}.npy',s);print(json.dumps(q),flush=True)
def ancestor_close(seeds):
    closed=set(map(int,seeds));stack=list(closed)
    while stack:
        b=stack.pop()
        for p in pred[b]:
            if p not in closed:closed.add(p);stack.append(p)
    return closed
def normalized(a):
    return a/max(float(np.max(a)),1e-12)
def scores(inc,lp,tau,dual):
    x=lp[:M].reshape(N,T);dummy=1-x.sum(1)
    fraction=1-(x*x).sum(1)-dummy**2
    disp=abs(periods(inc)-tau)/T
    lam=abs(dual[capacity_start:].reshape(2,T))
    pressure=ton*(x@lam[0]+lp[M:].reshape(N,T)@lam[1])
    return (normalized(fraction)+normalized(disp)+normalized(pressure))/3
def choose_neighborhood(kind,size,score,tau,dual,rng,random_mode=False):
    # Ancestor closure can exceed requested seed count; never truncate a closure.
    jitter=rng.uniform(0,.04,N);ss=score+jitter
    if random_mode:rank=rng.permutation(N)
    elif kind==0:
        seed=int(np.argmax(ss));rank=np.argsort(abs(tau-tau[seed])-.25*ss)
    elif kind==1:
        lam=abs(dual[capacity_start:].reshape(2,T));bottleneck=int(np.argmax(lam.sum(0)))
        rank=np.argsort(abs(tau-bottleneck)-ss)
    elif kind==2:
        seed=int(np.argmax(ss));rank=np.argsort(abs(z-z[seed])-.5*ss)
    else:
        seed=int(np.argmax(ss));visited=set();queue=[(-ss[seed],seed)];rank=[]
        while queue:
            _,b=heapq.heappop(queue)
            if b in visited:continue
            visited.add(b);rank.append(b)
            for p in pred[b]+succ[b]:
                if p not in visited:heapq.heappush(queue,(-ss[p],p))
        rank+=list(np.argsort(-ss))
    free=set()
    for b in rank:
        if int(b) in free:continue
        free |= ancestor_close([b])
        if len(free)>=size:break
    return free
def subproblem(inc,free,seconds,seed,k):
    h=model(True,seconds,seed);fixed=np.array(sorted(set(range(N))-free),dtype=int)
    ids=np.r_[(fixed[:,None]*T+np.arange(T)).ravel(),M+(fixed[:,None]*T+np.arange(T)).ravel()].astype(np.int32)
    h.changeColsBounds(len(ids),ids,inc[ids],inc[ids])
    # Hard aggregate budget over mining-time distance, with a penalized elastic excess.
    # Dummy period T for unmined blocks makes additions/removals measurable.
    tinc=periods(inc);idx=[];coef=[];constant=0.
    for b in sorted(free):
        dist0=abs(T-tinc[b]);constant+=dist0
        for t in range(T):idx.append(b*T+t);coef.append(abs(t-tinc[b])-dist0)
    # Scale slack penalty to $1000/block-period; start incumbent remains feasible.
    h.addCol(-1000.,0.,hs.kHighsInf,0,np.array([],np.int32),np.array([],float))
    idx.append(NV);coef.append(-1.)
    h.addRow(-hs.kHighsInf,k-constant,len(idx),np.array(idx,np.int32),np.array(coef,float))
    h.setSolution(NV+1,np.arange(NV+1,dtype=np.int32),np.r_[inc,0.])
    t=time.perf_counter();h.run();elapsed=time.perf_counter()-t
    sol=h.getSolution();s=np.array(sol.col_value[:NV]);info=h.getInfo()
    if not sol.value_valid or not valid(s):s=inc.copy()
    return s,elapsed,dict(status=str(h.getModelStatus()),nodes=info.mip_node_count,subproblem_bound=info.mip_dual_bound)
def search(mode,seed,budget,start_kind='initial'):
    d=np.load(OUT/'prepared.npz');inc=d[start_kind].copy();lp=d['lp'].copy();dual=d['dual'].copy();tau=d['tau'].copy()
    rng=np.random.default_rng(seed);counts=np.zeros(4);means=np.zeros(4);size=180;logs=[];refreshes=[]
    start=time.perf_counter();iteration=0
    while time.perf_counter()-start < budget-1:
        remain=budget-(time.perf_counter()-start)
        kind=iteration if iteration<4 else int(np.argmax(means+0.25*np.sqrt(2*np.log(iteration+1)/(counts+1e-9))))
        score=scores(inc,lp,tau,dual)
        free=choose_neighborhood(kind,size,score,tau,dual,rng,mode=='random_closure')
        old=float(c@inc);candidate,elapsed,info=subproblem(inc,free,min(12.,max(.5,remain-.3)),seed+iteration, max(20,.5*len(free)))
        obj=float(c@candidate);gain=max(0,obj-old)
        if gain>1e-5:inc=candidate
        reward=100*(gain/old)/max(elapsed,.01);counts[kind]+=1;means[kind]+=(reward-means[kind])/counts[kind]
        q=dict(iteration=iteration,kind=int(kind),requested_size=size,actual_size=len(free),seconds=elapsed,objective=float(c@inc),gain=gain,reward=reward,**info)
        logs.append(q);print(mode,seed,json.dumps(q),flush=True)
        # Size adaptation uses observed sub-MIP time, bounded to keep local scope.
        size=int(np.clip(size*(1.2 if elapsed<8 else .85 if elapsed>12.5 else 1),100,500))
        iteration+=1
        # Refresh full relaxation (valid global bound and coherent duals), deterministically every 8 iterations.
        if iteration%8==0 and budget-(time.perf_counter()-start)>6:
            ht=model(False,5);rt=time.perf_counter();ht.run()
            if ht.getModelStatus()==hs.HighsModelStatus.kOptimal:
                so=ht.getSolution();lp=np.array(so.col_value);dual=np.array(so.row_dual)
                xx=lp[:M].reshape(N,T);tau=xx@np.arange(T)+T*(1-xx.sum(1))
            refreshes.append(time.perf_counter()-rt)
    method=mode if start_kind=='initial' else mode+'_reference'
    q=dict(method=method,seed=seed,budget=budget,seconds=time.perf_counter()-start,iterations=iteration,refresh_seconds=refreshes,**check(inc))
    stem=f'{method}_{seed}';(OUT/f'{stem}.json').write_text(json.dumps(q,indent=2));pd.DataFrame(logs).to_csv(OUT/f'{stem}_trace.csv',index=False);np.save(OUT/f'{stem}.npy',inc)
    print('FINAL',json.dumps(q),flush=True)
def summarize():
    rows=[];p=json.loads((OUT/'prepare.json').read_text());ub=p['lp_objective'];start=p['initial']['objective'];ref=p['reference']['objective']
    for f in OUT.glob('*.json'):
        q=json.loads(f.read_text())
        if 'method' not in q:continue
        rows.append({k:q.get(k) for k in ['method','seed','objective','seconds','iterations','max_row_violation','integrality_error']} )
    df=pd.DataFrame(rows);df['lp_gap_pct']=100*(ub-df.objective)/ub;df['gain_vs_start_pct']=100*(df.objective-start)/start;df['gain_vs_reference_pct']=100*(df.objective-ref)/ref
    df.to_csv(OUT/'comparison.csv',index=False);print(df.to_string(index=False))
    best=max([(float(c@np.load(f)),f) for f in OUT.glob('*.npy')],key=lambda a:a[0]);s=np.load(best[1]);x=s[:M].reshape(N,T);y=s[M:].reshape(N,T)
    records=[]
    for b in range(N):
        t=int(periods(s)[b]);records.append(dict(block=b,x=blocks.iloc[b,1],y=blocks.iloc[b,2],z=blocks.iloc[b,3],period=t if t<T else -1,tonnes=ton[b],process_fraction=float(y[b].sum()),waste_fraction=float(x[b].sum()-y[b].sum()),discounted_value=float(c[b*T:(b+1)*T]@x[b]+c[M+b*T:M+(b+1)*T]@y[b])))
    pd.DataFrame(records).to_csv(OUT/'best_schedule.csv',index=False)
    with (OUT/'best_pcpsp.sol').open('w') as f:
        for b in range(N):
            for t in range(T):
                for dest,a in [(0,x[b,t]-y[b,t]),(1,y[b,t])]:
                    if a>1e-9:f.write(f'{b} {dest} {t} {a:.12f}\n')
    yearly=[]
    for t in range(T):yearly.append(dict(period_file=t,year=t+1,blocks=int(np.rint(x[:,t].sum())),mined_tonnes=float(ton@x[:,t]),processed_tonnes=float(ton@y[:,t]),npv=float(c[:M].reshape(N,T)[:,t]@x[:,t]+c[M:].reshape(N,T)[:,t]@y[:,t])))
    pd.DataFrame(yearly).to_csv(OUT/'annual_schedule.csv',index=False)
    (OUT/'best_validation.json').write_text(json.dumps(dict(source=best[1].name,**check(s)),indent=2))
    # Exact primal-dual gap identity for this all-<=, [0,1]-bounded LP.
    # Group totals depend on the chosen optimal dual and on model formulation.
    d=np.load(OUT/'prepared.npz');pi=d['dual'];reduced=c-A.T@pi
    cases={'initial':d['initial'],'reference':d['reference'],'best_found':s}
    gaps=[]
    for name,zz in cases.items():
        rowgap=pi*(up-A@zz);boundgap=np.maximum(reduced,0)*(1-zz)+np.maximum(-reduced,0)*zz
        q=dict(solution=name,lp_gap=float(ub-c@zz),block_once=float(rowgap[:N].sum()),
          precedence=float(rowgap[N:capacity_start-M].sum()),routing_coupling=float(rowgap[capacity_start-M:capacity_start].sum()),
          resource_capacity=float(rowgap[capacity_start:].sum()),variable_bounds=float(boundgap.sum()))
        q['identity_residual']=q['lp_gap']-sum(q[k] for k in ['block_once','precedence','routing_coupling','resource_capacity','variable_bounds']);gaps.append(q)
    pd.DataFrame(gaps).to_csv(OUT/'exact_gap_decomposition.csv',index=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--seconds',type=float,default=90);ap.add_argument('--seed',type=int,default=42);ap.add_argument('--mode',default='guided');ap.add_argument('--start',default='initial',choices=['initial','reference']);a=ap.parse_args()
    if a.stage=='prepare':prepare()
    elif a.stage=='baseline':baseline(a.seconds,a.start)
    elif a.stage=='search':search(a.mode,a.seed,a.seconds,a.start)
    elif a.stage=='summarize':summarize()
