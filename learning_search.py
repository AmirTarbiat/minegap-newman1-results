"""Exact-dual-gap ranked boundary neighborhoods, with online contextual LinUCB.
The learned model chooses a neighborhood; MIP and independent checks certify feasibility.
Non-learning comparator: uniform random choice of the same neighborhood candidates.
"""
import argparse,json,time
import numpy as np
import pandas as pd
import experiment as e

def exact_score(inc,pi):
    absA=abs(e.A);row_norm=np.asarray(absA.sum(1)).ravel()
    rowgap=np.maximum(0,pi*(e.up-e.A@inc));reduced=e.c-e.A.T@pi
    boundgap=np.maximum(reduced,0)*(1-inc)+np.maximum(-reduced,0)*inc
    vargap=np.asarray(absA.T@(rowgap/row_norm)).ravel()+boundgap
    return vargap[:e.M].reshape(e.N,e.T).sum(1)+vargap[e.M:].reshape(e.N,e.T).sum(1)

def candidates(inc,score,tau,rng,size):
    ts=e.periods(inc);nscore=e.normalized(score);seed=int(rng.choice(e.N,p=(score+1)/(score+1).sum()))
    # 0: adjacent periods, emphasizing their precedence frontiers.
    t=int(np.clip(ts[seed],0,e.T-2));inpair=(ts==t)|(ts==t+1)
    eligible=np.array([all(ts[p]<=t for p in e.pred[b]) or all(ts[q]>=t+1 for q in e.succ[b]) for b in range(e.N)])
    rank0=np.argsort(-(nscore+1.5*inpair+.7*eligible+rng.uniform(0,.1,e.N)))
    # 1: LP disagreement, with dummy time for blocks not extracted.
    rank1=np.argsort(-(abs(ts-tau)/e.T+nscore+rng.uniform(0,.1,e.N)))
    # 2: connected growth without transitive closure inflation.
    chosen=[];seen=set();queue=[seed]
    while queue and len(chosen)<size:
        queue.sort(key=lambda b:nscore[b]+rng.uniform(0,.08),reverse=True);b=int(queue.pop(0))
        if b in seen:continue
        seen.add(b);chosen.append(b);queue.extend([a for a in e.pred[b]+e.succ[b] if a not in seen])
    rank2=chosen+list(np.argsort(-nscore))
    # 3: spatial bench slab intersecting the current high-gap region.
    rank3=np.argsort(abs(e.z-e.z[seed])-.8*nscore+rng.uniform(0,.08,e.N))
    out=[]
    for rank in [rank0,rank1,rank2,rank3]:
        free=set()
        for b in rank:
            free.add(int(b))
            if len(free)>=size:break
        # Add immediate boundary partners, with a strict block-count ceiling.
        border=sorted(set(a for b in free for a in e.pred[b]+e.succ[b])-free,key=lambda b:-nscore[b])
        free.update(border[:max(0,min(300,int(1.3*size))-len(free))]);out.append(free)
    return out

def run(mode,seed,budget=60,start_kind='reference'):
    data=np.load(e.OUT/'prepared.npz');inc=np.load(e.OUT/('exchange_best.npy' if start_kind=='reference' else 'exchange_initial_best.npy'));pi=data['dual'];tau=data['tau'];lp=data['lp']
    rng=np.random.default_rng(seed);dim=7;cov=[np.eye(dim) for _ in range(4)];bvec=[np.zeros(dim) for _ in range(4)]
    trace=[];iteration=0;begin=time.perf_counter();size=160;ub=float(data['ub']);initial=float(e.c@inc)
    while time.perf_counter()-begin<budget-1:
        score=exact_score(inc,pi);free_sets=candidates(inc,score,tau,rng,size);ts=e.periods(inc);gap=ub-float(e.c@inc)
        contexts=[]
        for free in free_sets:
            ids=np.array(sorted(free));contexts.append(np.array([1.,score[ids].sum()/max(gap,1),np.mean(abs(ts[ids]-tau[ids]))/e.T,len(free)/e.N,iteration/30,gap/ub,np.mean((lp[:e.M].reshape(e.N,e.T)[ids]*(1-lp[:e.M].reshape(e.N,e.T)[ids])).sum(1))]))
        if iteration<4:arm=iteration
        elif mode=='uniform':arm=int(rng.integers(4))
        else:
            u=[]
            for k,x in enumerate(contexts):
                inv=np.linalg.inv(cov[k]);u.append(float((inv@bvec[k])@x+.2*np.sqrt(x@inv@x)))
            arm=int(np.argmax(u))
        old=float(e.c@inc);remain=budget-(time.perf_counter()-begin)
        sol,elapsed,info=e.subproblem(inc,free_sets[arm],max(.5,min(4.,remain-.2)),seed+iteration,60)
        # Optimize destinations globally after every accepted mining schedule.
        rerouted=e.route(e.periods(sol))
        if e.valid(rerouted) and e.c@rerouted>e.c@sol:sol=rerouted
        gain=max(0,float(e.c@sol)-old)
        if gain>1e-6:inc=sol
        reward=gain/max(elapsed,.02)/1000.;x=contexts[arm];cov[arm]+=np.outer(x,x);bvec[arm]+=reward*x
        q=dict(iteration=iteration,arm=arm,free_blocks=len(free_sets[arm]),solve_seconds=elapsed,wall_seconds=time.perf_counter()-begin,objective=float(e.c@inc),gain=gain,reward=reward,**info)
        trace.append(q);print(json.dumps(q),flush=True);iteration+=1
        size=int(np.clip(size*(1.1 if elapsed<2 else .9 if elapsed>5 else 1),100,230))
    name=f'linucb_{seed}' if mode=='linucb' else f'boundary_uniform_{seed}'
    method='boundary_linucb' if mode=='linucb' else 'boundary_uniform'
    if start_kind=='initial':name+='_initial';method+='_initial'
    result=dict(method=method,seed=seed,seconds=time.perf_counter()-begin,budget=budget,iterations=iteration,start_objective=initial,**e.check(inc))
    (e.OUT/f'{name}.json').write_text(json.dumps(result,indent=2));pd.DataFrame(trace).to_csv(e.OUT/f'{name}_trace.csv',index=False);np.save(e.OUT/f'{name}.npy',inc);np.savez(e.OUT/f'{name}_learning.npz',cov=np.array(cov),bvec=np.array(bvec))
    print('FINAL',json.dumps(result),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['linucb','uniform'],default='linucb');p.add_argument('--seed',type=int,default=42);p.add_argument('--seconds',type=float,default=60);p.add_argument('--start',choices=['initial','reference'],default='reference');a=p.parse_args();run(a.mode,a.seed,a.seconds,a.start)
