"""Boundary-feasible time exchanges with exact fractional destination optimization.
Starting from a native PCPSP incumbent; no relaxation of feasibility.
"""
import time,json,argparse
import numpy as np
import experiment as e

delta=e.v[:,1]-e.v[:,0];ratio=delta/e.ton
def value_period(ids,t):
    ids=np.asarray(ids,dtype=int)
    if not len(ids):return 0.
    base=e.v[ids,0].sum();ids=ids[delta[ids]>0]
    ids=ids[np.argsort(-ratio[ids])];tons=e.ton[ids]
    before=np.cumsum(tons)-tons;frac=np.clip((1.1e6-before)/tons,0,1)
    return float((base+delta[ids]@frac)*e.disc[t])
def run(seconds=90,start='reference'):
    data=np.load(e.OUT/'prepared.npz');s=data[start].copy();ts=e.periods(s)
    stem='exchange' if start=='reference' else 'exchange_initial'
    begin=time.perf_counter();log=[];iteration=0;evals=0
    while time.perf_counter()-begin<seconds:
        groups=[set(np.where(ts==t)[0]) for t in range(e.T)]
        mass=np.array([e.ton[list(g)].sum() for g in groups]);pv=np.array([value_period(list(g),t) for t,g in enumerate(groups)])
        best=(0.,None)
        # Single-block moves, including adjacent-year changes and all other years.
        for b in range(e.N):
            old=ts[b]
            if old>=e.T:continue
            low=max([ts[p] for p in e.pred[b]],default=0)
            high=min([ts[q] for q in e.succ[b]],default=e.T-1)
            for new in range(low,min(high,e.T-1)+1):
                if new==old or mass[new]+e.ton[b]>2e6+1e-7:continue
                gain=value_period(list(groups[old]-{b}),old)+value_period(list(groups[new]|{b}),new)-pv[old]-pv[new];evals+=1
                if gain>best[0]+1e-7:best=(gain,[(b,new)])
        # Swap two blocks between periods. Boundary feasibility uses immediate arcs.
        for early in range(e.T-1):
            for late in range(early+1,e.T):
                delays=[a for a in groups[early] if all(ts[q]>=late for q in e.succ[a])]
                advances=[b for b in groups[late] if all(ts[p]<=early for p in e.pred[b])]
                for a in delays:
                    if time.perf_counter()-begin>seconds:break
                    for b in advances:
                        if a in e.pred[b] or b in e.succ[a]:continue
                        diff=e.ton[b]-e.ton[a]
                        if mass[early]+diff>2e6+1e-7 or mass[late]-diff>2e6+1e-7:continue
                        gain=value_period(list((groups[early]-{a})|{b}),early)+value_period(list((groups[late]-{b})|{a}),late)-pv[early]-pv[late];evals+=1
                        if gain>best[0]+1e-7:best=(gain,[(a,late),(b,early)])
        if best[1] is None:break
        for b,t in best[1]:ts[b]=t
        s=e.route(ts);assert e.valid(s),e.check(s)
        log.append(dict(iteration=iteration,seconds=time.perf_counter()-begin,objective=float(e.c@s),gain=best[0],changes=[[int(b),int(t)] for b,t in best[1]],evaluations=evals))
        print(json.dumps(log[-1]),flush=True);iteration+=1
        np.save(e.OUT/f'{stem}_best.npy',s)
    result=dict(method='feasible_exchange_'+start,seconds=time.perf_counter()-begin,iterations=iteration,evaluations=evals,**e.check(s))
    np.save(e.OUT/f'{stem}_best.npy',s);(e.OUT/f'{stem}.json').write_text(json.dumps(result,indent=2));(e.OUT/f'{stem}_trace.json').write_text(json.dumps(log,indent=2));print('FINAL',json.dumps(result),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=90);p.add_argument('--start',choices=['initial','reference'],default='reference');a=p.parse_args();run(a.seconds,a.start)
