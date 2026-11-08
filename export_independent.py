"""Export the best independent-start LinUCB result, separate from warm starts."""
import json
import numpy as np,pandas as pd
import experiment as e
best=max([(float(e.c@np.load(p)),p) for p in e.OUT.glob('linucb_*_initial.npy')],key=lambda q:q[0]);s=np.load(best[1]);x=s[:e.M].reshape(e.N,e.T);y=s[e.M:].reshape(e.N,e.T);ts=e.periods(s)
with (e.OUT/'best_independent_pcpsp.sol').open('w') as f:
    for b in range(e.N):
        for t in range(e.T):
            for d,a in [(0,x[b,t]-y[b,t]),(1,y[b,t])]:
                if a>1e-9:f.write(f'{b} {d} {t} {a:.12f}\n')
records=[]
for b in range(e.N):records.append(dict(block=b,x=e.blocks.iloc[b,1],y=e.blocks.iloc[b,2],z=e.blocks.iloc[b,3],period=int(ts[b]) if ts[b]<e.T else -1,tonnes=e.ton[b],process_fraction=float(y[b].sum()),waste_fraction=float(x[b].sum()-y[b].sum())))
pd.DataFrame(records).to_csv(e.OUT/'best_independent_schedule.csv',index=False)
(e.OUT/'best_independent_validation.json').write_text(json.dumps(dict(source=best[1].name,**e.check(s)),indent=2))
print(best)
