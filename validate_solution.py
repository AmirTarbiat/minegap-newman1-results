"""Independent validator: native MineLib file parsing; no experiment module import."""
from pathlib import Path
import json,argparse
import numpy as np
root=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--solution',default='results/best_pcpsp.sol');ap.add_argument('--expected',default='results/best_validation.json');ap.add_argument('--output',default='results/independent_validation.json');args=ap.parse_args()
lines=(root/'upload/newman1(4).pcpsp').read_text().splitlines()
meta={l.split(':')[0]:l.split(':')[1].strip() for l in lines if ':' in l}
n=int(meta['NBLOCKS']);nt=int(meta['NPERIODS']);nd=int(meta['NDESTINATIONS']);nr=int(meta['NRESOURCE_SIDE_CONSTRAINTS']);rate=float(meta['DISCOUNT_RATE'])
obj_start=lines.index('OBJECTIVE_FUNCTION:')+1;profits=np.zeros((n,nd))
for line in lines[obj_start:obj_start+n]:
    b,*values=line.split();profits[int(b)]=list(map(float,values))
q=np.zeros((n,nd,nr));qstart=lines.index('RESOURCE_CONSTRAINT_COEFFICIENTS:')+1
for line in lines[qstart:]:
    if line=='EOF':break
    b,d,r,val=line.split();q[int(b),int(d),int(r)]=float(val)
lims=[]
for line in lines[lines.index('RESOURCE_CONSTRAINT_LIMITS:')+1:obj_start-1]:
    r,t,sense,*vals=line.split();lims.append((int(r),int(t),sense,list(map(float,vals))))
dest=np.zeros((n,nd,nt))
for line in (root/args.solution).read_text().splitlines():
    b,d,t,val=line.split();dest[int(b),int(d),int(t)]+=float(val)
x=dest.sum(1);cum=np.cumsum(x,axis=1);violations=[]
for line in (root/'upload/newman1(4).prec').read_text().splitlines():
    b,k,*pp=map(int,line.split());assert k==len(pp)
    for p in pp:violations.extend((cum[b]-cum[p]).tolist())
usage=np.einsum('bdr,bdt->rt',q,dest);resource_violation=0.
for r,t,sense,vals in lims:
    a=usage[r,t]
    if sense=='L':resource_violation=max(resource_violation,a-vals[0])
    elif sense=='G':resource_violation=max(resource_violation,vals[0]-a)
    else:resource_violation=max(resource_violation,vals[0]-a,a-vals[1])
obj=float(np.einsum('bd,bdt,t->',profits,dest,(1+rate)**(-np.arange(nt))))
result=dict(objective=obj,precedence_violation=float(max(0,max(violations))),resource_violation=float(resource_violation),
    extraction_integrality_error=float(abs(x-np.rint(x)).max()),unique_extraction_violation=float(max(0,x.sum(1).max()-1)),
    negative_destination_violation=float(max(0,-dest.min())),resource_usage=usage.tolist())
assert result['precedence_violation']<1e-7
assert resource_violation<.05
assert result['extraction_integrality_error']<1e-7
assert result['unique_extraction_violation']<1e-7
expected=json.loads((root/args.expected).read_text())['objective']
assert abs(obj-expected)<.01
result['objective_reconstruction_error']=obj-expected
(root/args.output).write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
