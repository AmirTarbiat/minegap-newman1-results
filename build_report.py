"""Build the readable pilot report and scientific convergence figure from saved results."""
from pathlib import Path
import json,platform
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent;out=root/'results'
p=json.loads((out/'prepare.json').read_text());ub=p['lp_objective'];ref=p['reference']['objective']
df=pd.read_csv(out/'comparison.csv')
def table(frame):
    f=frame.copy()
    for c in f.select_dtypes('number').columns:f[c]=f[c].map(lambda a:f'{a:,.4f}' if pd.notna(a) else '')
    headers=list(map(str,f.columns));rows=[[str(v).replace('|','/') for v in row] for row in f.itertuples(index=False,name=None)]
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(row)+' |' for row in rows])
metrics=['method','seed','objective','seconds','lp_gap_pct','gain_vs_reference_pct']
summary=df.groupby('method').agg(runs=('objective','size'),mean_objective=('objective','mean'),best_objective=('objective','max'),std_objective=('objective','std'),mean_seconds=('seconds','mean')).reset_index()
summary['mean_lp_gap_pct']=100*(ub-summary.mean_objective)/ub
summary.to_csv(out/'aggregate_comparison.csv',index=False)
best=json.loads((out/'best_validation.json').read_text());validation=json.loads((out/'independent_validation.json').read_text())
gap=pd.read_csv(out/'exact_gap_decomposition.csv');annual=pd.read_csv(out/'annual_schedule.csv')
original=df[df.method.isin(['guided','random_closure','full_mip'])]
independent=df[df.method.str.contains('initial')]
learning=df[df.method.isin(['boundary_linucb','boundary_uniform'])]
supp=df[df.method.isin(['guided_reference','random_closure_reference','full_mip_reference'])]
fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
colors={42:'#2864a0',43:'#cc7441',44:'#3b8b70'}
for seed in [42,43,44]:
    for kind,style in [('linucb','-'),('boundary_uniform','--')]:
        path=out/f'{kind}_{seed}_trace.csv'
        if path.exists():
            tr=pd.read_csv(path);axes[0].step(tr.wall_seconds,tr.objective/1e6,where='post',color=colors[seed],linestyle=style,label=f'{"LinUCB" if kind=="linucb" else "Uniform"}, seed {seed}')
    path=out/f'linucb_{seed}_initial_trace.csv'
    if path.exists():
        tr=pd.read_csv(path);axes[1].step(tr.wall_seconds,tr.objective/1e6,where='post',color=colors[seed],label=f'Seed {seed}')
for ax in axes:
    ax.axhline(ref/1e6,color='#b23c3c',ls=':',label='Published reference');ax.axhline(ub/1e6,color='#333333',ls='-.',label='Global LP bound')
    ax.set_xlabel('Neighborhood-search wall time (s)');ax.set_ylabel('Discounted objective (million)');ax.grid(alpha=.2);ax.legend(fontsize=7)
axes[0].set_title('Common start: reference + feasible exchanges')
axes[1].set_title('Independent LP construction + exchanges + LinUCB')
fig.savefig(out/'convergence.png',dpi=200);fig.savefig(out/'convergence.pdf');plt.close(fig)
text=f'''# Newman1 PCPSP: reproducible pilot and revised method

Date: 6 November 2026. All objectives below were computed, not synthesized.

## Main finding

Best feasible objective found: **{best['objective']:,.6f}**, from `{best['source']}`.
Improvement over the supplied published solution: **{100*(best['objective']-ref)/ref:.6f}%** ({best['objective']-ref:,.6f} objective units).
Gap to the independently solved global LP relaxation: **{100*(ub-best['objective'])/ub:.6f}%**.
This is a feasible improvement over the supplied benchmark, not an optimality certificate or a claim of a new global best-known solution.

## Data interpretation and validation

- 1,060 blocks, 6 periods, 2 destinations, 3,922 precedence arcs.
- The actual `.pcpsp` file contains upper limits only: mining 2,000,000 t/period and processing 1,100,000 t/period. No grade/general side constraints.
- The website also lists minimum throughputs, but adding mandatory mining of 1,400,000 t in all 6 periods is infeasible: the entire model contains only {p['total_tonnes']:,.2f} t, versus 8,400,000 t required. Those minima were not imposed.
- PCPSP uses binary extraction times and continuous destination fractions. A block may be split between waste and processing in its extraction period. It cannot be mined fractionally or processed in a later period in this model.
- Destination coefficients in the `.pcpsp` file already include extraction costs. They are used directly; costs are not subtracted twice.
- Native files index periods 0,...,5. Discounting is `(1.08)^(-t)`, exactly reproducing the published objective.
- Published solution reconstructed: {ref:,.8f}; global LP optimum: {ub:,.8f}, solved in {p['lp_seconds']:.4f} s.
- Uploaded `newman1.LPcpit` is the CPIT LP solution, not a usable PCPSP destination solution. PCPSP was solved independently.
- Equivalent reduced formulation: {p['variables']:,} variables, including {p['binary_variables']:,} binary extraction variables; {p['rows']:,} rows. Let x_bt be extraction and y_bt the processed fraction; 0 <= y_bt <= x_bt. Waste fraction is x_bt-y_bt. Objective coefficients are v_b0*x_bt + (v_b1-v_b0)*y_bt, discounted by period.

## Experiment A: original heuristic architecture

All runs in this table start from the same independent LP-priority topological construction, objective {p['initial']['objective']:,.6f}. The supplied published schedule was not used for those runs.

{table(original[metrics])}

The original method has four ancestor-closed neighborhood types and UCB selection. The random comparator uses a random seed order with the same ancestor-completion rule, adaptive size, trust penalty and LP refresh. It is not a full factorial ablation of all three score terms.

Original score: equal-weight, max-normalized sum of (i) Gini fractionality including the not-mined mass, (ii) displacement from the LP expected time, and (iii) resource-dual pressure weighted by block tonnage and LP routing. For unmined blocks, the dummy period is 6 in zero-based indexing. Requested seed-closure sizes range from 100 to 500; ancestor completion can exceed that request. Never truncate a closed set and silently call it closed. Nominal search budget 90 s, sub-MIP target 12 s, solver MIP tolerance 1e-6. A repeated full LP refresh occurs every 8 completed iterations; because the LP is unchanged, this is redundant overhead rather than a new bound.

Full-MIP baseline: same LP-construction start, one solver seed, nominal 120 s. Its result is solver- and formulation-dependent and does not establish superiority over commercial solvers. Actual durations are reported; solver time-limit checks can overrun, so nominal budgets must not be interpreted as equal measured wall times.

## Why the method was revised

Fractionality, displacement and resource pressure are useful priorities, but their sum is not an exact attribution of the LP gap. For the LP max c'z subject to Az <= b and 0 <= z <= 1, let pi >= 0 be an optimal row dual and r = c - A'pi. Then for any feasible incumbent z:

    U - c'z = pi'(b-Az) + sum_j max(r_j,0)*(1-z_j) + sum_j max(-r_j,0)*z_j.

This exact identity includes precedence, resource, extraction-routing and bound contributions. It requires a valid dual solution of the same global LP. Group totals depend on the selected optimal dual and on the formulation; they are not causal or unique ownership estimates.

{table(gap)}

The revised block score allocates each row contribution among its incident variables using fixed nonnegative weights `abs(A_ij)/sum_k abs(A_ik)`, adds each variable's bound contribution, and sums the variable allocations by block. Its sum equals the gap, up to numerical residual. Allocation of shared-row contributions remains a modeling convention, explicitly stated.

## Experiment B: revised method

Pipeline:
1. Solve the global LP and retain a valid upper bound, all row duals and reduced costs.
2. Construct an extraction schedule with a topological LP-time priority sweep; represent not-mined status by a dummy period.
3. Improve the schedule with exact feasible single-block time moves and two-block time exchanges. Destination fractions are optimized exactly by fractional knapsack within each period. Accept only feasible improving moves.
4. Use the full primal-dual gap identity to rank blocks, recomputing the score after accepted incumbent updates.
5. Generate four bounded neighborhoods: adjacent-period frontiers, LP disagreement, connected graph growth, and spatial bench slabs. Add immediate boundary partners up to a strict ceiling of 300 free blocks. Preserve every precedence and resource row, including rows crossing the fixed/free interface. Ancestor closure is no longer required.
6. Solve the restricted MIP, then reoptimize destination fractions globally for the candidate extraction schedule. Retain only a feasible incumbent with a better original objective.
7. Choose neighborhoods with disjoint online LinUCB models. The learning model predicts improvement rate from 7 contextual features: intercept, share of attributed gap, LP displacement, free-block share, iteration, relative LP gap, and LP mining fractionality. Ridge initialization I, exploration coefficient 0.2. Reward is objective improvement divided by observed sub-MIP solve time and 1,000. It is trained only on evaluations from the current run; no reference optimum or future reward labels are provided to it.

Restricted-MIP parameters: target 4 s, initial core size 160 adapted within 100..230, immediate-neighbor completion to <=300, soft aggregate time-shift budget 60 block-periods with a penalized excess variable at 1,000 objective units per block-period. Destination fractions outside the neighborhood are fixed during the sub-MIP and may be reoptimized in the explicit global-routing postprocessing step. The trust region is soft and does not itself guarantee a runtime bound.

### Independent construction: no reference schedule used

{table(independent[metrics])}

`feasible_exchange_initial` starts from the independent LP construction. `boundary_linucb_initial` starts from that independently improved schedule. Thus these runs demonstrate improvements without importing the published incumbent. Search times in the table exclude the earlier LP and exchange stages; add those times for an end-to-end total.

### Learning ablation from a common warm start

All runs below start from the same reference-plus-exchange schedule. This experiment isolates learned selection versus uniform selection of the same candidate generator; it does not test the benefit of importing a published schedule. Exchange stage: {json.loads((out/'exchange.json').read_text())['seconds']:.4f} s, objective {json.loads((out/'exchange.json').read_text())['objective']:,.6f}.

{table(learning[metrics])}

Three seeds are a pilot, not sufficient evidence for general statistical significance. If uniform selection matches or exceeds LinUCB, report that outcome; the feasibility and objective gain must not be attributed automatically to machine learning.

### Supplemental runs starting directly from the published schedule

{table(supp[metrics])}

## Aggregate results

{table(summary)}

## Best schedule and independent checks

Best native schedule: `results/best_pcpsp.sol`; one row per positive block-destination-period fraction. `best_schedule.csv` has all blocks, including the unmined block (period=-1). `annual_schedule.csv` uses both native zero-based period IDs and human year numbers.

{table(annual)}

An independent validator reads the native `.pcpsp`, `.prec` and output `.sol` files directly, without importing the optimization model. It checks unique binary extraction, cumulative precedence, native destination-resource coefficients, capacity limits and original discounted profit. Objective reconstruction residual: {validation['objective_reconstruction_error']:.12g}. Maximum resource violation: {validation['resource_violation']:.12g} t; maximum precedence violation: {validation['precedence_violation']:.12g}; extraction integrality residual: {validation['extraction_integrality_error']:.12g}.

## Reproduction and limits

Python {p['versions']['python']}, HiGHS {p['versions']['highs']}; each solve uses one thread. Shared Linux environment, AMD EPYC 9V74 host, 9 exposed CPUs. Some experimental batches ran concurrently with other single-threaded batches. Timing is exploratory and should be repeated in an isolated environment before publication; objective and feasibility checks remain directly verifiable. Actual stopping time changes the number of completed neighborhoods, so fixed random seeds alone do not guarantee bit-identical time-limited outcomes.

The algorithm guarantees retention of a feasible nondecreasing incumbent when starting from one. It does not guarantee a strict improvement on every instance or prove global optimality. No minimum-throughput, stockpile, blending or grade-control claims can be made from this instance. The experiments are not large-scale validation, a literature novelty assessment, or evidence that ML dominates non-learning alternatives.

The currently supported contribution is a precisely stated primal-dual gap allocation combined with bounded boundary neighborhoods and an auditable online selection policy. A paper still needs a multi-instance MineLib suite, stronger solvers/baselines, matched budgets, more seeds and separate ablations of score, topology, trust penalty and learning.

Files include all supplied inputs, scripts, seed-level results, solver traces, SHA-256 input checksums in `prepare.json`, the independent validation and convergence plots. See `RUN.md` for commands.

Sources: [MineLib Newman1](https://minelib.org/v1/newman1.xhtml), [MineLib format specification](https://minelib.org/v1/minelib_format.pdf).
'''
(root/'Newman1_Results_and_Method.md').write_text(text)
print(table(summary));print('BEST',best)
