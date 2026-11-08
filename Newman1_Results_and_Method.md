# Newman1 PCPSP: reproducible pilot and revised method

Date: 6 November 2026. All objectives below were computed, not synthesized.

## Main finding

Best feasible objective found: **24,120,482.966701**, from `random_closure_reference_42.npy`.
Improvement over the supplied published solution: **1.953879%** (462,253.084238 objective units).
Gap to the independently solved global LP relaxation: **1.494968%**.
This is a feasible improvement over the supplied benchmark, not an optimality certificate or a claim of a new global best-known solution.

The best independent-start learned run (without importing the published schedule) reached **24,102,782.997761**, from `linucb_42_initial.npy`. Its separately exported and independently validated schedule is `results/best_independent_pcpsp.sol`, with block data in `best_independent_schedule.csv` and native-file validation in `independent_start_native_validation.json`.

In the common-start learning ablation, LinUCB mean = 24,066,431.3521 and uniform mean = 24,060,765.1703; the difference is 5,666.1818 (0.02355%). These three seeds do not establish a reliable learning advantage. The highest overall result above was a reference-start random-closure run, with an actual 122.78-second search; do not present it as a LinUCB or independent-start result.

## Data interpretation and validation

- 1,060 blocks, 6 periods, 2 destinations, 3,922 precedence arcs.
- The actual `.pcpsp` file contains upper limits only: mining 2,000,000 t/period and processing 1,100,000 t/period. No grade/general side constraints.
- The website also lists minimum throughputs, but adding mandatory mining of 1,400,000 t in all 6 periods is infeasible: the entire model contains only 5,619,377.75 t, versus 8,400,000 t required. Those minima were not imposed.
- PCPSP uses binary extraction times and continuous destination fractions. A block may be split between waste and processing in its extraction period. It cannot be mined fractionally or processed in a later period in this model.
- Destination coefficients in the `.pcpsp` file already include extraction costs. They are used directly; costs are not subtracted twice.
- Native files index periods 0,...,5. Discounting is `(1.08)^(-t)`, exactly reproducing the published objective.
- Published solution reconstructed: 23,658,229.88246342; global LP optimum: 24,486,549.02161707, solved in 4.1371 s.
- Uploaded `newman1.LPcpit` is the CPIT LP solution, not a usable PCPSP destination solution. PCPSP was solved independently.
- Equivalent reduced formulation: 12,720 variables, including 6,360 binary extraction variables; 30,964 rows. Let x_bt be extraction and y_bt the processed fraction; 0 <= y_bt <= x_bt. Waste fraction is x_bt-y_bt. Objective coefficients are v_b0*x_bt + (v_b1-v_b0)*y_bt, discounted by period.

## Experiment A: original heuristic architecture

All runs in this table start from the same independent LP-priority topological construction, objective 22,607,937.039746. The supplied published schedule was not used for those runs.

| method | seed | objective | seconds | lp_gap_pct | gain_vs_reference_pct |
| --- | --- | --- | --- | --- | --- |
| random_closure | 43.0000 | 23,268,696.8199 | 89.8084 | 4.9736 | -1.6465 |
| guided | 43.0000 | 23,383,402.1885 | 120.2651 | 4.5051 | -1.1617 |
| guided | 42.0000 | 23,235,396.6770 | 89.7906 | 5.1095 | -1.7873 |
| full_mip |  | 22,607,937.0397 | 120.0617 | 7.6720 | -4.4394 |
| random_closure | 44.0000 | 23,167,255.4594 | 95.7600 | 5.3878 | -2.0753 |
| guided | 44.0000 | 23,531,604.3351 | 95.5463 | 3.8999 | -0.5352 |
| random_closure | 42.0000 | 23,366,858.0956 | 89.3232 | 4.5727 | -1.2316 |

The original method has four ancestor-closed neighborhood types and UCB selection. The random comparator uses a random seed order with the same ancestor-completion rule, adaptive size, trust penalty and LP refresh. It is not a full factorial ablation of all three score terms.

Original score: equal-weight, max-normalized sum of (i) Gini fractionality including the not-mined mass, (ii) displacement from the LP expected time, and (iii) resource-dual pressure weighted by block tonnage and LP routing. For unmined blocks, the dummy period is 6 in zero-based indexing. Requested seed-closure sizes range from 100 to 500; ancestor completion can exceed that request. Never truncate a closed set and silently call it closed. Nominal search budget 90 s, sub-MIP target 12 s, solver MIP tolerance 1e-6. A repeated full LP refresh occurs every 8 completed iterations; because the LP is unchanged, this is redundant overhead rather than a new bound.

Full-MIP baseline: same LP-construction start, one solver seed, nominal 120 s. Its result is solver- and formulation-dependent and does not establish superiority over commercial solvers. Actual durations are reported; solver time-limit checks can overrun, so nominal budgets must not be interpreted as equal measured wall times.

## Why the method was revised

Fractionality, displacement and resource pressure are useful priorities, but their sum is not an exact attribution of the LP gap. For the LP max c'z subject to Az <= b and 0 <= z <= 1, let pi >= 0 be an optimal row dual and r = c - A'pi. Then for any feasible incumbent z:

    U - c'z = pi'(b-Az) + sum_j max(r_j,0)*(1-z_j) + sum_j max(-r_j,0)*z_j.

This exact identity includes precedence, resource, extraction-routing and bound contributions. It requires a valid dual solution of the same global LP. Group totals depend on the selected optimal dual and on the formulation; they are not causal or unique ownership estimates.

| solution | lp_gap | block_once | precedence | routing_coupling | resource_capacity | variable_bounds | identity_residual |
| --- | --- | --- | --- | --- | --- | --- | --- |
| initial | 1,878,611.9819 | 0.0000 | 837,003.6681 | 1,011,106.0035 | 26,377.6816 | 4,124.6287 | -0.0000 |
| reference | 828,319.1392 | 0.0000 | 746,545.9761 | 0.0000 | 32,315.3727 | 49,457.7904 | -0.0000 |
| best_found | 366,066.0549 | 0.0000 | 326,313.2577 | 2,972.0480 | 12,514.4946 | 24,266.2547 | -0.0000 |

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

| method | seed | objective | seconds | lp_gap_pct | gain_vs_reference_pct |
| --- | --- | --- | --- | --- | --- |
| boundary_linucb_initial | 42.0000 | 24,102,782.9978 | 59.0080 | 1.5673 | 1.8791 |
| boundary_linucb_initial | 43.0000 | 24,069,979.6762 | 59.2363 | 1.7012 | 1.7404 |
| boundary_linucb_initial | 44.0000 | 24,080,353.3350 | 59.1577 | 1.6589 | 1.7843 |
| feasible_exchange_initial |  | 23,633,487.0198 | 12.1666 | 3.4838 | -0.1046 |

`feasible_exchange_initial` starts from the independent LP construction. `boundary_linucb_initial` starts from that independently improved schedule. Thus these runs demonstrate improvements without importing the published incumbent. Search times in the table exclude the earlier LP and exchange stages; add those times for an end-to-end total.

### Learning ablation from a common warm start

All runs below start from the same reference-plus-exchange schedule. This experiment isolates learned selection versus uniform selection of the same candidate generator; it does not test the benefit of importing a published schedule. Exchange stage: 5.9963 s, objective 23,736,719.732017.

| method | seed | objective | seconds | lp_gap_pct | gain_vs_reference_pct |
| --- | --- | --- | --- | --- | --- |
| boundary_linucb | 43.0000 | 24,059,125.5801 | 59.1883 | 1.7455 | 1.6945 |
| boundary_linucb | 42.0000 | 24,050,884.3411 | 59.6355 | 1.7792 | 1.6597 |
| boundary_uniform | 44.0000 | 24,099,543.7700 | 59.3581 | 1.5805 | 1.8654 |
| boundary_uniform | 43.0000 | 24,021,489.3495 | 59.2944 | 1.8992 | 1.5354 |
| boundary_uniform | 42.0000 | 24,061,262.3913 | 59.4692 | 1.7368 | 1.7036 |
| boundary_linucb | 44.0000 | 24,089,284.1351 | 59.3060 | 1.6224 | 1.8220 |

Three seeds are a pilot, not sufficient evidence for general statistical significance. If uniform selection matches or exceeds LinUCB, report that outcome; the feasibility and objective gain must not be attributed automatically to machine learning.

### Supplemental runs starting directly from the published schedule

| method | seed | objective | seconds | lp_gap_pct | gain_vs_reference_pct |
| --- | --- | --- | --- | --- | --- |
| full_mip_reference |  | 23,658,229.8825 | 118.4850 | 3.3828 | 0.0000 |
| guided_reference | 42.0000 | 24,082,703.0863 | 89.5240 | 1.6493 | 1.7942 |
| random_closure_reference | 42.0000 | 24,120,482.9667 | 122.7810 | 1.4950 | 1.9539 |

## Aggregate results

| method | runs | mean_objective | best_objective | std_objective | mean_seconds | mean_lp_gap_pct |
| --- | --- | --- | --- | --- | --- | --- |
| boundary_linucb | 3.0000 | 24,066,431.3521 | 24,089,284.1351 | 20,215.5082 | 59.3766 | 1.7157 |
| boundary_linucb_initial | 3.0000 | 24,084,372.0030 | 24,102,782.9978 | 16,766.8347 | 59.1340 | 1.6424 |
| boundary_uniform | 3.0000 | 24,060,765.1703 | 24,099,543.7700 | 39,029.5857 | 59.3739 | 1.7388 |
| feasible_exchange_initial | 1.0000 | 23,633,487.0198 | 23,633,487.0198 |  | 12.1666 | 3.4838 |
| feasible_exchange_reference | 1.0000 | 23,736,719.7320 | 23,736,719.7320 |  | 5.9963 | 3.0622 |
| full_mip | 1.0000 | 22,607,937.0397 | 22,607,937.0397 |  | 120.0617 | 7.6720 |
| full_mip_reference | 1.0000 | 23,658,229.8825 | 23,658,229.8825 |  | 118.4850 | 3.3828 |
| guided | 3.0000 | 23,383,467.7335 | 23,531,604.3351 | 148,103.8400 | 101.8673 | 4.5048 |
| guided_reference | 1.0000 | 24,082,703.0863 | 24,082,703.0863 |  | 89.5240 | 1.6493 |
| random_closure | 3.0000 | 23,267,603.4583 | 23,366,858.0956 | 99,805.8098 | 91.6305 | 4.9780 |
| random_closure_reference | 1.0000 | 24,120,482.9667 | 24,120,482.9667 |  | 122.7810 | 1.4950 |

## Best schedule and independent checks

Best native schedule: `results/best_pcpsp.sol`; one row per positive block-destination-period fraction. `best_schedule.csv` has all blocks, including the unmined block (period=-1). `annual_schedule.csv` uses both native zero-based period IDs and human year numbers.

| period_file | year | blocks | mined_tonnes | processed_tonnes | npv |
| --- | --- | --- | --- | --- | --- |
| 0.0000 | 1.0000 | 401.0000 | 1,996,713.3700 | 946,435.5600 | 8,558,778.7301 |
| 1.0000 | 2.0000 | 357.0000 | 1,998,390.6400 | 816,911.1400 | 7,679,043.2864 |
| 2.0000 | 3.0000 | 274.0000 | 1,475,861.1900 | 1,100,000.0000 | 7,486,664.5347 |
| 3.0000 | 4.0000 | 27.0000 | 145,047.1000 | 145,047.1000 | 395,996.4155 |
| 4.0000 | 5.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 5.0000 | 6.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |

An independent validator reads the native `.pcpsp`, `.prec` and output `.sol` files directly, without importing the optimization model. It checks unique binary extraction, cumulative precedence, native destination-resource coefficients, capacity limits and original discounted profit. Objective reconstruction residual: -7.45058059692e-09. Maximum resource violation: 1.39698386192e-09 t; maximum precedence violation: 0; extraction integrality residual: 0.

## Reproduction and limits

Python 3.12.14, HiGHS 1.15.1; each solve uses one thread. Shared Linux environment, AMD EPYC 9V74 host, 9 exposed CPUs. Some experimental batches ran concurrently with other single-threaded batches. Timing is exploratory and should be repeated in an isolated environment before publication; objective and feasibility checks remain directly verifiable. Actual stopping time changes the number of completed neighborhoods, so fixed random seeds alone do not guarantee bit-identical time-limited outcomes.

The algorithm guarantees retention of a feasible nondecreasing incumbent when starting from one. It does not guarantee a strict improvement on every instance or prove global optimality. No minimum-throughput, stockpile, blending or grade-control claims can be made from this instance. The experiments are not large-scale validation, a literature novelty assessment, or evidence that ML dominates non-learning alternatives.

The currently supported contribution is a precisely stated primal-dual gap allocation combined with bounded boundary neighborhoods and an auditable online selection policy. A paper still needs a multi-instance MineLib suite, stronger solvers/baselines, matched budgets, more seeds and separate ablations of score, topology, trust penalty and learning.

Files include all supplied inputs, scripts, seed-level results, solver traces, SHA-256 input checksums in `prepare.json`, the independent validation and convergence plots. See `RUN.md` for commands.

Sources: [MineLib Newman1](https://minelib.org/v1/newman1.xhtml), [MineLib format specification](https://minelib.org/v1/minelib_format.pdf).
