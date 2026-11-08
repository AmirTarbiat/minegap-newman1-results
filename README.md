# MineGap — Newman1

Research code and validated scheduling outputs for Newman1, a MineLib
benchmark for precedence-constrained open-pit mine production scheduling.

## Results

| Run | NPV | Solution |
|---|---:|---|
| Best run started from the reference schedule | 24,120,482.97 | [Solution](results/best_pcpsp.sol) |
| Best run started without the reference schedule | 24,102,783.00 | [Solution](results/best_independent_pcpsp.sol) |

Both exported schedules were checked against the original model's
precedence and capacity constraints.

- [Convergence plot](results/convergence.png)
- [Comparison of runs](results/comparison.csv)
- [Methods and results report](Newman1_Results_and_Method.md)

## Reproduction

Install the dependencies in `requirements.txt` and follow [RUN.md](RUN.md).
Download the original inputs from [MineLib Newman1](https://minelib.org/v1/newman1.xhtml)
and place them in `upload/` with the filenames expected by the scripts.
