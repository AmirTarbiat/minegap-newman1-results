# Reproduction

Use Python 3.12 and the versions in requirements.txt. Run from the directory containing the scripts and `upload/`.

```bash
python experiment.py --stage prepare
python experiment.py --stage baseline --seconds 120
```

For each seed 42, 43, 44:

```bash
python experiment.py --stage search --mode guided --seed 42 --seconds 90
python experiment.py --stage search --mode random_closure --seed 42 --seconds 90
```

Revised method, independent start (no reference schedule):

```bash
python exchange_search.py --start initial --seconds 30
python learning_search.py --start initial --mode linucb --seed 42 --seconds 60
```

Repeat the last command for seeds 43 and 44. Exchange search stops when no improving feasible single-block move or two-block exchange exists in its enumerated neighborhood.

Learning ablation from a common warm start:

```bash
python exchange_search.py --start reference --seconds 90
python learning_search.py --mode linucb --seed 42 --seconds 60
python learning_search.py --mode uniform --seed 42 --seconds 60
```

Repeat both learning commands for seeds 43 and 44.

Supplemental direct reference-start runs:

```bash
python experiment.py --stage search --mode guided --seed 42 --seconds 90 --start reference
python experiment.py --stage search --mode random_closure --seed 42 --seconds 90 --start reference
python experiment.py --stage baseline --seconds 90 --start reference
```

Summarize and independently validate:

```bash
python experiment.py --stage summarize
python validate_solution.py
python build_report.py
```

The reported values are from retained files, not assertions of deterministic repeatability under wall-clock limits. Keep original outputs before rerunning, since matching filenames are overwritten. `OPENBLAS_NUM_THREADS=1` was used for all experiment commands in this session.
