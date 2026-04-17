# loglikelihood

Minimal Fortran tool to:

1. Read observed lightcurve data from `nicer_profile.dat` (2 columns: phase/time, value).
2. Subtract a fixed background (`backgr = 13500.0`).
3. Compute Gaussian log-likelihood for each candidate lightcurve row in a model file.
4. Evaluate the first parseable model row and print its log-likelihood.
5. Write observed, first-model, and best-model lightcurves to separate files for plotting.

Current full mode assumptions:

- Each model row has exactly `76` numeric values.
- First value is the file loglikelihood.
- Next `11` values are model parameters.
- Next `64` values are model lightcurve bins.
- `nicer_profile.dat` must contain `64` valid rows.

The default run evaluates only the first case (`1` case, where each case is 76 values).

## Output files from full run

- `nicer_lightcurve.dat` (2 columns: phase/bin, background-subtracted NICER counts)
- `first_model_lightcurve.dat` (2 columns: same phase/bin, first evaluated model counts)
- `best_model_lightcurve.dat` (2 columns: same phase/bin, best model counts among evaluated cases)

## Build

```zsh
cmake -S . -B build
cmake --build build
```

## Run

Default model file path:

`/home7/jraynau1/nobackup/loglikelihood/640m_parameters_and_phase_amplitudes.dat`

```zsh
./build/loglikelihood
```

Or pass a model file path explicitly:

```zsh
./build/loglikelihood /path/to/your/models.dat
```

Pass number of cases to evaluate (each case = 76 values = 1 loglikelihood + 11 params + 64 bins):

```zsh
./build/loglikelihood /path/to/your/models.dat 10
```

Equivalent explicit full form:

```zsh
./build/loglikelihood --full /path/to/your/models.dat 10
```

Preview mode (prints first 152 numeric values as two rows of 76, without likelihood evaluation):

```zsh
./build/loglikelihood --preview /path/to/your/models.dat
```

## Plot all three lightcurves

Install plotting dependencies if needed:

```zsh
python3 -m pip install numpy matplotlib
```

Generate the comparison figure after a full run:

```zsh
python3 plot_lightcurves.py
```

Optional custom paths:

```zsh
python3 plot_lightcurves.py --nicer nicer_lightcurve.dat --first-model first_model_lightcurve.dat --model best_model_lightcurve.dat --out lightcurves_comparison.png
```
# loglikelihood
