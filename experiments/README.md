

**common.py** — useful constants and procedures to all the experiments.
    It is strongly recommended to use `RESULTS_PATH`, `get_population(city)`
    `get_incidence_filenames(city)`, and other procedures provided there,
    while developing an experiment.

***************

### DTGS conference
* **collect_data.py** — speedup metering script, depending on processes number
* **collect_data2.py** — (k, tpeak_bias, r^2) surface building script
* **draw_speedup.py** — speedup graph drawer

### Autumn 2016 experiments
* **fit_baroyan_rvachev.py** — fitting and graphs for Baroyan-Rvachev model,
    SLSQP implementation. All the parameters defined as constants at the top of
    the file, `Params` class in `main()` function. TODO factor out to **.cfg** file
* **fit_seir.py** — the same for SEIR model

* **methods_comparison.py** — optimization methods comparison in
    Baroyan-Rvachev model implementations
* **size_comparison.py** — optimization methods comparison, while
    (а) iteration number, (b) random seeds parameters are varying;
    Graphs drawing according to the previous iterations, and also
    without (`smooth=false` in `main()` method)
* **slsqp_comparison.py** — failed parameters impact comparison in **SLSQR**
    implementation, provided by `scipy` library
