# TurboGenius
`TurboGenius` is a python package realizing simple workflows of the SISSA ab-initio quantum Monte Carlo code, `TurboRVB` and also provide useful command-line tools.

`TurboRVB` software package is now composed of the 4 layers:

- `TurboWorkflows` (Complicated workflows for realizing high-throughput calculations with AiiDA) [ToBe]
- `TurboGenius` (Simple workflows and command-line tools)
- `pyturbo` (Python wrappers)
- `TurboRVB` (Quantum Monte Carlo kernel)

`TurboGenius` is the third layer package, which also contains the second layer package `pyturbo` inside as a submodule.

# Beta version
This is a **beta** version!!!! Contact the developers whenever you find bugs. Any suggestion is also welcome!

# Features of `turbogenius`
One can manage any job of `TurboRVB` on python scripts, or on your terminal using the provided command line tool `turbogenius`. 

For python users, several one-to-one corresponding python modules (classes) are provided, i.e., `makefort10.x` -> `Makefort10_genius` class in `makefort10_genius.py`, `turborvb.x` -> `vmc_genius` class in `vmc_genius.py`. `TurboGenius` is designed as a higher layer package that provide several complicated procedures and functions such as fully automatic workflows. `Turbo-Genius` is implemented based on the lower layer packages `pyturbo` and `TurboRVB`. You can see several examples of `TurboGenius` scripts in the `tests` directory. You can also see several simple workflows using `TurboGenius` in the `tests` directory.

`turbogenius` also implements several workflow methods. e.g., `workflow_vmc`, `workflow_lrdmc`. The workflows are implemebted based on the `genius` class memtioned above. You can also write your own workflows by piling up the `genius` blocks, or `workflows` as you want. The workflow method manages the `TurboRVB` jobs, i.e., automatically submits and collects jobs using `filemanager` pacakge. If you want to use `filemanager` pacakge. please contact Kosuke Nakano. You can see several workflow examples in `tests/workflows_examples` directory (Note: Since I am now refactoring file-manager and workflow again, I have not comitted the workflows to the repository.).

# Quick use of `turbogenius`

Installing from source

.. code-block:: console

    git clone git@git-scm.sissa.it:sorella/turbo_genius.git
    cd turbogenius
    pip install -e . or pip install .

Lauching the command line tool. You can easily see what commands are implemented in ``TurboGenius``.

.. code-block:: console

    % turbogenius --help
    Usage: turbogenius [OPTIONS] COMMAND [ARGS]...
    
    Options:
    --help  Show this message and exit.
    
    Commands:
    convertfort10        convertfort10_genius
    convertfort10mol     convertfort10mol_genius
    convertpfaff         readforward_genius
    convertwf            convert wavefunction
    correlated-sampling  correlated_sampling_genius
    lrdmc                lrdmc_genius
    lrdmcopt             lrdmc_genius
    makefort10           makefort10_genius
    prep                 prep_genius
    vmc                  vmc_genius
    vmcopt               vmcopt_genius

You can also see what options are implemented for each command, e.g., vmcopt:

.. code-block:: console

    % turbogenius vmcopt --help
    Usage: turbogenius vmcopt [OPTIONS]
    
    Options:
    -post                 Postprocess
    -r                    Run a program
    -g                    Generate an input file
    -vmcoptsteps INTEGER  Specify vmcoptsteps
    -optwarmup INTEGER    Specify optwarmupsteps
    -steps INTEGER        Specify steps per one iteration
    -bin INTEGER          Specify bin_block
    -warmup INTEGER       Specify warmupblocks
    -nw INTEGER           Specify num_walkers
    -maxtime INTEGER      Specify maxtime
    -optimizer TEXT       Specify optimizer, sr or lr
    -learn FLOAT          Specify learning_rate
    -reg FLOAT            Specify regularization
    -opt_onebody          flag for opt_onebody
    -opt_twobody          flag for opt_twobody
    -opt_det_mat          flag for opt_det_mat
    -opt_jas_mat          flag for opt_jas_mat
    -opt_det_basis_exp    flag for opt_det_basis_exp
    -opt_jas_basis_exp    flag for opt_jas_basis_exp
    -opt_det_basis_coeff  flag for opt_det_basis_coeff
    -opt_jas_basis_coeff  flag for opt_jas_basis_coeff
    -twist                flag for twist_average
    -kpts INTEGER...      kpts, Specify Monkhorst-Pack grids and shifts, [nkx,nky,nkz,kx,ky,kz]
    -plot                 flag for plotting graph
    -log TEXT             logger level, DEBUG, INFO, ERROR
    --help                Show this message and exit.

# Features of `pyturbo`
One can manage any job of `TurboRVB` on python scripts. There are one-to-one corresponding python modules (classes), i.e., `makefort10.x` -> `Makefort10` class in `makefort10.py`, `convertfort10.x` -> `Convertfort10mol` class in `convertfort10mol.py`. `pyturbo` is designed as a lower layer package such that the modules can be used as **components** of higher-level packages. Indeed, the classes are implemented as simple but flexible as possible. Other complicated methods and modules such as fully automatic workflows should be provided at the higher-level packages such as `TurboGenius`. You can see several examples of pyturbo scripts in the `tests` directory.

# Installation
- `git clone this repository`
- `cd turbogenius` 
- `pip install -e .` or `pip install .`

# Examples
Examples are in the `tests` directory.

# Documentation
[ToBe] There is a Read the Docs in the `docs` directory, but still in progress.
You can generate a html file using `sphinx`. Go to the `docs` directory, 
and type `make html`. The document is generated in `docs/_build/html`.
`index.html` is the main page.

# Reference
K. Nakano et. al in prepareation (2022).
