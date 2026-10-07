# amp-evomd

Physics-based design of antimicrobial peptides by embedding coarse-grained molecular dynamics
in a genetic algorithm.

Adapted from the Evo-MD code of van Hilten, Methorst, Verwei and Risselada
(Sci. Adv. 2023, 9(11), eade8839, DOI 10.1126/sciadv.ade8839), which optimised curvature
sensing peptides. The genetic algorithm, the GROMACS wrappers and the simulation scaffolding
are theirs. What changed is the axis of comparison and the objective.

![Pipeline overview](docs/pipeline.svg)

## What problem this addresses

The original compares one membrane composition held at two areas, and measures how much the
peptide reduces the work needed to stretch the leaflet. That difference is curvature sensing.

Here the axis is composition rather than tension: one bacterial membrane mimic and one host
mimic, both at constant area. The same measured quantity means something different on that
axis, and two objectives are built from it.

| | original | this fork |
|---|---|---|
| compared | two tension states of one membrane | two compositions |
| ensemble | constant area | constant area, unchanged |
| measured | surface tension, insertion depth | surface tension, insertion depth, peptide-lipid interaction energy |
| objective | curvature sensing free energy | induced leaflet strain (damage), or partitioning selectivity (binding) |

### The two objectives

Both read the same `observables.json`, produced by the same pair of production runs. Switching
between them costs no simulation, and re-scoring a finished campaign under the other is a
post-processing step rather than a second campaign. That is the reason they live in one
repository.

**damage** maximises `delta_gamma` on the target membrane, with a hinge penalty on the host.
`delta_gamma = gamma(free) - gamma(with peptide)` at constant area, so a positive value means
the peptide relieves tension in the leaflet it occupies, which is the same statement as the
peptide generating leaflet strain. This is a damage observable rather than a binding one: the
mechanical reading in the source paper is that a peptide which senses leaflet tension also
induces it, and that membrane-binding peptides therefore lower the energetic cost of the
curved interface of a toroidal pore. Curvature selectivity needed the derivative of this
quantity with respect to area; damage does not, so one area state suffices.

**binding** maximises a partitioning proxy on the target relative to the host. The proxy is
configurable because it should be chosen from the calibration run rather than assumed. None of
the options is the quantity one actually wants, which is the partitioning free energy by
thermodynamic integration; that does not fit inside the loop and is reserved for re-ranking the
survivors, exactly as in the source paper.

### What carries over from the source paper, and what does not

The transferable asset is not the fitness function, which is specific to curvature. It is:

1. Physics-based fitness inside a genetic algorithm is affordable only if the observable is a
   cheap readout. In the original it is a single `gmx energy` call on the surface tension, not
   an alchemical calculation. That ceiling constrains every choice of observable here.
2. Maximising a single physical observable returns its degenerate extreme, not a usable
   peptide. The source paper set out to optimise curvature sensing and obtained the optimal
   binder: a sequence so hydrophobic it would aggregate out of solution. The analogous failure
   here differs by objective. Maximising damage reproduces the hydrophobic extreme directly.
   Maximising binding selectivity does not, because hydrophobic insertion raises affinity for
   both membranes and largely cancels in the difference, leaving electrostatic discrimination;
   its extreme is instead a maximally cationic peptide, selective, weakly potent and
   salt-sensitive. A different failure, still a failure.
3. The useful regime is therefore a window on the objective, and the window has to be
   calibrated against peptides whose behaviour is known experimentally. Without that
   calibration the output is a ranking without a scale.

### Relation to Evo-MD-FE

Biriukov and Vacha (2026, reference 2 below) independently adapted the same Evo-MD lineage to a
bacterial membrane, optimising peptides that reduce the line tension of a membrane pore. Their
work overlaps this one and bounds what is left to do here, so the differences are worth stating
precisely rather than discovering later.

They model an effectively infinite toroidal pore as a lipid stripe, take the box size along the
rim as a collective variable, and obtain the line tension from the slope of the free-energy
profile over three umbrella-sampling windows. Their fitness is the relative reduction of that
line tension. The cost is comparable to the two production runs used here, which disposes of the
assumption that pore energetics are too expensive for a genetic algorithm: what was missing was
the collective variable, not the compute.

What does not transfer, and why:

- Their headline result is that the best line-tension reducers carry a negatively charged polar
  face, with lysine almost absent among top performers. That result was obtained with a hybrid
  Martini 2.2 model using polarisable water, chosen specifically because line tension requires
  it. Charge-related conclusions are the ones most sensitive to the water model, and this fork
  keeps Martini 3 with non-polarisable water, so the result is not imported as a prior.
- Their aromatic-terminal motif is strongly positional: cyclic shifts of their best sequence that
  preserve composition and helical-wheel arrangement all performed worse. The mechanism they
  propose is specific to the pore rim, where aromatic ends anchor both helix ends and align the
  peptide parallel to the membrane edge. The geometry here is a flat bilayer with no rim, so the
  selective pressure that produced the motif is absent and the motif is not encoded as a
  constraint. Whether it emerges anyway from an unconstrained run is a result either way.
- They report that two model curvature-sensing peptides showed no appreciable line-tension
  reduction. This is a direct challenge to the premise of the damage objective here, which treats
  induced leaflet strain as a damage proxy. Pore formation and pore stabilisation may simply be
  different steps, weakly correlated or anticorrelated. The calibration run settles it cheaply.

What this fork adds is the step their conclusions leave open: they simulate peptides already
placed at the pore rim, so membrane recruitment and selectivity are not part of the calculation,
and they flag explicitly that whether their anionic peptides can target an anionic membrane
remains to be established. Recruitment against a host membrane is what the binding objective here
measures. The tension between the two steps, cationic for recruitment and anionic for pore-rim
stabilisation, is the interesting part rather than a contradiction to resolve.

Adopted from their methodology, all of it force-field independent: neutral peptide termini,
explicit physiological salt, replicate evaluations per candidate, and the exclusion of isoleucine
on the grounds documented in `user/usersettings.py`.

## Order of work

The calibration run comes first, and it uses no genetic algorithm.

1. **Build and calibrate the membranes.** Each composition needs its own reference
   measurements, taken on the peptide-free membrane: box area, free-membrane surface tension,
   and bilayer thickness. These are not transferable between compositions. The procedure is
   documented at the top of `user/membranes.py`. Until they are filled in, every entry point
   refuses to start.
2. **Measure the calibration set.** Run a few dozen peptides with known MIC and hemolysis
   through the pipeline with `STATE = 1` in `run.py`. Prioritise coverage of the regimes over
   potency: an inactive peptide and a non-selective lytic one are what fix the ends of the
   scale, and a binding-without-lysis control is what tests whether the damage observable is
   measuring what it claims.
3. **Check whether the observables separate the classes.** `calibrate.py` scores those results
   under either objective and reports pairwise concordance against the labels. This is the
   go/no-go: an optimiser cannot recover a signal that is not in the observable.
4. **Only then run a campaign.** `run.py` with `STATE = 0`.

## Installation on an Ubuntu workstation

Two situations, with different first steps. Pick the one that matches the machine.

- **A: you have `sudo`.** Install the system packages, then continue from step 2.
- **B: you do not have `sudo`, and GROMACS comes from a module.** The common lab and cluster
  case. Skip to "Installing without root" below, then continue from step 2.

Everything after step 1 is identical in both cases and needs no root.

### 1A. System packages, with sudo

```bash
sudo apt update
sudo apt install -y build-essential git python3 python3-pip python3-venv \
                    gromacs libopenmpi-dev openmpi-bin
```

Ubuntu's `gromacs` package installs the binary as `gmx` and is fine for this pipeline. If the
lab has its own GROMACS build, use that instead and skip the `gromacs` package: a source build
is typically 2 to 3 times faster because it is compiled for the machine's own instruction set,
which matters a great deal across thousands of evaluations. A source build with MPI usually
installs the binary as `gmx_mpi` rather than `gmx`; see the `GMX` variable below.

### 1B. Installing without root

Nothing in this pipeline needs to be installed system-wide. The three things `apt` would have
provided are a GROMACS binary, an MPI library, and a Python with `venv`, and on a machine with
environment modules the first two normally come from modules while the third is usually already
present.

**Find and load the modules.** The exact names vary by site, so list them rather than guessing:

```bash
module avail 2>&1 | grep -i -E "gromacs|gmx"
module avail 2>&1 | grep -i -E "openmpi|mpich|mpi|python"
```

Load GROMACS **with an explicit version**, and note whether it brings MPI with it:

```bash
module load gromacs/2024.3     # the exact name and version 'module avail' printed
which gmx gmx_mpi
gmx --version | head -20       # or gmx_mpi --version
```

Name the version. A bare `module load gromacs` resolves to whatever the site has set as the
default, and administrators change that default without announcing it. Intending to use one
version is not the same as having used one, which is what the pinning below makes checkable.

The lines to read are `GROMACS version` (anything from 2021 onward works) and `GPU support`,
since CUDA support cuts the per-candidate cost substantially. If the binary is `gmx_mpi`, set
`export GMX=gmx_mpi` later; the pipeline reads that variable rather than assuming a name.

**Put the module loads in a file, not in `.bashrc`.** A campaign that silently picks up a
different GROMACS build partway through is not reproducible, and the difference will not show up
anywhere in the output:

```bash
cat > ~/amp-env.sh <<'EOF'
module purge
module load gromacs/2024.3     # exact version, not a bare 'gromacs'
export GMX=gmx                 # or gmx_mpi
export GMX_VERSION_EXPECT=2024.3

# Ranks and threads sized to this machine, leaving one core free.
eval "$(python3 "$HOME/AMP-evo-MD/tools/cores.py" --export)"
EOF
```

That `eval` sets `AMP_NP` and `OMP_NUM_THREADS`. Run `python3 tools/cores.py` without
`--export` to see what it chose and why.

Source it before every run: `source ~/amp-env.sh`.

`GMX_VERSION_EXPECT` pins the build. Every entry point checks the running GROMACS against it and
stops with a clear message on a mismatch, rather than producing numbers that silently do not
belong with the rest of the campaign. Matching is a prefix test, so `2024` accepts any 2024.x
and `2024.3` does not accept `2024.4`. Leave it unset and the version is recorded but not
enforced.

**If `python3 -m venv` is unavailable.** Some minimal installs ship Python without `venv`. Check
first, since this is usually a non-problem:

```bash
python3 -m venv --help >/dev/null 2>&1 && echo "venv ok" || echo "venv missing"
```

If it is missing, install a complete Python into your home directory with micromamba. It is a
single static binary, needs no root, and does not touch anything outside `~/micromamba`:

```bash
cd ~
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
./bin/micromamba shell init -s bash -p ~/micromamba
source ~/.bashrc
micromamba create -y -n amp python=3.11
micromamba activate amp
```

With micromamba active, use `micromamba activate amp` wherever the steps below say
`source .venv/bin/activate`, and skip the `python3 -m venv` line.

**MPI without root.** `pip install mpi4py` now ships a binary wheel, and that wheel does not
bundle an MPI library: it links dynamically against `libmpi.so.40` (OpenMPI) or `libmpi.so.12`
(MPICH) and picks whichever it finds at import time. So as long as a module provides one of
those, the wheel works with no compiler and no root. Verify after installing:

```bash
python -c "from mpi4py import MPI; print('MPI ok, size =', MPI.COMM_WORLD.Get_size())"
```

If that fails with `libmpi.so.* not found`, load an MPI module (`module load openmpi`) and try
again. If it still fails, the site's MPI has an incompatible soname, and the fallback is to
compile mpi4py against it:

```bash
export MPICC=$(which mpicc)
pip install --no-binary=mpi4py --force-reinstall mpi4py
```

That needs a compiler, which a `gcc` module usually provides. And if none of this works, the
pipeline still runs: `mpi4py` is optional throughout, and both the campaign and the calibration
set fall back to a single process. It is slower, not blocked.

### 2. Python environment

A virtual environment, never the system Python. Nothing here is pinned, but `vermouth` and
`PeptideBuilder` pull in versions that are better kept isolated, and on a shared machine you
would not be able to install them system-wide anyway.

```bash
git clone https://github.com/madsondeluna/AMP-evo-MD.git
cd AMP-evo-MD

python3 -m venv .venv            # or: micromamba activate amp
source .venv/bin/activate

pip install --upgrade pip
pip install numpy matplotlib vermouth PeptideBuilder
pip install mpi4py               # optional, see 1B; everything runs serially without it
```

`vermouth` provides the `martinize2` command that converts an atomistic peptide into its
coarse-grained topology. Confirm everything is reachable:

```bash
python -c "import numpy, matplotlib, PeptideBuilder; print('python deps ok')"
python -c "from mpi4py import MPI; print('MPI ok')" || echo "no MPI: serial only, fine"
martinize2 --version
```

If `pip` cannot reach the network on a locked-down machine, download the wheels on a machine
that can and copy them over:

```bash
# on a machine with network access, matching the target's Python version:
pip download numpy matplotlib vermouth PeptideBuilder mpi4py -d wheels/
# then on the lab machine:
pip install --no-index --find-links wheels/ numpy matplotlib vermouth PeptideBuilder
```

### 3. insane, for building membranes

`insane.py` builds the initial lipid bilayer. It is not on PyPI, and it goes in your home
directory, so this step needs no root either:

```bash
mkdir -p ~/bin
curl -L https://raw.githubusercontent.com/Tsjerk/Insane/master/insane.py -o ~/bin/insane
chmod +x ~/bin/insane
export PATH="$HOME/bin:$PATH"
insane --help | head -5
```

Add that `PATH` line to `~/amp-env.sh` alongside the module loads, so a run always picks up the
same set of tools.

### 4. Verify the installation

The test suite runs without GROMACS, because the scoring layer is a pure function of the
observables file. It should take under a second:

```bash
python -m unittest discover -s tests
```

Then exercise the reporting tools against synthetic data, which checks matplotlib and the
parsers without spending any simulation time:

```bash
python tests/make_fake_progress.py --out /tmp/fake_progress.txt
python tools/ga_report.py /tmp/fake_progress.txt --outdir /tmp/ga_report
```

That writes `performance.png`, `composition.png`, `positional.png`, `progress.csv` and
`summary.json` to `/tmp/ga_report`. If those appear, everything except GROMACS itself is
working.

### 5. Environment variables

| Variable | Purpose | Default |
|---|---|---|
| `GMX` | GROMACS binary name | `gmx` |
| `GMX_VERSION_EXPECT` | pin the GROMACS version; runs stop on a mismatch | unset, version recorded but not enforced |
| `MARTINIZE2` | martinize2 binary | `martinize2` |
| `OMP_NUM_THREADS` | threads per GROMACS process | set it explicitly, see below |

A complete `~/amp-env.sh` for a machine with modules and no root:

```bash
module purge
module load gromacs                 # the exact module name for this machine

export GMX=gmx                      # or gmx_mpi, whichever 'which' found
export GMX_VERSION_EXPECT=2024.3    # stop if a different build gets loaded
export MARTINIZE2=martinize2
export PATH="$HOME/bin:$PATH"       # for insane

source "$HOME/AMP-evo-MD/.venv/bin/activate"
eval "$(python3 "$HOME/AMP-evo-MD/tools/cores.py" --export)"
```

Source it before every run, rather than relying on `.bashrc`.

### Provenance: what actually ran

Every campaign records the GROMACS version and banner, the martinize2 and Python package
versions, and the loaded module list in the header of its progress file. Every candidate records
the GROMACS version in its own `observables.json`, so a single result is self-describing and a
later re-score does not depend on trusting a note about which build was loaded.

`calibrate.py` warns when a set of results was produced by more than one build and reports how
many candidates came from each. This matters because surface tension is sensitive to the build:
a version change shifts the values in a way that looks exactly like a real effect of the
sequences, and nothing else in the output would reveal it.

Check what a finished run used:

```bash
head -40 campaign_damage/output.txt                     # the campaign header
python -c "import json; print(json.load(open('campaign_damage/KWLRK.../observables.json'))['_provenance'])"
```

## Running it

### Step 1: prepare and calibrate the membranes

Nothing else runs until this is done. `user/membranes.py` raises `UncalibratedMembrane` on any
entry point while the reference values are missing, which is deliberate: a plausible placeholder
would propagate into every fitness value in the campaign while looking entirely reasonable.

```bash
source ~/amp-env.sh
cd ~/AMP-evo-MD

# First pass. Prints the plan and the insane command, then stops.
./tools/protocol.sh bacterial "POPE:70,POPG:20,CDL2:10" target
```

Read the plan it prints. It reports the lipid counts, the water count, the box dimensions, and
whether a 30-residue helix fits with enough clearance from its own periodic image. Stopping here
on the first pass is intentional: lipid naming and the `insane` flags differ between Martini
lipid sets, and a wrong `-l` / `-u` split produces an asymmetric bilayer that looks fine until
the surface tension refuses to reach zero.

Once the plan looks right:

```bash
AUTO=1 ./tools/protocol.sh bacterial "POPE:70,POPG:20,CDL2:10" target
```

This builds, minimises, and equilibrates at zero tension under a semi-isotropic barostat, then
measures the equilibrium area per lipid and prints a second plan with corrected lipid counts.
Rebuild with those counts, re-equilibrate, and then run the production reference:

```bash
PASS=2 AUTO=1 ./tools/protocol.sh bacterial "POPE:70,POPG:20,CDL2:10" target
```

Repeat the whole sequence for the host membrane:

```bash
./tools/protocol.sh host "POPC:60,DPSM:10,CHOL:30" host
```

Why two passes. The production runs use a constant-area ensemble, so the box must sit at the
composition's own equilibrium area. The equilibrium area per lipid is only known after
equilibrating, so the lipid count is corrected on a second pass while the box stays fixed at
8.0 x 8.0 nm. The lipid count is the free parameter, not the box, because the two membranes also
need comparable areas: `delta_gamma` is a per-unit-area quantity, and a large area mismatch
would dilute the per-peptide perturbation differently in each membrane.

Cholesterol is the composition where the first pass will be furthest off. Sterols condense their
neighbours, so a linear sum over the seed area-per-lipid values underestimates the real packing.
Expect to need the second pass more than for the bacterial membrane.

Useful overrides:

```bash
BOX=8.0 MAX_RESIDUES=30 SALT=0.15 EQUIL_NS=200 PROD_NS=400 \
    AUTO=1 ./tools/protocol.sh bacterial "POPE:70,POPG:20,CDL2:10" target
```

Check the result before trusting it:

```bash
python tools/membrane_analysis.py single --name bacterial --dir prep/bacterial
python tools/membrane_analysis.py compare --dir prep --names bacterial host
```

Three things must hold, and all three fail quietly if ignored. The surface tension must be zero
within error: the Evo-MD POPC reference sits at 1.93 +/- 0.99 bar, and a value far from zero
means the box area and the equilibrium area disagree, so every peptide would be scored against a
pre-stressed bilayer. The area must have stopped drifting. The leaflets must be symmetric,
because an area mismatch between them shows up as spontaneous curvature and contaminates the
surface tension. The script states each verdict in words rather than leaving it to be read off
a number.

Finally, paste the emitted `MembraneDef` blocks into `MEMBRANES` in `user/membranes.py`, and
copy the equilibrated structures into place:

```bash
cp prep/bacterial/bacterial_eq.gro user/base/insert_peptide/membrane_bacterial.gro
cp prep/bacterial/bacterial.top    user/base/insert_peptide/membrane_bacterial.top
cp prep/host/host_eq.gro           user/base/insert_peptide/membrane_host.gro
cp prep/host/host.top              user/base/insert_peptide/membrane_host.top

python -c "import user.membranes as m; m.require_all_calibrated(); print('membranes calibrated')"
```

### Step 2: measure the calibration set

Put one sequence per line in `seq.txt`, set `STATE = 1` in `run.py`, and set
`WRITE_XTC = True` in `user/usersettings.py`. The trajectory is needed here because the
insertion-depth analysis reads it; the GA campaign later turns it back off, since storing
trajectories for thousands of candidates fills a disk quickly.

With MPI, using the layout `amp-env.sh` computed:

```bash
source ~/amp-env.sh
mpirun -np $AMP_NP python run.py --output_dir calibration_run
```

Without MPI, the same command minus `mpirun`. It runs the sequences one after another in a
single process, which is slow but not blocked, and for a calibration set of a few dozen peptides
that is a perfectly reasonable way to work:

```bash
export OMP_NUM_THREADS=$(python3 tools/cores.py --threads-per-rank 99 --export \
                         | sed -n 's/export OMP_NUM_THREADS=//p')
python run.py --output_dir calibration_run
```

### How many ranks and threads

`tools/cores.py` works this out, but the reasoning is worth knowing, because the obvious answer
is wrong in two ways.

**Count physical cores, not `nproc`.** On a machine with hyperthreading, `nproc` reports twice
the physical core count. GROMACS will place two threads on one physical core where they contend
for the same execution units, which for small coarse-grained systems is usually a net loss. The
script parses `lscpu` for the real figure and falls back to logical CPUs only when it cannot.

**Prefer many ranks over many threads.** There is no single number to set: occupancy is ranks
times threads per rank, and the split matters more than the product. This pipeline evaluates
many independent candidates, so throughput is what counts, not how fast any one simulation
finishes. A roughly 5000-particle Martini system stops scaling after a few threads, because
domain decomposition and communication overhead grow faster than the work does. Fifteen ranks of
two threads therefore completes fifteen candidates in about the time four ranks of eight threads
completes four, on the same thirty cores. The default is two threads per rank for that reason.

On a 32-core workstation the script returns fifteen ranks of two threads, leaving one core free
for the session and one idle because thirty-one does not divide by two. To use that last core,
pass `--threads-per-rank 1`.

```bash
python3 tools/cores.py                       # explain the choice
python3 tools/cores.py --reserve 2           # leave two cores free on a shared machine
python3 tools/cores.py --threads-per-rank 1  # maximum rank count, nothing idle
python3 tools/cores.py --cores 64            # plan for a machine you are not on
```

The one core held back costs a few percent of throughput and keeps the machine usable while a
campaign runs, which on a shared workstation is not optional.

### Step 3: the go/no-go check

```bash
python calibrate.py calibration_run --objective damage --labels labels.csv
python calibrate.py calibration_run --objective binding --labels labels.csv
```

`labels.csv` has a header and two columns, `sequence,class`, where class is one of `inactive`,
`selective` or `lytic`:

```csv
sequence,class
GIGKFLHSAKKFGKAFVGEIMNS,selective
GIGAVLKVLTTGLPALISWIKRKRQQ,lytic
```

The script reports a pairwise concordance against those labels and states the reading. This is
the decision point for the whole approach: if the observable does not order peptides whose
behaviour is already known, no optimiser setting downstream recovers a signal that is not there.

### Step 4: run a campaign

Set `STATE = 0` in `run.py`, `WRITE_XTC = False` in `user/usersettings.py`, and pick the
objective:

```python
# user/usersettings.py
OBJECTIVE = "damage"      # or "binding"
```

A campaign runs for days, so it must survive a dropped SSH session. Use `tmux`, which lets you
reattach and watch, rather than `nohup`, which does not:

```bash
tmux new -s campaign
# inside tmux:
source ~/amp-env.sh
cd ~/AMP-evo-MD
mpirun -np $AMP_NP python run.py --output_dir campaign_damage   # drop mpirun if no MPI
# detach with Ctrl-b then d
```

Reattach later from any session with `tmux attach -t campaign`. If `tmux` is not installed and
you cannot install it, `screen` is usually present, and failing both:

```bash
nohup mpirun -np $AMP_NP python run.py --output_dir campaign_damage > campaign.log 2>&1 &
```

Resuming after a crash, a reboot, or a scheduler kill:

```bash
mpirun -np $AMP_NP python run.py --output_dir campaign_damage \
    --restart campaign_damage/checkpoint.dat
```

Without MPI, a full campaign in a single process is likely impractical: the default is 144
candidates over 100 iterations, each candidate two 400 ns runs. Measure one candidate first
(see "Resource planning") and multiply before committing to it. If the serial rate is too slow,
the realistic options are getting MPI working, reducing the population and iteration count, or
moving the campaign to a cluster with a scheduler. The calibration set, which is what actually
has to happen first, is fine serially.

### Step 5: watch it while it runs

The report tool reads the progress file the GA appends to as it goes, so it can be pointed at a
campaign in flight:

```bash
python tools/ga_report.py campaign_damage/output.txt --outdir report_damage
```

It writes `progress.csv`, `summary.json` and three figures, and prints a summary table with a
plain-language reading of the state of the search.

Read fitness and diversity together, never fitness alone. Rising fitness on its own does not
distinguish a search that found a motif from one that collapsed onto a few near-identical
survivors and stopped exploring: the two look identical on a fitness plot. Positional entropy
separates them. A real motif shows entropy falling at the positions that matter while others
stay high; premature convergence pulls every position down together. The script says which it
sees.

Re-score a finished campaign under the other objective without running anything:

```bash
python calibrate.py campaign_damage --objective binding
```

This works because `compute_observables` writes every observable from each pair of production
runs, and the two objectives are pure functions of that file. Comparing the two rankings on one
dataset is a post-processing step, not a second campaign.

### Resource planning

Per candidate, two production runs of 400 ns each on a roughly 8 x 8 x 9 nm coarse-grained
system. The exact wall time depends heavily on the build and whether a GPU is available, so
measure it rather than taking an estimate: run a single sequence through `STATE = 1` and time
it.

```bash
printf 'GIGKFLHSAKKFGKAFVGEIMNS\n' > seq.txt
time mpirun -np 1 python run.py --output_dir timing_test
```

Multiply that by the population size and the iteration count. A default campaign is 144
candidates over 100 iterations, which is roughly 14 400 evaluations before deduplication; the GA
reuses the stored fitness of sequences it has already simulated, so the real count is lower, but
the order of magnitude is what matters when deciding whether the workstation is the right place
to run this at all.

Disk: with `WRITE_XTC = False` each candidate leaves a few megabytes. With it on, budget
hundreds of megabytes per candidate, which is why it is only enabled for the calibration set.

## Troubleshooting

**`UncalibratedMembrane` on startup.** Expected until step 1 is complete. The exception names
which reference values are missing.

**`martinize2: command not found`.** The virtual environment is not active, or `vermouth` is not
installed in it. Run `source .venv/bin/activate` and `pip install vermouth`. If martinize2 lives
somewhere else, set `export MARTINIZE2=/full/path/to/martinize2`.

**`gmx: command not found`.** Either GROMACS is not loaded, or the binary is named something
else. Run `module list` to see what is loaded and `which gmx gmx_mpi` to find the binary, then
`export GMX=` whichever exists. On a machine with modules, `module load gromacs` comes first and
has to be repeated in every new shell, which is why it belongs in `~/amp-env.sh`.

**`ImportError: libmpi.so.40: cannot open shared object file`.** The mpi4py wheel found no MPI
library to link against. Load an MPI module (`module load openmpi`) and retry. If the site's MPI
has a different soname, rebuild mpi4py against it with `export MPICC=$(which mpicc)` followed by
`pip install --no-binary=mpi4py --force-reinstall mpi4py`, which needs a compiler from a `gcc`
module. If none of that is available, run without MPI: both the campaign and the calibration set
fall back to a single process.

**`mpirun: command not found`.** No MPI in the environment. Drop `mpirun -np N` from the command
and run `python run.py ...` directly. Raise `OMP_NUM_THREADS` to give the single GROMACS process
the cores the MPI ranks would have used.

**`Permission denied` writing anywhere.** Everything this pipeline writes goes under the clone
and the directories passed with `--output_dir`, so run it from a path you own, such as
`~/AMP-evo-MD`, and not from a shared or read-only location.

**Surface tension will not reach zero.** The box area and the composition's equilibrium area
disagree. Re-run the second planning pass with the measured area per lipid and rebuild with the
corrected lipid count. Do not proceed: every peptide would be scored against a pre-stressed
bilayer.

**Leaflets reported as asymmetric.** The `insane` `-l` and `-u` flags set the lower and upper
leaflets separately, and a wrong split is the usual cause. Rebuild with matching counts, unless
the asymmetry is intentional for an outer-membrane model.

**`gmx energy` finds no interaction-energy terms.** `energygrps` must be set in
`user/base/production/npt.mdp`. It is enabled in this fork, so a missing term means the mdp was
edited or an older `.edr` is being read. The damage objective does not need those terms, but the
binding objective does, which means a run configured for damage with `energygrps` disabled
cannot be re-scored for binding after the fact.

**The campaign stops improving early.** Look at the diversity panels before changing the
objective. If positional entropy is near zero across all positions, the population has collapsed
and the fix is a higher mutation rate or a larger parent pool. If entropy is still high, the
observable may simply not discriminate, which is what the calibration set exists to establish.

## Limitations to state rather than discover

**Imposed helix.** Martini cannot fold, so secondary structure is assigned, not predicted
(`run.py` writes an all-helix assignment). This was a mild assumption for curvature sensors,
which are helical by construction, and is a strong one for AMPs, many of which are disordered
in solution and fold only on binding. It is the first thing to revisit if the calibration set
is misranked.

**Peptide density.** One peptide in a roughly 200-lipid patch is about P:L 1:200. AMPs act
cooperatively at ratios an order of magnitude higher. `delta_gamma` here is dilute-limit,
per-peptide tension generation; extrapolating it to lysis is an assumption. Reference 2 measured
this dependence directly in their geometry and found the effect saturating at four peptides per
pore rim, which is a reason to take the single-peptide limit as a lower bound rather than a
proxy for the cooperative regime. Vestigial two-peptide support exists upstream
(`ENERGY_select_groups_Protein1` and `Protein2`) if cooperativity becomes central, at double the
cost and with a different interpretation of gamma.

**Threshold.** Generating leaflet strain is necessary for rupture, not sufficient, and nothing
in either objective knows where the rupture threshold lies. Establishing it requires long
unbiased runs or biased sampling on the top candidates, outside the loop.

**Electrostatics.** The inherited settings are `epsilon_r = 15` with reaction-field
electrostatics and non-polarisable water, and the systems carry no salt beyond the counterions
added for neutralisation. For a cationic peptide against an anionic membrane that combination
is the least favourable case. Physiological ionic strength and polarisable water should be
evaluated before any quantitative claim.

**Force field.** This fork keeps Martini 3, inherited from reference 1. Reference 2 chose a
hybrid Martini 2.2 model with polarisable water instead, citing a report that Martini 3 captures
peptide curvature sensing less accurately, and noting that polarisable water is needed to
reproduce membrane line tension. The observable used here is derived from exactly that physics,
so the criticism applies. Keeping Martini 3 is a deliberate choice for continuity with the method
this fork is built on, not evidence that the concern does not hold, and any quantitative claim
should address it.

**Ensemble.** The constant-area coupling in `npt.mdp` is load-bearing. Under a
constant-tension barostat the area would respond instead and gamma would be pinned near its
reference value; the damage observable would vanish. Do not change `Pcoupltype` without
replacing the observable.

**Scope.** Neither membrane includes LPS, so "bacterial" means an inner or cytoplasmic
membrane. Any claim about Gram-negative potency has to account for the outer membrane
separately. Nothing here models the cell wall, efflux, proteolysis or serum binding.

## Changes to the upstream code

Fixes to defects that would otherwise surface silently:

- `charge_correction` could only move the net charge toward zero, so a positive lower bound was
  unsatisfiable; it raised `NameError` on an exactly neutral sequence outside the window, and
  had no iteration guard. Rewritten to move in both directions and give up cleanly when the
  alphabet makes the window unreachable. Covered by `tests/test_charge_correction.py`.
- The `prod-NDX` index file was registered once and shared across membranes, with the resulting
  `ValueError` swallowed. Harmless while both membranes were identical POPC systems, wrong as
  soon as compositions differ: group names survive, atom indices do not. Now one per membrane.
- `mpi4py` was imported unconditionally in `core/system.py` and `core/genetic_algorithm.py`
  despite both guarding on `if MPI is not None`, so the package could not be imported without a
  working MPI stack. Guarded, matching the intent already in the code.
- The `Constraint` class was defined and never assigned to `ptr_constraint_function`, so the
  amphipathicity and charge constraints were dead code. Wired, and extended with a net-charge
  window and a hydrophobic-fraction cap.
- `martinize2` was invoked through an absolute path inside one author's home directory. Now
  resolved from `$MARTINIZE2`, defaulting to the one on `PATH`.
- `energygrps` was commented out in `npt.mdp`, which is why the interaction-energy terms
  referenced in the settings were never computed. Enabled; it costs no extra simulation.

Structural changes:

- Membrane reference quantities moved from four bare globals to a registry of `MembraneDef`
  entries. An entry missing a measurement raises rather than defaulting, because a plausible
  placeholder would propagate into every fitness value in the run while looking reasonable.
- Measurement split from scoring: `compute_observables` writes `observables.json`,
  `fitness_damage` and `fitness_binding` are pure functions of it.
- Arginine enabled in the residue alphabet. Excluding it is defensible for curvature sensors
  and not for AMPs.

One upstream discrepancy is carried forward rather than resolved, and flagged in
`user/membranes.py`: the code derives the insertion reference depth as 92 percent of the
monolayer thickness, a factor that does not appear in the main text of the paper, which quotes
b as the monolayer thickness itself. The code value is kept, since it is what produced the
published results, but it is exposed as `B_FROM_THICKNESS` so it can be checked against
section S1.

Note also that `delta_gamma` is positive for a peptide that relieves tension, following the
internal convention of the original code, where fitness is maximised. The paper plots the
corresponding quantity as negative and inverts the y axis of Fig. 2B. Do not compare values
across the two without accounting for that.

## Requirements

| Component | Needed for | Notes |
|---|---|---|
| Python 3.8 or newer | everything | |
| GROMACS 2021 or newer | all simulation | from `apt`, from a module, or a source build; the binary name is read from `$GMX` |
| `numpy` | observables, analysis | |
| `matplotlib` | figures | analysis only, not needed to run a campaign |
| `vermouth` (`martinize2`) | coarse-graining peptides | |
| `PeptideBuilder` | atomistic helix construction | |
| `insane` | building membranes | not on PyPI, fetched separately |
| `mpi4py` | parallel evaluation | optional; without it the campaign and the calibration set both run single-process |

The Martini 3 force field files are included under `user/base/`, so no separate download is
needed for those.

## Layout

```
core/                        genetic algorithm and GROMACS wrappers, from upstream
user/
  membranes.py               membrane registry and reference measurements
  usersettings.py            run configuration, objective selection
  modules/
    generate_peptide_PB.py   structure generation and coarse-graining
    insert_peptide_nmem.py   insertion into each registered membrane
    production.py            production runs
    compute_observables.py   measurement, writes observables.json
    fitness_damage.py        objective B
    fitness_binding.py       objective A
    compute_fitness.py       dispatch
  base/                      mdp files and Martini 3 force field
tools/
  protocol.sh                membrane preparation, end to end
  membrane.py                plan lipid counts, calibrate the reference constants
  membrane_analysis.py       equilibration verdict and cross-membrane comparison
  ga_report.py               campaign performance, diversity and composition
  mdtools.py                 shared parsing, statistics and plot style
docs/pipeline.svg            the diagram at the top of this file
run.py                       campaign and calibration-set entry points
calibrate.py                 score existing results, check against known classes
tests/                       regression tests, no GROMACS required
calibration/reference_popc/  the original POPC membranes and reference energies, as provenance
```

### Which script to reach for

| Task | Command |
|---|---|
| Size ranks and threads for this machine | `python tools/cores.py` |
| Prepare a membrane | `./tools/protocol.sh <name> "<composition>" <role>` |
| Check a membrane is usable | `python tools/membrane_analysis.py single --name X --dir prep/X` |
| Compare membranes | `python tools/membrane_analysis.py compare --dir prep --names bacterial host` |
| Measure known peptides | `run.py` with `STATE = 1` |
| The go or no-go check | `python calibrate.py <dir> --labels labels.csv` |
| Run a campaign | `run.py` with `STATE = 0` |
| Watch a campaign | `python tools/ga_report.py <dir>/output.txt --outdir report` |
| Re-score under the other objective | `python calibrate.py <dir> --objective binding` |

## References

1. van Hilten, N.; Methorst, J.; Verwei, N.; Risselada, H. J. Physics-based generative model of
   curvature sensing peptides; distinguishing sensors from binders. *Science Advances* **2023**,
   9 (11), eade8839. DOI [10.1126/sciadv.ade8839](https://doi.org/10.1126/sciadv.ade8839)

   The method and the code this fork is built on. Cite this when using the pipeline. Their Evo-MD
   and CNN code and datasets are at [zenodo.org/record/7380110](https://zenodo.org/record/7380110);
   the Evo-MD source is at
   [github.com/nvanhilten/Evo-MD_curvature_sensing](https://github.com/nvanhilten/Evo-MD_curvature_sensing).

2. Biriukov, D.; Vacha, R. Peptide Sequence Features that Reduce Membrane Pore Line Tension.
   *bioRxiv* **2026**, posted 10 September 2026, v1.
   DOI [10.64898/2026.09.06.749768](https://doi.org/10.64898/2026.09.06.749768)

   Preprint, not peer reviewed. The nearest prior work: the same Evo-MD lineage applied to a
   bacterial membrane with a pore line-tension objective. Their design decisions and their open
   questions shaped several choices here; see "Relation to Evo-MD-FE" above. Their code and files
   are to be released on Zenodo following peer review and were not available when this was
   written.
