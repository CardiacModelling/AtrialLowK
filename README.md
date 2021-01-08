# Atrial Low-K code

This repository contains the code needed to reproduce the figures and simulations for the paper _"Small Decreases in Plasma K+ Prolong Action Potential Repolarization and May Promote Rhythm Disturbances in Human Atria"_.

## Requirements

- Python 3.4+
- Myokit 1.30.1+
- A working OpenCL installation is required for the tissue simulations
- [Installation help for Myokit and OpenCL](http://myokit.org/install)

## Shared code and models

- Shared code is stored in [`shared.py`](./shared.py) and used by several figures.

Models:

- The **Courtemanche et al. 1998** model.
  Original paper doi: [10.1152/ajpheart.1998.275.1.h301](https://doi.org/10.1152/ajpheart.1998.275.1.h301).
  Myokit [model code](models/courtemanche-1998.mmt) was adapted from
  an existing [CellML model](https://models.physiomeproject.org/exposure/0e03bbe01606be5811691f9d5de10b65) in PMR.
- The **Grandi, Pandit, Voigt et al. 2011** model.
  Original paper doi: [10.1161/CIRCRESAHA.111.253955](https://doi.org/10.1161/CIRCRESAHA.111.253955)
  Myokit [model code](models/grandi-2011.mmt) was newly implemented
  and a [CellML model](https://models.physiomeproject.org/e/596) was uploaded to PMR.
- The **Koivumäki, Korhonen et al. 2011** model.
  Original paper doi: [10.1371/journal.pcbi.1001067](https://doi.org/10.1371/journal.pcbi.1001067)
  Myokit [model code](models/koivumaki-2011.mmt) was newly implemented
  and a [CellML model](https://models.physiomeproject.org/e/632) was uploaded to PMR.
- The **Maleckar et al. 2008** model.
  Original paper doi: [10.1016/j.pbiomolbio.2009.01.010](https://doi.org/10.1016/j.pbiomolbio.2009.01.010)
  Myokit [model code](models/maleckar-2008.mmt) was adapted from
  an existing [CellML model](https://models.physiomeproject.org/exposure/bbd802c6a6d6e69b746244f83b4fb89b) in PMR.
- The **Ni et al. 2017** model.
  Original paper doi: [10.3389/fphys.2017.00946](https://doi.org/10.3389/fphys.2017.00946)
  Myokit [model code](models/ni-2017.mmt) was newly implemented
  and a [CellML model](https://models.physiomeproject.org/e/606) was uploaded to PMR.
- The **Nygren et al. 1998** model.
  Original paper doi: [10.1161/01.RES.82.1.63](https://doi.org/10.1161/01.RES.82.1.63)
  Myokit [model code](models/nygren-1998.mmt) was adapted from
  an existing [CellML model](https://models.physiomeproject.org/exposure/ad761ce160f3b4077bbae7a004c229e3).
- The **Voigt, Heijman et al. 2013** model.
  Original paper doi: [10.1016/j.yjmcc.2013.03.011](https://doi.org/10.1016/j.yjmcc.2013.03.011)
  Myokit [model code](models/voigt-heijman-2013.mmt)
  was adapted from [this code](https://github.com/JordiHeijman/MANTA)
  and a [CellML model](https://models.physiomeproject.org/e/634) was uploaded to PMR.

## Main figures

To re-run the single cell simulations, and create the main figures, use:

- Figure 1: [`f1.py`](./f1.py), stores figures in [`figure-1`](./figure-1).
- Figure 2: [`f2.py`](./f2.py), stores figures in [`figure-2`](./figure-2).
- Figure 3: [`f3.py`](./f3.py), stores figures in [`figure-3`](./figure-3).
- Figure 4: [`f4.py`](./f4.py), stores figures in [`figure-4`](./figure-4).

The final figures require multiple scripts to run.
For the fiber simulations:

- [`m6a-strand-run.py <model> <ko> <cvt>`](./m6a-strand-run.py) runs the "1 second" experiments for the given model, potassium level, and CV threshold (the minimum CV at which propagation is said to have occurred).
  For the simulations in the main body of the paper, at baseline potasisum, use `m6a-strand-run.py voigt 5.4 40`.
- For the "15 minute" experiments, run the same script but with the added argument "long", e.g. `m6a-strand-run.py voigt 5.4 40 long`
- These scripts store generated data in [`strand-data`](./strand-data)
- Supplemental figure results are generated and stored in [`sup-figure-m6a-strand-run`](./sup-figure-m6a-strand-run).
- Figure 5: [`f5.py`](./f5.py), stores figures in [`figure-5`](./figure-5).

Finally, to generate the non-schematic parts of Figure 6, use:

- Figure 6: [`f6.py`](./f6.py), stores figures in [`figure-6`](./figure-6).

## Supplementary figures (fixed models)

Several supplementary figures are provided, some of which can run with different models, some of which are fixed to a single model or compare models.

- [`a1-ap-cat.py`](./a1-ap-cat.py`) creates a model comparison figure, stored in [`sup-figure-a1-ap-cat`](./sup-figure-a1-ap-cat).
- [`a2-relative-contributions.py`](./a2-relative-contributions.py`) creates a figure showing relative contributions of all currents, stored in [`sup-figure-a2-relative-contributions`](./sup-figure-a2-relative-contributions).
- [`a3-rm-calculations.py`](./a3-rm-calculations.py) creates a figure showing details of membrane resistance calculations, as do
- [`a4-rm-ko.py`](./a4-rm-ko.py), and
- [`a5-rm-low.py`](./a5-rm-low.py) (which can be customised by changing the variables `k`, `dt`, and `dv`` inside the script).
- [`a6-rm.py`](./a6-rm.py) creates a membrane resistance figure specifically for the Voigt-Heijman model.
- [`a7-sd.py`](./a7-sd.py) creates a strength-duration curve specifically for the Voigt-Heijman model. It requires `m5-sd` to have been run (see below).

## Supplementary figures (model agnostic)

- [`d1-restitution.py <model>`](./d1-restitution.py) creates a restitution curve for the given model, and plots it alongside data stored in `literature-data`, digitised from [Van Wagoner et al.](https://doi.org/10.1161/01.RES.85.5.428).
- [`m1-ko.py <model>`](./m1-ko.py) generates a figure showing short-term changes with low K, for the selected model.
- [`m2-iv.py <model>`](./m2-iv.py) generates a figure showing instantaneous IV relationships at low K.
- [`m3-inak.py <model>`](./m3-inak.py) generates a figure showing short and long-term changes in concentrations, AP, IK1, INaK, and INaCa.
- [`m4-ina.py <model>`](./m4-ina.py) generates a figure showing the `h_inf * j_inf` curve for INa, and indicates the resting potential at various K concentrations.
- [`m5-sd.py`](./m5-sd.py) generates an SD curve (see "main figures" section above).
- [`m6a-strand-run.py`](./m6a-strand-run.py) performs fiber simulations (see "main figures" section above).
- [`m6b-strand-frp.py`](./m6b-strand-frp.py) generates a supplemental figure showing details of the procedure to determine the functional refractory period.
- [`m6c-strand-wl.py`](./m6c-strand-wl.py) generates a supplemental figure showing wavelength of conduction.
- Linux or OS/X users can run supplementary code for all models (and concentrations) using the `all-x` scripts.

