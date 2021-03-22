# Atrial Low-K code

This repository contains the code needed to reproduce the figures and simulations for the paper _"Immediate and delayed response of simulated human atrial myocytes to clinically-relevant hypokalemia"_.

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

- Figure 1: [`f1.py`](./f1.py)
- Figure 2: [`f2.py`](./f2.py)
- Figure 3: [`f3.py`](./f3.py)
- Figure 4: [`f4.py`](./f4.py)
- Figure 5: [`f5.py`](./f5.py)
- Finally, to generate the non-schematic parts of Figure 6, use [`f6.py`](./f6.py).

These store their figures in directories `figure-1`, `figure-2`, etc.

## Supplementary figures

Several supplementary figures are provided:

- [`s1-vr.py`](./s1-vr.py) creates supplementary figure 1, showing digitised data from [Gelband et al. 1972](https://doi.org/10.1161/01.RES.30.3.293) and [Ten Eick & Singer, 1979](https://doi.org/10.1161/01.RES.44.4.545), and comparing with model predictions.
- [`s2-ap-cat.py`](./s2-ap-cat.py) creates a model comparison figure showing AP and CaT.
- [`s3-relative-contributions.py`](./s3-relative-contributions.py) creates a figure showing relative contributions of all currents.
- [`s4-restitution.py <model>`](./s4-restitution.py) creates a restitution curve for the given model, and plots it alongside data stored in `literature-data`, digitised from [Van Wagoner et al.](https://doi.org/10.1161/01.RES.85.5.428).
- [`s5-ik1-ko.py`](./s5-ik1-ko.py) generates a figure showing the three ways in which external K+ affects IK1.
- [`s6-ik1-x.py`](./s6-ik1-x.py) generates a figure showing the K-dependence of the rectifiation variable in four models.
- [`s7-fast-slow.py <model>`](./s7-fast-slow.py) generates a figure showing short and long-term changes in concentrations, AP, IK1, INaK, and INaCa.
- [`s8-strand.py`](./s8-strand.py) performs strand simulations and creates an overview figure.
- [`s9-dvdtmax.py <model>`](./s9-dvdtmax.py) generates a figure showing maximum dV/dt for all models.
- [`s10-ina-sd.py <model>`](./s10-ina-sd.py) generates a figure with single-cell excitability data (INa availability and SD curves).
- [`s11a-rm-calculations.py`](./s11a-rm-calculations.py) creates a figure showing details of membrane resistance calculations, as do
- [`s11b-rm-ko.py`](./s11b-rm-ko.py), and
- [`s11c-rm-low.py`](./s11c-rm-low.py) (which can be customised by changing the variables `k`, `dt`, and `dv`` inside the script).
- [`s11d-rm.py`](./s11d-rm.py) creates a membrane resistance figure for a chosen model.
- [`s12-iv.py <model>`](./s12-iv.py) generates a figure showing instantaneous IV relationships at low K.

Linux or OS/X users can use the `all-x` scripts to generate figures for all models.
Similarly, the `all` script regenerates all figures in this study.


