# AL_hiPSC_ionic_model

## Introduction

This MATLAB (The MathWorks, Natick, MA) source code provides an implementation of a phenotype-specific ionic model for matured and paced atrial-like hiPSC-CMs integrating the atrial-specific IKur and IKCa currents.

## Software Requirements

- MATLAB R2022a or later versions
- Any additional Toolbox required

## Main Contributors

- Dr. Sofia Botti (Euler Institute, USI - Department of Mathematics, UNIPV)
- Dr. Chiara Bartolucci (Department of Engineering, UNIBO)
- Dr. Rolf Krause (Euler Institute, USI - Department of Mathematics, UNIDISTANCE)
- Dr. Luca Franco Pavarino (Department of Mathematics, UNIPV)
- Dr. Stefano Severi (Department of Engineering, UNIBO)

## License

The software is provided with NO WARRANTY and is licensed under the BSD 2-Clause "Simplified" License.

## How to Run

- Run the main file `main_Botti2024.m` for the simulation of 10 seconds of the presented ionic model with pacing at 1 Hertz.
- To obtain a spontaneous beating behavior, change the `stimFlag` at line 50 of `main_Botti2024.m`, and change the IK1 conductance setting `gK1 = 0.05` at line 214 of the function file `Botti2024.m` to run the simulation for 10 seconds of the presented ionic model without pacing.

## Citing This Code

When using this ionic model, please cite the following manuscript:

```tex
@article{Botti2024,
  author = {Sofia Botti and Chiara Bartolucci and Claudia Altomare and Michelangelo Paci and Lucio Barile and Rolf Krause and Luca Franco Pavarino and Stefano Severi}
  title = {A novel ionic model for matured and paced atrial-like human iPSC-CMs integrating IKur and IKCa currents},
  journal = {Computers in Biology and Medicine},
  volume = {180},
  pages = {108899},
  year = {2024},
  issn = {0010-4825},
  doi = {https://doi.org/10.1016/j.compbiomed.2024.108899},
}
}


