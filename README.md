# AL_hiPSC_ionic_model

# Introduction

This MATLAB (The MathWorks, Natick, MA) source code provides an implementation of a phenotype specific ionic model for matured and paced atrial-like hiPSC-CMs integrating the atrial specific IKur and IKCa currents

# Software requirements

Matlab R2022a or later versions, any additional Toolbox required

# Main contributors

Dr. Sofia Botti (Euler institute, USI - Department of Mathematics, UNIPV)
Dr. Chiara Bartolucci (Department of Engineering, UNIBO)
Dr. Rolf Krause (Euler institute, USI - Department of Mathematics, UNIDISTANCE)
Dr. Luca Franca Pavarino (Department of Mathematics, UNIPV)
Dr. Stefano Severi (Department of Engineering, UNIBO)

# License
The software is realized with NO WARRANTY and it is licenzed under BSD 3-Clause license


# How to run
- run the main file main_Botti2024.m for the simulation of 10 seconds of the presented ionic model with pacing at 1 Hertz. 
- To obtain a spontaneous beating behaviour, change the stimFlag at line 50 of main_Botti2024.m, and change the IK1 conductance setting gK1 = 0.05 at line 214 of the function file Botti2024.m to run the simulation the 10 second of the presented ionic model without pacing


# Citing this code
You are kindly requested to cite the following manuscript when using this ionic model:
- S. BOTTI, C. BARTOLUCCI, C. ALTOMARE, M. PACI, R. KRAUSE, L. F. PAVARINO, and S. SEVERI (2024) A novel ionic model for matured and paced atrial-like hiPSC-CMs integrating IKur and IKCa currents, _Computers in Biology and Medicine_, To appear;
  
or, in addition, the preprint version:
-S. BOTTI, C. BARTOLUCCI, C. ALTOMARE, M. PACI, R. KRAUSE, L. F. PAVARINO, and S. SEVERI (2024) A novel ionic model for matured and paced atrial-like hiPSC-CMs integrating IKur and IKCa currents. bioRxiv preprint, DOI: 10.1101/2024.01.12.574782

