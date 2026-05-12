# MetamorphicSulfur
Data and code to reproduce the analyses in Stewart et al. (2026), Metamorphic sulfur release as a driver of sustained cooling and mass extinction, Science Advances

doi:10.5281/zenodo.20146614
[![DOI](https://zenodo.org/badge/1235941055.svg)](https://doi.org/10.5281/zenodo.20146613)

Scripts included are:
--
*.mat or *.m: MATLAB scripts for generating metamorphic C and S fluxes.

SimpleClimateModel.ipynb: Jupyter Notebook for converting metamorphic C and S fluxes into CO2 concentration (carbon cycle model) and effective radiative forcings (climate model emulator).

Data included are:
--
TT_Xflux_MM.csv: Annual fluxes of X (C for carbon, S for sulfur) for a sill exmplacement of TT (Instantaneous or 10yr) using model MM (AR (default), eqm, LR, Ro for carbon; AR, eqm, Fegley, HF (default), Linert, LR for sulfur). For the volcanic-like emissions in Figure 9, use MM = AR_volceq for carbon and MM = HF_volceq for sulfur.
