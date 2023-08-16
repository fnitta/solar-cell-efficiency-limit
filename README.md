# Efficiency Limit of Solar Cells
This repository contains MATLAB code to calculate the fundamental performance limits of single-junction solar cells with a realistic analysis based on the Tiedje-Yablonovitch model and including defect-assisted Shockley-Read-Hall (SRH) recombination. It calculates key photovoltaic parameters, such as the open circuit voltage, short circuit current density, fill factor, and power conversion efficiency, and While our work focuses on the efficiency limits of multilayer MoS_2, MoSe_2, WS_2, and WSe_2 solar cells under AM 1.5 G illumination, this code can be used to calculate the efficiency limits of any material under any spectrum.

This MATLAB script models the efficiency limit of single-junction solar cells as a function of the absorber layer's thickness. It calculates key photovoltaic parameters, such as the open circuit voltage, short circuit current density, fill factor, and power conversion efficiency, based on the Shockley-Queisser limit. The script also allows for visual representation of these parameters through plots to better understand the efficiency trends with varying thicknesses.

This README serves two purposes:
1. It provides an overview of the MATLAB code we employed to produce the results presented in our paper.
2. It guides users on the necessary modifications and data files required to calculate the efficiency limits of (a) their materials of interest and/or (b) under other spectra.

## Overview of the MATLAB Code
The script is divided into four main parts:

0. Parameter Initialization and Directory Setup
  - Here, the script sets up directory paths for storing output data and figures and initializes various parameters. The code stores the calculated parameters in .txt format, and the plots in both .fig and .jpg format. For the modeling parameters, see Table I: this includes band gap, effective electron mass, and effective hole mass (Lines 41 - 43, respectively). Also, the material's optical constants (n and k) are defined in Lines 55 - 59. To minimize the amount of code to edit, define your material name with the variable Material_Name in Line 38 (this will be used for the plots); then, put your n and k data files in the folder "Data/Material/[Material_Name]" (without the square brackets) and name them "[Material_Name]-n.txt" and "[Material_Name]-k.txt", respectively. See how it is done for the example data (MoS<sub>2</sub>, MoSe<sub>2</sub>, WS<sub>2</sub>, and WSe<sub>2</sub>). Finally, the spectrum (AM 1.5 G illumination in our case) is defined in Lines 72 - 74.
  - Note: Make sure that your k data is defined for all **wavelengths** of your spectrum. If your k data is defined for energy (such as eV or J), use Planck's equation (λ = hc / E) to convert energy into wavelength. For AM 1.5 G, this means having k data for all wavelengths between 280 nm and 4000 nm, but for other spectra, this can be redefined in Lines 78 - 79. This may mean that you will have to extrapolate some of your k data. Make sure that your n data is defined at the band gap energy of your material, converted to wavelength. Also, make sure your n data is defined in terms of wavelengths, not energy.

1. Calculation of Efficiency Limits with Shockley-Queisser (SQ) Model
  - For comparison, the efficiency limit from the Shockley-Queisser (SQ) model is calculated here. The You should not need to edit anything here. This part exports the current density-voltage (or J-V) characteristics in "Shockley-Queisser Model JV Curve.txt" and the open-circuit voltage (Voc), short-circuit current (Jsc), fill factor (FF), and power conversion efficiency (eff) in "Shockley-Queisser Parameters.txt".
  
2. Calculation of Efficiency Limits with Extended Tiedje-Yablonovitch (TY) Model
  - This section performs calculations to determine Voc, Jsc, FF, and eff for a range of thicknesses (with the L_Range variable in Line 103) with our extended Tiedje-Yablonovitch (TY) model, which includes Auger recombination and defect-assisted Shockley-Reed-Hall (SRH) recombination. It employs Equation 11 in Supplementary Note 1 to generate J-V characteristics, and from that, extracts Voc, Jsc, FF, and eff. While it saves these parameters at all thicknesses, it only saves the J-V data for thicknesses of interests (with the L_Interest variable in Line 104). If any thickness in L_Interest is not part of the L_Range list, the code will not run.
  - This section also analyzes the magnitude of each recombination mechanism. It exports the magnitude of the recombination in units of mA cm<sup>-2</sup>, as well as all of the lifetimes for each recombination mechanism in units of seconds.
  - The data tables that are exported for the photovoltaic parameters and recombination loss mechanism magnitude and lifetimes have number of rows corresponding to the number of thicknesses, and number of columns corresponding to the number of 𝜏_SRH values (and in the order of the 𝜏<sub>SRH</sub> values defined in the variable SRH_tau_range). For the J-V data, the first half of the columns correspond to voltages for each 𝜏_SRH value (and in the order of the 𝜏<sub>SRH</sub> values defined in the variable `SRH_tau_range`); the second half of the columns correspond to current densities. The data is exported in Lines 414 - 430.

3. Data Visualization
- The last section is short, and is used to visualize the computed parameters. Data plots are generated for each parameter against the thickness of the material. Each parameter plot is saved in both .fig and .jpg formats for convenience, and the plots are color-coded to distinguish between different parameters.

## Calculating the Efficiency Limits Under Other Conditions
We hope this code will be useful for calculating the efficiency limits under other conditions. Here is how you can modify the code to calculate the efficiency limits for other materials of interest and for other spectra.

### For Other Materials of Interest
(fill in here)

### For Other Spectra
(fill in here)

## Citing Our Work
If you find this model useful for your research, please cite our paper: K. Nassiri Nazif, F.U. Nitta, A. Daus, K.C. Saraswat, E. Pop, ["Efficiency Limit of Transition Metal Dichalcogenide Solar Cells,"](https://arxiv.org/abs/2307.13166) pre-print arXiv:2307.13166 (2023)
