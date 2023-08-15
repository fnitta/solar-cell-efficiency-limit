# Efficiency Limit of Solar Cells
This repository contains MATLAB code to calculate the fundamental performance limits of single-junction solar cells with a realistic analysis based on the Tiedje-Yablonovitch model and including defect-assisted Shockley-Read-Hall (SRH) recombination. It calculates key photovoltaic parameters, such as the open circuit voltage, short circuit current density, fill factor, and power conversion efficiency, and While our work focuses on the efficiency limits of multilayer MoS_2, MoSe_2, WS_2, and WSe_2 solar cells under AM 1.5 G illumination, this code can be used to calculate the efficiency limits of any material under any spectrum.

This MATLAB script models the efficiency limit of single-junction solar cells as a function of the absorber layer's thickness. It calculates key photovoltaic parameters, such as the open circuit voltage, short circuit current density, fill factor, and power conversion efficiency, based on the Shockley-Queisser limit. The script also allows for visual representation of these parameters through plots to better understand the efficiency trends with varying thicknesses.

This README serves two purposes:
1. It provides an overview of the MATLAB code we employed to produce the results presented in our paper.
2. It guides users on the necessary modifications and data files required to calculate the efficiency limits of (a) their materials of interest and/or (b) under other spectra.

## Overview of the MATLAB Code
The script is divided into four main parts:

0. Parameter Initialization and Directory Setup
- Here, the script sets up directory paths for storing output data and figures and initializes various parameters. The code stores the calculated parameters in .txt format, and the plots in both .fig and .jpg format. For the modeling parameters, see Table I: this includes band gap, effective electron mass, and effective hole mass (Lines 41 - 43, respectively). Also, the material's optical constants (n and k) are defined in Lines 55 - 59. To minimize the amount of code to edit, define your material name with the variable Material_Name in Line 38 (this will be used for the plots); then, put your n and k data files in the folder "Data/Material/[Material_Name]" (without the square brackets) and name them "[Material_Name]-n.txt" and "[Material_Name]-k.txt", respectively. See how the Finally, the spectrum (AM 1.5 G illumination in our case) is defined in Lines 72 - 74.
- Note: Make sure that your k data is defined for all *wavelengths* of your spectrum. If your k data is defined for energy (such as eV or J), use Planck's equation (λ = hc / E) to convert energy into wavelength. For AM 1.5 G, this means having k data for all wavelengths between 280 nm and 4000 nm, but for other spectra, this can be redefined in Lines 78 - 79. This may mean that you will have to extrapolate some of your k data. Make sure that your n data is defined at the band gap energy of your material, converted to wavelength. Also, make sure your n data is defined in terms of wavelengths, not energy.

2. Calculation of Efficiency Limits with Shockley-Queisser (SQ) Model
For comparison, the efficiency limit from the Shockley-Queisser (SQ) model is calculated here.
  
4.
5. 


In this section, the script performs calculations to determine various photovoltaic parameters, including the open circuit voltage (Voc), short circuit current density (Jsc), fill factor (FF), and power conversion efficiency (eff).
The calculations are based on the Shockley-Queisser limit for single-junction solar cells, a theoretical limit for the maximum efficiency.
The results are written into a file named "Shockley-Queisser Parameters.txt" within the specified output directory.
Data Visualization

This section is dedicated to visualizing the computed parameters.
Data plots are generated for each parameter against the thickness of the absorber layer.
Each parameter plot is saved in both .fig and .jpg formats for convenience.
The plots are color-coded to distinguish between different parameters.

## Calculating the Efficiency Limits Under Other Conditions


### For Other Materials of Interest

### For Other Spectra
