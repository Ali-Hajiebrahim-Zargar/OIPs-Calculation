# Estimation Enhancing in Optoelectronic Property: A Novel Approach Using Orbital Interaction Parameters (OIPs) and Tight-Binding

This project is designed to optimize OIPs using the Empirical Tight-Binding Method (ETBM) for zinc-blend crystal structures. The project focuses on optimizing OIPs to predict optoelectronic properties such as bandgap energy and effective mass in materials like GaAs, AlAs, and other related semiconductor compounds.

## Introduction

The accurate estimation of optoelectronic properties is crucial in the design and optimization of semiconductor materials and quantum structures, such as super-lattices and quantum wells. Traditional methods like Density Functional Theory (DFT) face limitations in precisely predicting bandgap energy and electron transfer properties, especially in complex materials. This MATLAB code provides a powerful tool for researchers to optimize OIPs using the ETBM, allowing for more precise predictions of optoelectronic properties. By leveraging a combination of genetic algorithms and customized cost functions, this approach optimizes the OIPs for a variety of zinc-blend crystal materials, enhancing the accuracy of simulations in semiconductor research.

## Usage

To optimize OIPs for different materials, follow these steps:

1. **Update Input Values**: Open the `parameterscalc_GaAs_AlAs` function and update the input values in lines 33 to 52 based on the specified material:

    - `inputss(1,1)`: Substrate material (Use the provided material code)
    - `inputss(1,2)`: Direction of substrate
    - `inputss(2,1)`: First material
    - `inputss(3,1)`: First interface
    - `inputss(4,1)`: Second material
    - `inputss(5,1)`: Second interface

    **Material codes:**
    ```
    1 -> Si
    2 -> Ge
    3 -> GaAs
    4 -> AlAs
    5 -> InAs
    6 -> GaP
    7 -> AlP
    8 -> InP
    9 -> GaSb
    10 -> AlSb
    11 -> InSb
    12 -> ZnSe
    13 -> ZnS
    14 -> ZnTe
    15 -> CdTe
    16 -> HgTe
    ```

    **Elastic constants for three directions (line 36):**
    ```
    2 -> D001
    3 -> D110
    4 -> D111
    ```

2. **Set Experimental Values**: Use experimental values for your specific material in lines 22 to 69 in both the `parameters_calculate_AlAs` and `parameters_calculate_GaAs` functions.

3. **Adjust Fitness Function (Optional)**: Optionally, adjust the fitness function in line 8 of either `parameters_calculate_AlAs` or `parameters_calculate_GaAs` to change the importance of optoelectronic features.

4. **Run the Code**: After setting up the parameters, run the script in MATLAB to perform the optimization:
   ```
   run('parameterscalc_GaAs_AlAs.m')
   ```
Feel free to customize and experiment with the code according to your specific needs.

## Output

The output consists of the optimized OIPs for the chosen material and structure, which can be used for further simulations or research in optoelectronics.

## Citation

If you use this code in your research, please cite the following paper:
```
@article{HAJIEBRAHIMZARGAR2024207817,
title = {Estimation enhancing in optoelectronic property: A novel approach using orbital interaction parameters and tight-binding},
journal = {Micro and Nanostructures},
volume = {189},
pages = {207817},
year = {2024},
issn = {2773-0123},
doi = {https://doi.org/10.1016/j.micrna.2024.207817},
url = {https://www.sciencedirect.com/science/article/pii/S2773012324000669},
author = {Ali {Haji Ebrahim Zargar} and Ali Amini and Ahmad Ayatollahi},
keywords = {Orbital interaction parameter (OIP), Tight-binding (TB)_ genetic algorithm (GA)},
abstract = {This paper advocates for an innovative approach designed for estimating optoelectronic properties of quantum structures utilizing Tight-Binding (TB) theory. Predicated on the comparative analysis between estimated and actual properties, the study strives to validate the efficacy of this proposed technique; focusing notably on the computation of bandgap energy. It is observed that preceding methodologies offered a restricted accuracy when predicting complex structures like super-lattices and quantum wells. To address this gap, we propose a methodology involving three distinct phases using orbital interaction parameters (OIPs) and the TB theory. The research employed Aluminium Arsenide (AlAs) and Gallium Arsenide (GaAs) as the primary bulk materials. Our novel approach introduces a computation framework that first focuses on bulk computation, subsequently expanding to super-lattice structures. The findings of this research demonstrate promising results regarding the accuracy of predicated optoelectronic properties, particularly the cut-off wavelength. This study paves the way for future research, potentially enhancing the precision of the proposed methodology and its application scope within the field of quantum optoelectronics.}
}
```
