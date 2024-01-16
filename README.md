# OIPs Calculation

This MATLAB code is designed for calculating OIPs (Optical, Ionization, and Piezoelectric Constants) in ETBM (Empirical Tight-Binding Method) with a focus on zinc-blend crystal structures.

## Usage

To obtain OIPs for different materials, follow these steps:

1. In the `parameterscalc_GaAs_AlAs` function, update the input values in lines 33 to 52 based on the specified material. Modify the following parameters:
    - `inputss(1,1)`: Substrate material (Use the provided material code)
    - `inputss(1,2)`: Direction of substrate
    - `inputss(2,1)`: First material
    - `inputss(3,1)`: First interface
    - `inputss(4,1)`: Second material
    - `inputss(5,1)`: Second interface

    Material codes:
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

    Specify elastic constants for three directions in line 36 using the numbers:
    ```
    2 -> D001
    3 -> D110
    4 -> D111
    ```

2. Utilize experimental values for your specific material in lines 22 to 69 in both the `parameters_calculate_AlAs` and `parameters_calculate_GaAs` functions.

3. Optionally, adjust the fitness function in line 8 of either `parameters_calculate_AlAs` or `parameters_calculate_GaAs` to change the importance of optoelectronic features.

Feel free to customize and experiment with the code according to your specific needs.

Note: This code assumes zinc-blend crystal structures.
