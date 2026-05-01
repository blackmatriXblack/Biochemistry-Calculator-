# Biochemistry Calculator Suite: Technical Documentation

## 1. Overview

The **Biochemistry Calculator Suite** is a comprehensive command-line utility written in C, designed to assist researchers, students, and laboratory professionals with a wide array of biochemical and chemical calculations. The tool integrates three distinct modules into a unified interface (or separate executables, depending on compilation), covering:

1.  **General Chemistry**: Periodic table lookup, molar mass calculation, ideal gas law, and solution chemistry.
2.  **Bioinformatics & Sequence Analysis**: DNA/RNA/Protein sequence analysis, including molecular weight, melting temperature ($T_m$), translation, and isoelectric point (pI) calculation.
3.  **Advanced Biochemistry**: Enzyme kinetics, buffer preparation, radioactive decay, and spectrophotometry (Beer-Lambert Law).

This document details the installation, usage, algorithms, and data sources for the provided `main.c` implementations.

---

## 2. Compilation and Installation

The program requires a standard C compiler (GCC or Clang) and links against the math library (`libm`).

### Prerequisites
*   GCC or Clang compiler
*   Standard C Library (`stdio.h`, `stdlib.h`, `string.h`, `ctype.h`, `math.h`)

### Build Instructions
Compile the source code into an executable. Depending on which version of `main.c` is being used, the binary name may vary.

```bash
gcc -o biochem_calc main.c -lm
```

*Note: The `-lm` flag is mandatory for mathematical functions such as `log()`, `pow()`, and `exp()`.*

---

## 3. Module 1: General Chemistry Calculator

*(Based on the second provided `main.c` file)*

### 3.1. Features
*   **Periodic Table Database**: Contains data for the first 118 elements, including atomic weight, electron configuration, state at STP, and melting/boiling points.
*   **Molar Mass Calculator**: Parses chemical formulas (e.g., `H2SO4`, `Ca(OH)2`) using a recursive descent parser to calculate molecular weight. Supports nested parentheses.
*   **Mass/Moles Conversion**: Interconverts between grams and moles for any valid chemical formula.
*   **Molarity Calculator**: Computes solution molarity given solute mass/moles and solution volume.
*   **Ideal Gas Law**: Solves for Pressure ($P$), Volume ($V$), Moles ($n$), or Temperature ($T$) using $PV = nRT$.
*   **Temperature Converter**: Converts between Celsius, Kelvin, and Fahrenheit.

### 3.2. Usage
Run the executable and follow the interactive menu:

```text
=== Chemistry Calculator ===
1. Periodic Table
2. Element Search
3. Molar Mass Calculator
4. Mass/Moles Conversion
5. Molarity Calculator
6. Ideal Gas Law Calculator (PV=nRT)
7. Temperature Converter
8. Help
0. Exit
Enter your choice: 
```

### 3.3. Algorithmic Details
*   **Formula Parsing**: The `parseFormula` function uses recursion to handle groups within parentheses. It iterates through the string, identifying elements, accumulating subscripts, and recursively calculating the mass of grouped atoms.
*   **Element Lookup**: Linear search is used for symbol/name matching ($O(N)$ where $N=118$), which is computationally negligible.

---

## 4. Module 2: BioSequence Analyzer

*(Based on the first provided `main.c` file)*

### 4.1. Features
This module operates via command-line arguments rather than an interactive menu.

**Syntax:**
```bash
./bioanalyzer <molecule_type> <calculation> <sequence> [options]
```

**Supported Molecule Types:** `DNA`, `RNA`, `Protein`

**Calculations:**
| Calculation | Supported Types | Description |
| :--- | :--- | :--- |
| `length` | All | Returns sequence length. |
| `mw` | All | Calculates Molecular Weight (Da). |
| `gc_content` | DNA, RNA | Percentage of G and C bases. |
| `tm` | DNA | Melting Temperature using Nearest-Neighbor model. |
| `revcomp` | DNA | Reverse complement sequence. |
| `transcribe` | DNA | Converts DNA to RNA (T $\to$ U). |
| `translate` | RNA | Translates RNA to Protein (Standard Genetic Code). |
| `composition` | Protein | Count of each amino acid. |
| `pi` | Protein | Isoelectric Point (iterative pH scan). |
| `gravy` | Protein | Grand Average of Hydropathy Index. |
| `extinction_coefficient` | All | Molar extinction coefficient ($\epsilon$). |

### 4.2. Usage Examples

1.  **Calculate DNA Melting Temperature:**
    ```bash
    ./bioanalyzer DNA tm ATCGATCG -conc 50e-9
    ```
2.  **Translate RNA to Protein:**
    ```bash
    ./bioanalyzer RNA translate AUGGCCUAA
    ```
3.  **Calculate Protein pI:**
    ```bash
    ./bioanalyzer Protein pi MKWVTFISLLFLFSSAYSR
    ```

### 4.3. Algorithmic Details
*   **Melting Temperature ($T_m$)**: Uses the **Nearest-Neighbor (NN) Thermodynamic Model** (SantaLucia, 1998).
    $$ T_m = \frac{\Delta H \cdot 1000}{\Delta S + R \cdot \ln(C_t / 4)} - 273.15 $$
    *   $\Delta H$ and $\Delta S$ are summed from dinucleotide parameters.
    *   Default concentration ($C_t$) is 50 nM unless specified via `-conc`.
*   **Isoelectric Point (pI)**: Calculated by iterating pH from 0.0 to 14.0 in steps of 0.01. The net charge is computed based on the Henderson-Hasselbalch equation for all ionizable groups (N-term, C-term, D, E, C, Y, H, K, R). The pH with the net charge closest to zero is returned.
*   **Translation**: Uses a 3D lookup table `GENETIC_CODE[4][4][4]` mapping RNA codons to amino acids. Translation stops at the first stop codon (`UAA`, `UAG`, `UGA`).

---

## 5. Module 3: Advanced Biochemistry Calculator

*(Based on the third provided `main.c` file)*

### 5.1. Features
This interactive module focuses on laboratory-specific calculations.

1.  **Protein/Peptide MW**: Calculates MW from 1-letter amino acid codes, including terminal H and OH groups.
2.  **Nucleic Acid MW**: Calculates MW for DNA or RNA, including 5' phosphate and 3' hydroxyl ends.
3.  **Concentration Conversion**: Converts between Molarity (M), g/L, % (w/v), and Molality (m). Requires solute MW and solution density for molality conversions.
4.  **Dilution Calculator**: Solves $C_1V_1 = C_2V_2$ for any variable.
5.  **Beer-Lambert Law**: Solves $A = \epsilon \cdot b \cdot c$ for Absorbance, Molar Absorptivity, Path Length, or Concentration.
6.  **pH/pOH Calculations**: Interconverts $[H^+]$, $[OH^-]$, pH, and pOH. Assumes $25^\circ C$ ($pH + pOH = 14$).
7.  **Henderson-Hasselbalch**: Calculates pH, pKa, or the ratio $[A^-]/[HA]$ for buffer systems.
8.  **Radioactive Decay**: Calculates remaining activity, time elapsed, or half-life using $A = A_0 e^{-\lambda t}$.
9.  **Amino Acid pI**: Calculates the isoelectric point for a **single** amino acid based on its $\alpha$-carboxyl, $\alpha$-amino, and side chain pKa values.
10. **Buffer Preparation**: Calculates the mass of weak acid and conjugate base salt required to prepare a specific volume and concentration of buffer at a target pH.
11. **Michaelis-Menten Kinetics**: Calculates reaction velocity $V$ given $V_{max}$, $K_m$, and substrate concentration $[S]$.

### 5.2. Usage
Run the executable and select from the menu:

```text
--------------------------------------------
       Biochemistry Calculator Menu
--------------------------------------------
 1. Calculate Protein/Peptide MW
 2. Calculate Nucleic Acid MW
 3. Concentration Unit Conversion
 4. Dilution Calculation
 5. Beer-Lambert Law Calculation
 6. pH and pOH Calculations
 7. Henderson-Hasselbalch Equation
 8. Radioactive Decay Calculation
 9. Calculate Single Amino Acid pI
10. Basic Buffer Preparation
11. Michaelis-Menten Kinetics
 0. Exit
--------------------------------------------
```

### 5.3. Algorithmic Details
*   **Buffer Prep**: Uses the Henderson-Hasselbalch equation to determine the ratio of base to acid, then solves a system of linear equations to find the individual concentrations required to meet the total buffer strength.
*   **Single AA pI**: Sorts the three relevant pKa values (if applicable) and averages the two that bracket the neutral zwitterionic species.
    *   Acidic AA: Average of two lowest pKas.
    *   Basic AA: Average of two highest pKas.
    *   Neutral AA: Average of $\alpha$-COOH and $\alpha$-NH$_3^+$ pKas.

---

## 6. Data Sources and Constants

*   **Atomic Weights**: Standard IUPAC values.
*   **Nearest-Neighbor Parameters**: SantaLucia, J. (1998). *PNAS*, 95(4), 1460-1465.
*   **Extinction Coefficients**: Pace et al. (1995) for proteins; Standard biochemistry texts for nucleotides.
*   **Hydropathy Scale**: Kyte & Doolittle (1982).
*   **pKa Values**: Standard biochemical averages (e.g., IPC 2004 for protein pI, standard texts for single AA pI).
*   **Gas Constant**: $R = 0.0821 \, \text{L·atm/(mol·K)}$ for Ideal Gas Law; $R = 1.987 \, \text{cal/(mol·K)}$ for Tm calculations.

---

## 7. Limitations and Assumptions

1.  **Ideal Gas Law**: Assumes ideal behavior; no Van der Waals corrections.
2.  **Tm Calculation**: Uses the basic Nearest-Neighbor model without salt correction adjustments beyond the standard formula. Assumes linear, single-stranded DNA.
3.  **Protein pI**: The iterative method (Module 2) approximates pI to 0.01 pH units. It does not account for local structural effects on pKa values.
4.  **Sequence Validation**: Input sequences are cleaned (whitespace removed, uppercased). Invalid characters result in error messages.
5.  **Memory Management**: Dynamic memory allocation is used for sequence transformations (reverse complement, translation). Users must ensure proper freeing if modifying the code, though the current implementation handles this within the scope.
6.  **Single Amino Acid pI**: Module 3 calculates pI for *free* amino acids, not residues within a protein chain (which is handled by Module 2).

---

## 8. Error Handling

*   **Input Validation**: All numeric inputs are checked for validity (e.g., positive concentrations, non-zero volumes).
*   **Sequence Validation**: Characters are checked against allowed sets (ATCG for DNA, AUCG for RNA, ARNDCEQGHILKMFPSTWYV for Protein).
*   **Division by Zero**: Protected in dilution, Beer-Lambert, and Ideal Gas Law calculations.
*   **Unknown Elements/Formulas**: The molar mass parser returns errors for unrecognized symbols or invalid syntax.
