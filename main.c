#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>

// --- Constants ---

// Average Molecular Weights for Amino Acid Residues (in a peptide chain, subtracts H2O)
// Source: https://www.thermofisher.com/us/en/home/references/ambion-tech-support/rna-modifications/rna-molecular-weight-calculations.html
#define MW_A 71.0788
#define MW_R 156.1875
#define MW_N 114.1038
#define MW_D 115.0886
#define MW_C 103.1388
#define MW_Q 128.1292
#define MW_E 129.1155
#define MW_G 57.0519
#define MW_H 137.1411
#define MW_I 113.1594
#define MW_L 113.1594
#define MW_K 128.1420
#define MW_M 131.1926
#define MW_F 147.1766
#define MW_P 97.1167
#define MW_S 87.0782
#define MW_T 101.1051
#define MW_W 186.2132
#define MW_Y 163.1760
#define MW_V 99.1357

// Additional weights for N- and C-termini (for a complete peptide/protein)
#define MW_N_TERM_H 1.008
#define MW_C_TERM_OH 17.007

// Average Molecular Weights for Nucleotide Residues (in a nucleic acid chain)
// Source: https://www.thermofisher.com/us/en/home/references/ambion-tech-support/rna-modifications/rna-molecular-weight-calculations.html
// DNA
#define MW_dA 313.2
#define MW_dT 304.2
#define MW_dC 289.2
#define MW_dG 329.2
// RNA
#define MW_A_RNA 329.2
#define MW_U_RNA 306.2
#define MW_C_RNA 289.2
#define MW_G_RNA 345.2

// Additional weights for 5' monophosphate and 3' hydroxyl ends
#define MW_DNA_5_PRIME_MONOPHOSPHATE 79.0 // Phosphate group minus H2O
#define MW_DNA_3_PRIME_HYDROXYL 17.0 // Hydroxyl group
#define MW_RNA_5_PRIME_MONOPHOSPHATE 79.0 // Phosphate group minus H2O
#define MW_RNA_3_PRIME_HYDROXYL 17.0 // Hydroxyl group

// pKa values for amino acids (used for pI calculation)
// Source: https://www.chem.ucalgary.ca/courses/351/Carey5th/Ch27/ch27-1-4.html
// These are approximate values and can vary slightly depending on the source.
#define PK1_ALPHA_CARBOXYL 2.3
#define PK2_ALPHA_AMINO 9.6
// Side chain pKa values (PKR)
#define PKR_D 3.9 // Aspartic Acid
#define PKR_E 4.1 // Glutamic Acid
#define PKR_C 8.3 // Cysteine
#define PKR_Y 10.5 // Tyrosine
#define PKR_H 6.0 // Histidine
#define PKR_K 10.5 // Lysine
#define PKR_R 12.5 // Arginine
#define PKR_S 13.0 // Serine (sometimes considered neutral or very high) - Using a typical value
#define PKR_T 13.0 // Threonine (sometimes considered neutral or very high) - Using a typical value

// Standard Temperature for some calculations (e.g., pH) - 25°C
#define STANDARD_TEMP 298.15 // Kelvin

// --- Function Declarations ---
void displayMenu();
void calculateProteinMW();
void calculateNucleicAcidMW();
void convertConcentration();
void calculateDilution();
void calculateBeerLambert();
void calculatePH();
void calculateHendersonHasselbalch();
void calculateRadioactiveDecay();
void calculateAminoAcidPI();
void calculateBufferPrep();
void calculateMichaelisMenten();


// --- Main Function ---
int main() {
    int choice;

    printf("Welcome to the Biochemistry Calculator!\n");

    do {
        displayMenu();
        printf("Enter your choice: ");
        // Basic input validation to handle non-integer input
        if (scanf("%d", &choice) != 1) {
            printf("Invalid input. Please enter a number.\n");
            while (getchar() != '\n'); // Clear the input buffer
            choice = -1; // Set choice to an invalid value
        }


        switch (choice) {
            case 1:
                calculateProteinMW();
                break;
            case 2:
                calculateNucleicAcidMW();
                break;
            case 3:
                convertConcentration();
                break;
            case 4:
                calculateDilution();
                break;
            case 5:
                calculateBeerLambert();
                break;
            case 6:
                calculatePH();
                break;
            case 7:
                calculateHendersonHasselbalch();
                break;
            case 8:
                calculateRadioactiveDecay();
                break;
            case 9:
                calculateAminoAcidPI();
                break;
            case 10:
                calculateBufferPrep();
                break;
            case 11:
                calculateMichaelisMenten();
                break;
            case 0:
                printf("Thank you for using the calculator. Goodbye!\n");
                break;
            default:
                // Only print error if choice wasn't -1 (handled by the if)
                if (choice != -1) {
                    printf("Invalid choice. Please try again.\n");
                }
                break;
        }
        printf("\n"); // For cleaner output

    } while (choice != 0);

    return 0;
}

// --- Display Menu ---
void displayMenu() {
    printf("--------------------------------------------\n");
    printf("       Biochemistry Calculator Menu\n");
    printf("--------------------------------------------\n");
    printf(" 1. Calculate Protein/Peptide MW (1-letter code)\n");
    printf(" 2. Calculate Nucleic Acid MW (DNA/RNA, 1-letter code)\n");
    printf(" 3. Concentration Unit Conversion\n");
    printf(" 4. Dilution Calculation (C1V1 = C2V2)\n");
    printf(" 5. Beer-Lambert Law Calculation\n");
    printf(" 6. pH and pOH Calculations\n");
    printf(" 7. Henderson-Hasselbalch Equation\n");
    printf(" 8. Radioactive Decay Calculation\n");
    printf(" 9. Calculate Single Amino Acid pI\n");
    printf("10. Basic Buffer Preparation (Mass Calculation)\n");
    printf("11. Michaelis-Menten Kinetics (Calculate V)\n");
    printf(" 0. Exit\n");
    printf("--------------------------------------------\n");
}

// --- Calculate Protein/Peptide Molecular Weight ---
// Input: Amino acid 1-letter sequence
void calculateProteinMW() {
    char sequence[1000]; // Assume max length 999
    double total_mw = 0.0;
    int seq_len;

    printf("\n--- Protein/Peptide Molecular Weight Calculation ---\n");
    printf("Enter the amino acid 1-letter sequence (e.g.: GAVK): ");
    scanf("%s", sequence);

    seq_len = strlen(sequence);

    if (seq_len == 0) {
        printf("Sequence cannot be empty.\n");
        return;
    }

    // Calculate total MW of amino acid residues in the sequence
    for (int i = 0; i < seq_len; i++) {
        char aa = toupper(sequence[i]); // Convert to uppercase
        switch (aa) {
            case 'A': total_mw += MW_A; break;
            case 'R': total_mw += MW_R; break;
            case 'N': total_mw += MW_N; break;
            case 'D': total_mw += MW_D; break;
            case 'C': total_mw += MW_C; break;
            case 'Q': total_mw += MW_Q; break;
            case 'E': total_mw += MW_E; break;
            case 'G': total_mw += MW_G; break;
            case 'H': total_mw += MW_H; break;
            case 'I': total_mw += MW_I; break;
            case 'L': total_mw += MW_L; break;
            case 'K': total_mw += MW_K; break;
            case 'M': total_mw += MW_M; break;
            case 'F': total_mw += MW_F; break;
            case 'P': total_mw += MW_P; break;
            case 'S': total_mw += MW_S; break;
            case 'T': total_mw += MW_T; break;
            case 'W': total_mw += MW_W; break;
            case 'Y': total_mw += MW_Y; break;
            case 'V': total_mw += MW_V; break;
            default:
                printf("Warning: Invalid 1-letter code '%c' in sequence. Ignoring.\n", sequence[i]);
                // Continue without adding MW
        }
    }

    // Add N-terminal H and C-terminal OH for a complete peptide chain (not cyclic)
    total_mw += MW_N_TERM_H;
    total_mw += MW_C_TERM_OH;

    printf("Approximate molecular weight of peptide/protein sequence: %.2f g/mol\n", total_mw);
}

// --- Calculate Nucleic Acid Molecular Weight ---
// Input: DNA or RNA sequence (1-letter code)
void calculateNucleicAcidMW() {
    char sequence[1000]; // Assume max length 999
    double total_mw = 0.0;
    int seq_len;
    int type_choice;

    printf("\n--- Nucleic Acid Molecular Weight Calculation ---\n");
    printf("Select nucleic acid type:\n");
    printf("1. DNA\n");
    printf("2. RNA\n");
    printf("Enter choice (1 or 2): ");
    // Basic input validation
    if (scanf("%d", &type_choice) != 1 || (type_choice != 1 && type_choice != 2)) {
        printf("Invalid nucleic acid type selection.\n");
        while (getchar() != '\n'); // Clear the input buffer
        return;
    }


    printf("Enter the nucleic acid sequence (e.g.: ATGC or AUGC): ");
    scanf("%s", sequence);

    seq_len = strlen(sequence);

    if (seq_len == 0) {
        printf("Sequence cannot be empty.\n");
        return;
    }

    // Calculate total MW of nucleotide residues in the sequence
    for (int i = 0; i < seq_len; i++) {
        char base = toupper(sequence[i]); // Convert to uppercase

        if (type_choice == 1) { // DNA
            switch (base) {
                case 'A': total_mw += MW_dA; break;
                case 'T': total_mw += MW_dT; break;
                case 'C': total_mw += MW_dC; break;
                case 'G': total_mw += MW_dG; break;
                default:
                    printf("Warning: Invalid 1-letter code '%c' in DNA sequence. Ignoring.\n", sequence[i]);
                    // Continue without adding MW
            }
        } else { // RNA
            switch (base) {
                case 'A': total_mw += MW_A_RNA; break;
                case 'U': total_mw += MW_U_RNA; break;
                case 'C': total_mw += MW_C_RNA; break;
                case 'G': total_mw += MW_G_RNA; break;
                default:
                    printf("Warning: Invalid 1-letter code '%c' in RNA sequence. Ignoring.\n", sequence[i]);
                    // Continue without adding MW
            }
        }
    }

    // Add end weights if sequence is not empty
    if (seq_len > 0) {
        if (type_choice == 1) { // DNA
            total_mw += MW_DNA_5_PRIME_MONOPHOSPHATE;
            total_mw += MW_DNA_3_PRIME_HYDROXYL;
        } else { // RNA
            total_mw += MW_RNA_5_PRIME_MONOPHOSPHATE;
            total_mw += MW_RNA_3_PRIME_HYDROXYL;
        }
    }


    printf("Approximate molecular weight of %s sequence: %.2f g/mol\n", (type_choice == 1) ? "DNA" : "RNA", total_mw);
}

// --- Concentration Unit Conversion ---
// Supports: M <-> g/L <-> % (w/v) <-> molality (mol/kg solvent)
// Requires solute molecular weight and solution density for M <-> molality
void convertConcentration() {
    int from_unit, to_unit;
    double concentration, molecular_weight = 0.0, density = 0.0; // density in g/mL or kg/L

    printf("\n--- Concentration Unit Conversion ---\n");
    printf("Select input concentration unit:\n");
    printf("1. Molarity (M = mol/L solution)\n");
    printf("2. Grams per Liter (g/L)\n");
    printf("3. Percentage (%% w/v = g/100mL solution)\n");
    printf("4. Molality (mol/kg solvent)\n");
    printf("Enter choice (1-4): ");
     if (scanf("%d", &from_unit) != 1 || from_unit < 1 || from_unit > 4) {
        printf("Invalid input unit selection.\n");
        while (getchar() != '\n');
        return;
    }


    printf("Enter the concentration value: ");
    if (scanf("%lf", &concentration) != 1 || concentration < 0) {
        printf("Invalid concentration value. Must be non-negative.\n");
        while (getchar() != '\n');
        return;
    }


    printf("Select output concentration unit:\n");
    printf("1. Molarity (M = mol/L solution)\n");
    printf("2. Grams per Liter (g/L)\n");
    printf("3. Percentage (%% w/v = g/100mL solution)\n");
    printf("4. Molality (mol/kg solvent)\n");
    printf("Enter choice (1-4): ");
     if (scanf("%d", &to_unit) != 1 || to_unit < 1 || to_unit > 4) {
        printf("Invalid output unit selection.\n");
        while (getchar() != '\n');
        return;
    }


    if (from_unit == to_unit) {
        printf("Input and output units are the same. No conversion needed.\n");
        return;
    }

    // Determine if molecular weight is needed
    if (from_unit == 1 || to_unit == 1 || from_unit == 4 || to_unit == 4) {
        printf("Enter the solute molecular weight (g/mol): ");
        if (scanf("%lf", &molecular_weight) != 1 || molecular_weight <= 0) {
            printf("Invalid molecular weight. Must be positive.\n");
             while (getchar() != '\n');
            return;
        }
    }

    // Determine if density is needed (for M <-> molality conversion)
    if ((from_unit == 1 && to_unit == 4) || (from_unit == 4 && to_unit == 1)) {
        printf("Enter the solution density (g/mL or kg/L) (e.g., pure water is ~1.0 g/mL or 1.0 kg/L): ");
         if (scanf("%lf", &density) != 1 || density <= 0) {
            printf("Invalid density. Must be positive.\n");
             while (getchar() != '\n');
            return;
        }
        // Convert density to kg/L if entered in g/mL for molality calculation consistency
        // Assume user enters g/mL or kg/L and they are numerically equivalent for this context.
        // A more robust version would ask for units. Let's assume kg/L for molality calculation.
        // If user enters g/mL, the number is the same as kg/L (1 g/mL = 1 kg/L).
    }


    double conc_in_g_per_l; // Intermediate conversion to g/L

    // Convert input concentration to g/L
    switch (from_unit) {
        case 1: // M to g/L
            conc_in_g_per_l = concentration * molecular_weight;
            break;
        case 2: // g/L to g/L
            conc_in_g_per_l = concentration;
            break;
        case 3: // % (w/v) to g/L (g/100mL * 10)
            conc_in_g_per_l = concentration * 10.0;
            break;
        case 4: // Molality (mol/kg solvent) to g/L (requires density)
            // molality = (mol solute / kg solvent)
            // kg solvent = (kg solution - kg solute)
            // kg solution = Volume (L) * Density (kg/L)
            // mol solute = Molarity * Volume (L)
            // Let's work with a fixed volume, say 1 L solution.
            // Mass of solute in 1 L solution (if we knew M) = M * MW
            // Mass of 1 L solution = Density (kg/L) * 1 L = Density (kg)
            // Mass of solvent in 1 L solution = Density - Mass of solute
            // Molality = Moles solute / (Density - Mass of solute)
            // Molality = M / (Density - M * (MW/1000))  -- if MW is in g/mol, convert to kg/mol
            // From Molality to Molarity is complex algebraically: M = molality * Density / (1 + molality * MW/1000)
            // Let's convert Molality -> Molarity -> g/L
            double molarity_from_molality = (concentration * density) / (1 + concentration * molecular_weight / 1000.0);
            conc_in_g_per_l = molarity_from_molality * molecular_weight;
            break;
    }

    double result;

    // Convert g/L to output concentration unit
    switch (to_unit) {
        case 1: // g/L to M
            result = conc_in_g_per_l / molecular_weight;
            printf("%.4f %s is equal to %.4f M\n", concentration,
                   (from_unit == 1) ? "M" : ((from_unit == 2) ? "g/L" : ((from_unit == 3) ? "%%" : "mol/kg")), result);
            break;
        case 2: // g/L to g/L
            result = conc_in_g_per_l;
             printf("%.4f %s is equal to %.4f g/L\n", concentration,
                   (from_unit == 1) ? "M" : ((from_unit == 2) ? "g/L" : ((from_unit == 3) ? "%%" : "mol/kg")), result);
            break;
        case 3: // g/L to % (w/v) (g/L / 10)
            result = conc_in_g_per_l / 10.0;
             printf("%.4f %s is equal to %.4f %%\n", concentration,
                   (from_unit == 1) ? "M" : ((from_unit == 2) ? "g/L" : ((from_unit == 3) ? "%%" : "mol/kg")), result);
            break;
        case 4: // g/L to Molality (mol/kg solvent) (requires density)
            // g/L -> M -> Molality
            double molarity_for_molality_conv = conc_in_g_per_l / molecular_weight;
             // Molarity -> Molality
            result = (molarity_for_molality_conv * density) / (density - molarity_for_molality_conv * molecular_weight / 1000.0); // MW in kg/mol
            printf("%.4f %s is equal to %.4f mol/kg\n", concentration,
                   (from_unit == 1) ? "M" : ((from_unit == 2) ? "g/L" : ((from_unit == 3) ? "%%" : "mol/kg")), result);
            break;
    }
}


// --- Dilution Calculation C1V1 = C2V2 ---
// Calculate any of C1, V1, C2, or V2 given the other three
void calculateDilution() {
    int unknown;
    double c1 = 0.0, v1 = 0.0, c2 = 0.0, v2 = 0.0;

    printf("\n--- Dilution Calculation (C1V1 = C2V2) ---\n");
    printf("Calculate one variable given the other three.\n");
    printf("1. Calculate C1 (Starting Concentration)\n");
    printf("2. Calculate V1 (Starting Volume)\n");
    printf("3. Calculate C2 (Final Concentration)\n");
    printf("4. Calculate V2 (Final Volume)\n");
    printf("Enter the variable you want to calculate (1-4): ");
     if (scanf("%d", &unknown) != 1 || unknown < 1 || unknown > 4) {
        printf("Invalid choice.\n");
        while (getchar() != '\n');
        return;
    }


    // Prompt for known variables based on the unknown
    if (unknown != 1) {
        printf("Enter C1 (Starting Concentration): ");
         if (scanf("%lf", &c1) != 1 || c1 < 0) { printf("Invalid input. Must be non-negative.\n"); while(getchar() != '\n'); return; }
    }
    if (unknown != 2) {
        printf("Enter V1 (Starting Volume): ");
         if (scanf("%lf", &v1) != 1 || v1 < 0) { printf("Invalid input. Must be non-negative.\n"); while(getchar() != '\n'); return; }
    }
    if (unknown != 3) {
        printf("Enter C2 (Final Concentration): ");
         if (scanf("%lf", &c2) != 1 || c2 < 0) { printf("Invalid input. Must be non-negative.\n"); while(getchar() != '\n'); return; }
    }
    if (unknown != 4) {
        printf("Enter V2 (Final Volume): ");
         if (scanf("%lf", &v2) != 1 || v2 < 0) { printf("Invalid input. Must be non-negative.\n"); while(getchar() != '\n'); return; }
    }

    // Perform calculation and check for division by zero
    double result;
    switch (unknown) {
        case 1: // Calculate C1 = (C2 * V2) / V1
            if (v1 == 0) {
                printf("Error: Starting volume V1 cannot be zero.\n");
                return;
            }
            result = (c2 * v2) / v1;
            printf("Starting Concentration C1 = %.4f\n", result);
            break;
        case 2: // Calculate V1 = (C2 * V2) / C1
            if (c1 == 0) {
                printf("Error: Starting concentration C1 cannot be zero.\n");
                return;
            }
            result = (c2 * v2) / c1;
            printf("Starting Volume V1 = %.4f\n", result);
            break;
        case 3: // Calculate C2 = (C1 * V1) / V2
            if (v2 == 0) {
                printf("Error: Final volume V2 cannot be zero.\n");
                return;
            }
            result = (c1 * v1) / v2;
            printf("Final Concentration C2 = %.4f\n", result);
            break;
        case 4: // Calculate V2 = (C1 * V1) / C2
             if (c2 == 0) {
                printf("Error: Final concentration C2 cannot be zero.\n");
                return;
                }
            result = (c1 * v1) / c2;
            printf("Final Volume V2 = %.4f\n", result);
            break;
    }
}

// --- Beer-Lambert Law Calculation A = epsilon * b * c ---
// Calculate any of A, epsilon, b, or c
// A: Absorbance
// epsilon: Molar Absorptivity
// b: Path Length
// c: Concentration
void calculateBeerLambert() {
    int unknown;
    double absorbance = 0.0, molar_absorptivity = 0.0, path_length = 0.0, concentration = 0.0;

    printf("\n--- Beer-Lambert Law Calculation (A = epsilon * b * c) ---\n");
    printf("Calculate one variable given the other three.\n");
    printf("1. Calculate Absorbance (A)\n");
    printf("2. Calculate Molar Absorptivity (epsilon)\n");
    printf("3. Calculate Path Length (b)\n");
    printf("4. Calculate Concentration (c)\n");
    printf("Enter the variable you want to calculate (1-4): ");
     if (scanf("%d", &unknown) != 1 || unknown < 1 || unknown > 4) {
        printf("Invalid choice.\n");
        while (getchar() != '\n');
        return;
    }


    // Prompt for known variables based on the unknown
    if (unknown != 1) {
        printf("Enter Absorbance (A): ");
        // Absorbance can technically be negative in some instruments, but typically non-negative
        if (scanf("%lf", &absorbance) != 1) { printf("Invalid input.\n"); while(getchar() != '\n'); return; }
    }
    if (unknown != 2) {
        printf("Enter Molar Absorptivity (epsilon): ");
         if (scanf("%lf", &molar_absorptivity) != 1 || molar_absorptivity < 0) { printf("Invalid input. Must be non-negative.\n"); while(getchar() != '\n'); return; }
    }
    if (unknown != 3) {
        printf("Enter Path Length (b): ");
         if (scanf("%lf", &path_length) != 1 || path_length <= 0) { printf("Invalid input. Must be positive.\n"); while(getchar() != '\n'); return; }
    }
    if (unknown != 4) {
        printf("Enter Concentration (c): ");
         if (scanf("%lf", &concentration) != 1 || concentration < 0) { printf("Invalid input. Must be non-negative.\n"); while(getchar() != '\n'); return; }
    }


    // Perform calculation and check for division by zero
    double result;
    switch (unknown) {
        case 1: // Calculate A = epsilon * b * c
            result = molar_absorptivity * path_length * concentration;
            printf("Absorbance A = %.4f\n", result);
            break;
        case 2: // Calculate epsilon = A / (b * c)
            if (path_length == 0 || concentration == 0) {
                printf("Error: Path length and concentration cannot be zero for this calculation.\n");
                return;
            }
            result = absorbance / (path_length * concentration);
            printf("Molar Absorptivity epsilon = %.4f\n", result);
            break;
        case 3: // Calculate b = A / (epsilon * c)
             if (molar_absorptivity == 0 || concentration == 0) {
                printf("Error: Molar absorptivity and concentration cannot be zero for this calculation.\n");
                return;
            }
            result = absorbance / (molar_absorptivity * concentration);
            printf("Path Length b = %.4f\n", result);
            break;
        case 4: // Calculate c = A / (epsilon * b)
             if (molar_absorptivity == 0 || path_length == 0) {
                printf("Error: Molar absorptivity and path length cannot be zero for this calculation.\n");
                return;
            }
            result = absorbance / (molar_absorptivity * path_length);
            printf("Concentration c = %.4f\n", result);
            break;
    }
}

// --- pH and pOH Calculations ---
// Convert between [H+], [OH-], pH, pOH
void calculatePH() {
    int calc_type;
    double value = 0.0;

    printf("\n--- pH and pOH Calculations ---\n");
    printf("Select the calculation you want to perform:\n");
    printf("1. Calculate pH from [H+]\n");
    printf("2. Calculate [H+] from pH\n");
    printf("3. Calculate pOH from [OH-]\n");
    printf("4. Calculate [OH-] from pOH\n");
    printf("5. Calculate pOH and [OH-] from pH\n");
    printf("6. Calculate pH and [H+] from pOH\n");
    printf("7. Calculate pOH and [OH-] from [H+]\n");
    printf("8. Calculate pH and [H+] from [OH-]\n");
    printf("Enter choice (1-8): ");
     if (scanf("%d", &calc_type) != 1 || calc_type < 1 || calc_type > 8) {
        printf("Invalid choice.\n");
        while (getchar() != '\n');
        return;
    }


    printf("Enter the known value: ");
     if (scanf("%lf", &value) != 1) {
        printf("Invalid input.\n");
        while (getchar() != '\n');
        return;
    }

    // Input validation for concentrations (must be positive)
    if ((calc_type == 1 || calc_type == 7) && value <= 0) {
        printf("Error: [H+] must be positive.\n");
        return;
    }
    if ((calc_type == 3 || calc_type == 8) && value <= 0) {
        printf("Error: [OH-] must be positive.\n");
        return;
    }


    double pH, pOH, h_plus, oh_minus;

    switch (calc_type) {
        case 1: // From [H+] to pH
            pH = -log10(value);
            printf("If [H+] = %.4e M, then pH = %.4f\n", value, pH);
            break;
        case 2: // From pH to [H+]
            h_plus = pow(10, -value);
            printf("If pH = %.4f, then [H+] = %.4e M\n", value, h_plus);
            break;
        case 3: // From [OH-] to pOH
            pOH = -log10(value);
            printf("If [OH-] = %.4e M, then pOH = %.4f\n", value, pOH);
            break;
        case 4: // From pOH to [OH-]
            oh_minus = pow(10, -value);
            printf("If pOH = %.4f, then [OH-] = %.4e M\n", value, oh_minus);
            break;
        case 5: // From pH to pOH and [OH-]
             if (value < -2 || value > 16) { // Warning for extreme pH values
                 printf("Warning: pH value is outside the typical 0-14 range, but calculation will proceed.\n");
            }
            pOH = 14.0 - value; // At 25°C
            oh_minus = pow(10, -pOH);
            printf("If pH = %.4f, then pOH = %.4f and [OH-] = %.4e M (at 25°C)\n", value, pOH, oh_minus);
            break;
        case 6: // From pOH to pH and [H+]
             if (value < -2 || value > 16) { // Warning for extreme pOH values
                 printf("Warning: pOH value is outside the typical 0-14 range, but calculation will proceed.\n");
            }
            pH = 14.0 - value; // At 25°C
            h_plus = pow(10, -pH);
            printf("If pOH = %.4f, then pH = %.4f and [H+] = %.4e M (at 25°C)\n", value, pH, h_plus);
            break;
        case 7: // From [H+] to pOH and [OH-]
             if (value <= 0) { printf("Error: [H+] must be positive.\n"); return; }
            pH = -log10(value);
            pOH = 14.0 - pH; // At 25°C
            oh_minus = pow(10, -pOH);
             printf("If [H+] = %.4e M, then pH = %.4f, pOH = %.4f, and [OH-] = %.4e M (at 25°C)\n", value, pH, pOH, oh_minus);
            break;
        case 8: // From [OH-] to pH and [H+]
            if (value <= 0) { printf("Error: [OH-] must be positive.\n"); return; }
            pOH = -log10(value);
            pH = 14.0 - pOH; // At 25°C
            h_plus = pow(10, -pH);
             printf("If [OH-] = %.4e M, then pOH = %.4f, pH = %.4f, and [H+] = %.4e M (at 25°C)\n", value, pOH, pH, h_plus);
            break;
    }
     printf("Note: pH + pOH = 14.0 calculation assumes a temperature of 25°C.\n");
}

// --- Henderson-Hasselbalch Equation pH = pKa + log([A-]/[HA]) ---
// Calculate pH, pKa, or the ratio [A-]/[HA]
void calculateHendersonHasselbalch() {
    int unknown;
    double pH = 0.0, pKa = 0.0, ratio_A_HA = 0.0;

    printf("\n--- Henderson-Hasselbalch Equation Calculation ---\n");
    printf("pH = pKa + log([A-]/[HA])\n");
    printf("Calculate one variable given the other two.\n");
    printf("1. Calculate pH\n");
    printf("2. Calculate pKa\n");
    printf("3. Calculate the ratio [A-]/[HA]\n");
    printf("Enter the variable you want to calculate (1-3): ");
     if (scanf("%d", &unknown) != 1 || unknown < 1 || unknown > 3) {
        printf("Invalid choice.\n");
        while (getchar() != '\n');
        return;
    }


    // Prompt for known variables based on the unknown
    if (unknown != 1) {
        printf("Enter pH: ");
         if (scanf("%lf", &pH) != 1) { printf("Invalid input.\n"); while(getchar() != '\n'); return; }
    }
    if (unknown != 2) {
        printf("Enter pKa: ");
         if (scanf("%lf", &pKa) != 1) { printf("Invalid input.\n"); while(getchar() != '\n'); return; }
    }
    if (unknown != 3) {
        printf("Enter the ratio [A-]/[HA]: ");
         if (scanf("%lf", &ratio_A_HA) != 1 || ratio_A_HA <= 0) { printf("Invalid input. Ratio must be positive.\n"); while(getchar() != '\n'); return; }
    }

    double result;
    switch (unknown) {
        case 1: // Calculate pH = pKa + log10([A-]/[HA])
             if (ratio_A_HA <= 0) {
                 printf("Error: The ratio [A-]/[HA] must be positive to calculate the logarithm.\n");
                 return;
             }
            result = pKa + log10(ratio_A_HA);
            printf("Given pKa = %.4f and [A-]/[HA] = %.4f, pH = %.4f\n", pKa, ratio_A_HA, result);
            break;
        case 2: // Calculate pKa = pH - log10([A-]/[HA])
             if (ratio_A_HA <= 0) {
                 printf("Error: The ratio [A-]/[HA] must be positive to calculate the logarithm.\n");
                 return;
             }
            result = pH - log10(ratio_A_HA);
            printf("Given pH = %.4f and [A-]/[HA] = %.4f, pKa = %.4f\n", pH, ratio_A_HA, result);
            break;
        case 3: // Calculate [A-]/[HA] = 10^(pH - pKa)
            result = pow(10, pH - pKa);
            printf("Given pH = %.4f and pKa = %.4f, the ratio [A-]/[HA] = %.4f\n", pH, pKa, result);
            break;
    }
}

// --- Radioactive Decay Calculation ---
// Calculate remaining activity or time elapsed
void calculateRadioactiveDecay() {
    int calc_type;
    double initial_amount = 0.0, final_amount = 0.0, half_life = 0.0, time_elapsed = 0.0;
    double lambda; // Decay constant

    printf("\n--- Radioactive Decay Calculation ---\n");
    printf("Select the calculation you want to perform:\n");
    printf("1. Calculate remaining amount/activity\n");
    printf("2. Calculate time elapsed\n");
    printf("3. Calculate half-life\n");
    printf("Enter choice (1-3): ");
     if (scanf("%d", &calc_type) != 1 || calc_type < 1 || calc_type > 3) {
        printf("Invalid choice.\n");
        while (getchar() != '\n');
        return;
    }


    // A = A0 * exp(-lambda * t)
    // A0 = initial amount/activity
    // A = final amount/activity
    // lambda = decay constant = ln(2) / half-life
    // t = time elapsed

    // Prompt for known variables based on the unknown
    if (calc_type != 1) {
        printf("Enter Initial Amount/Activity (A0): ");
         if (scanf("%lf", &initial_amount) != 1 || initial_amount <= 0) { printf("Invalid input. Must be positive.\n"); while(getchar() != '\n'); return; }
    }
     if (calc_type != 1 && calc_type != 2) { // Need final amount for calculating time or half-life
         printf("Enter Final Amount/Activity (A): ");
         if (scanf("%lf", &final_amount) != 1 || final_amount < 0) { printf("Invalid input. Must be non-negative.\n"); while(getchar() != '\n'); return; }
         if (final_amount > initial_amount && calc_type != 3) { // Cannot have more activity than starting, unless calculating half-life from 0 -> A
             printf("Warning: Final amount/activity is greater than initial amount/activity. Check inputs.\n");
             // Allow calculation to proceed, might be intended for specific cases.
         }
     }

    if (calc_type != 3) { // Need half-life or decay constant for calculating remaining amount or time
         printf("Enter Half-life (in desired time units, e.g., days): ");
          if (scanf("%lf", &half_life) != 1 || half_life <= 0) { printf("Invalid input. Must be positive.\n"); while(getchar() != '\n'); return; }
         lambda = log(2) / half_life; // Calculate decay constant
    }

    if (calc_type != 2) { // Need time elapsed for calculating remaining amount or half-life
        printf("Enter Time Elapsed (in the same units as half-life): ");
         if (scanf("%lf", &time_elapsed) != 1 || time_elapsed < 0) { printf("Invalid input. Must be non-negative.\n"); while(getchar() != '\n'); return; }
    }


    double result;
    switch (calc_type) {
        case 1: // Calculate remaining amount A = A0 * exp(-lambda * t)
            // Need A0, half-life (-> lambda), t
            printf("Enter Initial Amount/Activity (A0): ");
             if (scanf("%lf", &initial_amount) != 1 || initial_amount <= 0) { printf("Invalid input. Must be positive.\n"); while(getchar() != '\n'); return; }
             printf("Enter Half-life (in desired time units): ");
              if (scanf("%lf", &half_life) != 1 || half_life <= 0) { printf("Invalid input. Must be positive.\n"); while(getchar() != '\n'); return; }
             lambda = log(2) / half_life;
             printf("Enter Time Elapsed (in the same units as half-life): ");
              if (scanf("%lf", &time_elapsed) != 1 || time_elapsed < 0) { printf("Invalid input. Must be non-negative.\n"); while(getchar() != '\n'); return; }

            result = initial_amount * exp(-lambda * time_elapsed);
            printf("Remaining Amount/Activity (A) = %.4f\n", result);
            break;
        case 2: // Calculate time elapsed t = -ln(A / A0) / lambda
            // Need A0, A, half-life (-> lambda)
             if (initial_amount <= 0) { printf("Error: Initial amount A0 must be positive.\n"); return; }
             if (final_amount < 0) { printf("Error: Final amount A cannot be negative.\n"); return; }
             if (final_amount > initial_amount) { printf("Error: Final amount A cannot be greater than Initial amount A0 for positive time elapsed.\n"); return; }
             if (final_amount == 0 && initial_amount > 0) {
                  printf("Warning: Calculating time for complete decay is theoretically infinite. Returning a large value.\n");
                  result = -log(1e-15 / initial_amount) / lambda; // Calculate time to a very small amount
             } else {
                 result = -log(final_amount / initial_amount) / lambda;
             }

            printf("Time Elapsed (t) = %.4f (in half-life units)\n", result);
            break;
        case 3: // Calculate half-life = ln(2) / lambda = (t * ln(2)) / -ln(A / A0)
            // Need A0, A, t
             if (initial_amount <= 0) { printf("Error: Initial amount A0 must be positive.\n"); return; }
             if (final_amount <= 0) { printf("Error: Final amount A must be positive to calculate half-life from decay.\n"); return; }
             if (final_amount >= initial_amount) { printf("Error: Final amount A must be less than Initial amount A0 to calculate half-life from decay.\n"); return; }
             if (time_elapsed <= 0) { printf("Error: Time elapsed must be positive.\n"); return; }

            lambda = -log(final_amount / initial_amount) / time_elapsed;
            if (lambda <= 0) { // Should not happen if A < A0 and t > 0, but check
                 printf("Error: Could not calculate a positive decay constant.\n");
                 return;
            }
            result = log(2) / lambda;
            printf("Half-life = %.4f (in time elapsed units)\n", result);
            break;
    }
}

// --- Calculate Single Amino Acid Isoelectric Point (pI) ---
// Uses pKa values for alpha-carboxyl, alpha-amino, and side chain (if applicable)
void calculateAminoAcidPI() {
    char aa_code[2]; // To read the 1-letter code + null terminator
    char aa_upper;
    double pK1, pK2, pKR = 0.0; // pK1: alpha-carboxyl, pK2: alpha-amino, pKR: side chain
    int has_side_chain = 0; // Flag to indicate if it has a titratable side chain

    printf("\n--- Single Amino Acid Isoelectric Point (pI) Calculation ---\n");
    printf("Enter the 1-letter amino acid code (e.g., A, D, K): ");
    scanf("%s", aa_code);

    if (strlen(aa_code) != 1) {
        printf("Invalid input. Please enter a single letter.\n");
        return;
    }

    aa_upper = toupper(aa_code[0]);

    // Assign standard pKa values
    pK1 = PK1_ALPHA_CARBOXYL;
    pK2 = PK2_ALPHA_AMINO;

    // Assign side chain pKa if applicable
    switch (aa_upper) {
        case 'D': pKR = PKR_D; has_side_chain = 1; break; // Aspartic Acid (-)
        case 'E': pKR = PKR_E; has_side_chain = 1; break; // Glutamic Acid (-)
        case 'C': pKR = PKR_C; has_side_chain = 1; break; // Cysteine (-)
        case 'Y': pKR = PKR_Y; has_side_chain = 1; break; // Tyrosine (-)
        case 'H': pKR = PKR_H; has_side_chain = 1; break; // Histidine (+)
        case 'K': pKR = PKR_K; has_side_chain = 1; break; // Lysine (+)
        case 'R': pKR = PKR_R; has_side_chain = 1; break; // Arginine (+)
        // For S and T, side chains are typically considered neutral or have very high pKa
        // We can include them but note they don't significantly affect pI unless at very high pH
        case 'S': pKR = PKR_S; has_side_chain = 1; break; // Serine (Neutral/Very High)
        case 'T': pKR = PKR_T; has_side_chain = 1; break; // Threonine (Neutral/Very High)
        // Other amino acids (A, G, I, L, M, F, P, W, V, N, Q) have non-titratable side chains in this pI context
        case 'A': case 'G': case 'I': case 'L': case 'M': case 'F': case 'P': case 'W': case 'V': case 'N': case 'Q':
            has_side_chain = 0; // Explicitly set for clarity
            break;
        default:
            printf("Error: '%c' is not a recognized standard amino acid 1-letter code.\n", aa_upper);
            return;
    }

    double pI;

    // Calculate pI based on the pKa values surrounding the zwitterionic form (net charge 0)
    // For neutral amino acids (no titratable side chain): pI = (pK1 + pK2) / 2
    // For acidic amino acids (acidic side chain): pI = (pK1 + pKR) / 2 (The two lowest pKa values)
    // For basic amino acids (basic side chain): pI = (pK2 + pKR) / 2 (The two highest pKa values)

    if (!has_side_chain) {
        // Neutral amino acids: A, G, I, L, M, F, P, W, V, N, Q
         switch (aa_upper) {
            case 'A': case 'G': case 'I': case 'L': case 'M': case 'F': case 'P': case 'W': case 'V': case 'N': case 'Q':
                 pI = (pK1 + pK2) / 2.0;
                 printf("The amino acid %c is neutral.\n", aa_upper);
                 printf("Using pKa values: alpha-carboxyl (%.2f), alpha-amino (%.2f)\n", pK1, pK2);
                 printf("Calculated pI = %.2f\n", pI);
                 break;
             default:
                // This case should theoretically not be reached due to the switch above,
                // but included for robustness.
                printf("Internal error: Could not classify amino acid %c as neutral or charged.\n", aa_upper);
                return;
         }

    } else {
        // Amino acids with charged or polar titratable side chains: D, E, C, Y, H, K, R, S, T
        // Need to sort the three pKa values (pK1, pK2, pKR) to find the two that straddle the pI
        double pKas[3] = {pK1, pK2, pKR};
        // Sort the pKa values (Bubble sort is simple for 3 elements)
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2 - i; j++) {
                if (pKas[j] > pKas[j + 1]) {
                    double temp = pKas[j];
                    pKas[j] = pKas[j + 1];
                    pKas[j + 1] = temp;
                }
            }
        }
        // The pI is the average of the two pKa values surrounding the point where net charge is zero.
        // For acidic AAs (D, E, C, Y), the charges go +1 (low pH) -> 0 (between pK1 and pKR) -> -1 (between pKR and pK2) -> -2 (high pH). pI is average of pK1 and pKR.
        // For basic AAs (H, K, R), the charges go +2 (low pH) -> +1 (between pK1 and pKR) -> 0 (between pKR and pK2) -> -1 (high pH). pI is average of pKR and pK2.
        // For S and T (very high pKR), they behave more like neutral AAs at physiological pH.
        // We need to determine which two pKa values are involved in the transition from +1 to 0 or 0 to -1 net charge.

        printf("The amino acid %c has a titratable side chain.\n", aa_upper);
        printf("Using pKa values: alpha-carboxyl (%.2f), alpha-amino (%.2f), side chain (%.2f)\n", pK1, pK2, pKR);
        printf("Sorted pKa values: %.2f, %.2f, %.2f\n", pKas[0], pKas[1], pKas[2]);


        // Determine which average to use based on which pKa values define the zwitterionic state.
        // The zwitterionic state has a net charge of 0.
        // Consider the charges as pH increases:
        // Very low pH: Net charge determined by +1 (alpha-amino) +1 (basic side chain if any) -1 (alpha-carboxyl) = +1 or +2
        // As pH crosses lowest pKa: Lose one proton (e.g., alpha-carboxyl), net charge decreases by 1.
        // As pH crosses middle pKa: Lose another proton (e.g., side chain or alpha-amino), net charge decreases by 1.
        // As pH crosses highest pKa: Lose the last proton, net charge decreases by 1.

        // Find the pKa values that result in a net charge of 0 when averaged.
        // This occurs when the number of positive charges equals the number of negative charges.
        // For a single amino acid:
        // Possible titratable groups: alpha-carboxyl (-COOH), alpha-amino (-NH3+), and side chain (e.g., -COOH, -NH3+, -SH, imidazole, phenol).
        // At very low pH, charge is +1 (alpha-amino) [+1 if basic side chain]. -1 (alpha-carboxyl) = +1 or +2
        // Transition 1 (lowest pKa): Usually alpha-carboxyl (loss of H+ -> -COO-), charge drops by 1. Net charge becomes 0 or +1.
        // Transition 2 (middle pKa): Either side chain or alpha-amino. Charge drops by 1. Net charge becomes -1 or 0.
        // Transition 3 (highest pKa): The remaining group. Charge drops by 1. Net charge becomes -2 or -1.

        // The pI is the pH where the net charge is 0. This pH lies between the two pKa values that define the transition across the zero net charge state.

        // General rule:
        // Acidic side chain (D, E, C, Y): pK1 (COOH), pKR (side chain COOH/SH/phenol), pK2 (NH3+). Order is usually pK1 < pKR < pK2. Net charges: +1 -> 0 -> -1 -> -2. pI is average of pK1 and pKR.
        // Basic side chain (H, K, R): pK1 (COOH), pK2 (NH3+), pKR (side chain NH3+/imidazole/guanidino). Order is usually pK1 < pK2 < pKR or pK1 < pKR < pK2. Net charges: +2 -> +1 -> 0 -> -1. pI is average of pK2 and pKR.
        // Ambiguous (S, T - high pKR): Behaves like neutral AAs at relevant pH ranges. pI is average of pK1 and pK2.

        // Let's use the sorted pKa values to determine the pI.
        // pKas[0] is the lowest pKa, pKas[1] is the middle, pKas[2] is the highest.

        // If the amino acid is acidic (D, E, C, Y), the pI is between pKas[0] (pK1) and pKas[1] (pKR).
        if (aa_upper == 'D' || aa_upper == 'E' || aa_upper == 'C' || aa_upper == 'Y') {
             pI = (pKas[0] + pKas[1]) / 2.0;
             printf("Calculated pI (average of lowest two pKa values) = %.2f\n", pI);
        }
        // If the amino acid is basic (H, K, R), the pI is between pKas[1] (usually pK2 or pKR) and pKas[2] (usually pKR or pK2).
        else if (aa_upper == 'H' || aa_upper == 'K' || aa_upper == 'R') {
             pI = (pKas[1] + pKas[2]) / 2.0;
             printf("Calculated pI (average of highest two pKa values) = %.2f\n", pI);
        }
        // For S and T, their pKR is very high, so the transition to net charge 0 is between pK1 and pK2.
        else if (aa_upper == 'S' || aa_upper == 'T') {
             pI = (pKas[0] + pKas[1]) / 2.0; // Average of pK1 and pK2 (which are the two lowest)
             printf("Calculated pI (average of alpha-carboxyl and alpha-amino pKa values) = %.2f\n", pI);
        }
        else {
             // This case should not be reached if logic above is correct for 'has_side_chain'
             printf("Could not determine pI for amino acid %c with side chain.\n", aa_upper);
             return;
        }
    }
}

// --- Basic Buffer Preparation (Mass Calculation) ---
// Calculates the mass of weak acid and conjugate base needed for a buffer
// Assumes starting from pure solid weak acid and pure solid conjugate base salt.
void calculateBufferPrep() {
    double target_volume = 0.0, target_concentration = 0.0, desired_pH = 0.0, weak_acid_pKa = 0.0;
    double weak_acid_mw = 0.0, conjugate_base_mw = 0.0;

    printf("\n--- Basic Buffer Preparation (Mass Calculation) ---\n");
    printf("Calculates mass of weak acid and conjugate base needed for a buffer.\n");
    printf("Assumes starting from pure solids.\n");

    printf("Enter Target Buffer Volume (in Liters): ");
     if (scanf("%lf", &target_volume) != 1 || target_volume <= 0) { printf("Invalid input. Must be positive.\n"); while(getchar() != '\n'); return; }

    printf("Enter Target Total Buffer Concentration (Weak Acid + Conjugate Base, in M): ");
     if (scanf("%lf", &target_concentration) != 1 || target_concentration <= 0) { printf("Invalid input. Must be positive.\n"); while(getchar() != '\n'); return; }

    printf("Enter Desired pH of the buffer: ");
     if (scanf("%lf", &desired_pH) != 1) { printf("Invalid input.\n"); while(getchar() != '\n'); return; }

    printf("Enter the pKa of the weak acid: ");
     if (scanf("%lf", &weak_acid_pKa) != 1) { printf("Invalid input.\n"); while(getchar() != '\n'); return; }

    printf("Enter the Molecular Weight of the Weak Acid (g/mol): ");
     if (scanf("%lf", &weak_acid_mw) != 1 || weak_acid_mw <= 0) { printf("Invalid input. Must be positive.\n"); while(getchar() != '\n'); return; }

    printf("Enter the Molecular Weight of the Conjugate Base Salt (g/mol): ");
     if (scanf("%lf", &conjugate_base_mw) != 1 || conjugate_base_mw <= 0) { printf("Invalid input. Must be positive.\n"); while(getchar() != '\n'); return; }


    // Use Henderson-Hasselbalch to find the ratio [A-]/[HA]
    // pH = pKa + log([A-]/[HA])
    // pH - pKa = log([A-]/[HA])
    // [A-]/[HA] = 10^(pH - pKa)

    double ratio_A_HA = pow(10, desired_pH - weak_acid_pKa);

    // Total concentration = [HA] + [A-] = Target Concentration
    // [A-] = ratio_A_HA * [HA]
    // [HA] + ratio_A_HA * [HA] = Target Concentration
    // [HA] * (1 + ratio_A_HA) = Target Concentration
    // [HA] = Target Concentration / (1 + ratio_A_HA)

    double conc_HA = target_concentration / (1.0 + ratio_A_HA);

    // [A-] = Target Concentration - [HA]
    double conc_A = target_concentration - conc_HA;

    // Moles needed = Concentration (M) * Volume (L)
    double moles_HA = conc_HA * target_volume;
    double moles_A = conc_A * target_volume;

    // Mass needed = Moles * Molecular Weight
    double mass_HA = moles_HA * weak_acid_mw;
    double mass_A = moles_A * conjugate_base_mw;

    printf("\n--- Results ---\n");
    printf("Ratio [A-]/[HA] calculated from H-H: %.4f\n", ratio_A_HA);
    printf("Concentration of Weak Acid [HA]: %.4f M\n", conc_HA);
    printf("Concentration of Conjugate Base [A-]: %.4f M\n", conc_A);
    printf("To prepare %.4f L of %.4f M buffer at pH %.4f (pKa = %.4f):\n",
           target_volume, target_concentration, desired_pH, weak_acid_pKa);
    printf("You need to weigh out:\n");
    printf("Weak Acid (%.2f g/mol): %.4f grams\n", weak_acid_mw, mass_HA);
    printf("Conjugate Base Salt (%.2f g/mol): %.4f grams\n", conjugate_base_mw, mass_A);
    printf("\nDissolve these amounts in less than the target volume of water,\n");
    printf("then adjust the volume to %.4f L with water.\n", target_volume);
    printf("Note: You may need to adjust the pH slightly with acid or base.\n");
}

// --- Michaelis-Menten Kinetics (Calculate V) ---
// Calculates reaction velocity (V) given Vmax, Km, and substrate concentration [S]
// V = (Vmax * [S]) / (Km + [S])
void calculateMichaelisMenten() {
    double vmax = 0.0, km = 0.0, substrate_conc = 0.0;
    double velocity;

    printf("\n--- Michaelis-Menten Kinetics Calculation ---\n");
    printf("Calculate reaction velocity (V) using V = (Vmax * [S]) / (Km + [S])\n");

    printf("Enter Vmax (Maximum reaction velocity): ");
     if (scanf("%lf", &vmax) != 1 || vmax < 0) { printf("Invalid input. Must be non-negative.\n"); while(getchar() != '\n'); return; }

    printf("Enter Km (Michaelis Constant - in the same units as [S]): ");
     if (scanf("%lf", &km) != 1 || km <= 0) { printf("Invalid input. Must be positive.\n"); while(getchar() != '\n'); return; }

    printf("Enter Substrate Concentration [S] (in the same units as Km): ");
     if (scanf("%lf", &substrate_conc) != 1 || substrate_conc < 0) { printf("Invalid input. Must be non-negative.\n"); while(getchar() != '\n'); return; }

    // Check for division by zero (Km + [S]) - this will only be zero if Km and [S] are both zero,
    // which is handled by the positive/non-negative checks above.
    if (km + substrate_conc == 0) {
         printf("Error: Km and Substrate Concentration cannot both be zero.\n");
         return;
    }

    velocity = (vmax * substrate_conc) / (km + substrate_conc);

    printf("\n--- Result ---\n");
    printf("Given Vmax = %.4f, Km = %.4f, and [S] = %.4f,\n", vmax, km, substrate_conc);
    printf("The reaction velocity (V) = %.4f\n", velocity);
}