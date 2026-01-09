# Genomic Sequence Motif Analysis: PWM & ICM Calculator

This project implements a sequence motif analysis pipeline using **R**. It calculates the position-dependent preferences of nucleotides in a set of genomic sequences, generating key matrices used in computational genomics.

### Overview
The script processes a set of DNA sequences to characterize a binding site or motif by calculating:
1. **Position Weight Matrix (PWM):** Frequency of each nucleotide at each position.
2. **Position Probability Matrix (PPM):** Normalized probabilities including **pseudocounts** to handle zero-frequency events.
3. **Log-Likelihood Ratio Matrix:** Log-odds scores comparing the motif model against a null background model.
4. **Information Content Matrix (ICM):** Measures the conservation (in bits) at each position using Shannon Entropy.

### Example Output
The script generates matrices like the Information Content Matrix (ICM), which identifies the strength of the motif:

|  | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| **A** | 1.2743811 | 0.1820544 | 0.0000000 | 0 | 2 | 1.2743811 | 0.8915414 | 0.1173402 |
| **C** | 0.0000000 | 0.0000000 | 0.0000000 | 2 | 0 | 0.0000000 | 0.0000000 | 0.0000000 |
| **G** | 0.0000000 | 0.0000000 | 1.2743811 | 0 | 0 | 0.0000000 | 0.0000000 | 0.1173402 |
| **T** | 0.1820544 | 1.2743811 | 0.1820544 | 0 | 0 | 0.1820544 | 0.2971805 | 0.7040414 |

### 2. Motif Visualization
The following sequence logo shows the information expressed in the ICM matrix, and helps to visualize the consensus sequence:

![Motif Logo](.motif_logo.png)
