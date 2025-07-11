# Peptide Sequencing Algorithms

This repository contains implementations of core algorithms used in cyclo-peptide sequencing based on mass spectrometry (MS) data. The focus is on reconstructing amino acid sequences (peptides) from experimental spectra, using approaches such as branch-and-bound algorithms and leaderboard-based scoring.

---

## 🧬 What is a Peptide?

A peptide is a short chain of amino acids linked by peptide bonds. Peptides are the building blocks of proteins and play critical roles in biological processes such as signaling, immune responses, and metabolism.

In proteomics, scientists often need to determine the exact sequence of a peptide based on data produced by mass spectrometry. This technique breaks peptides into fragments and measures their masses, generating a spectrum that reflects the chemical composition of the original peptide. Since the spectrum alone does not explicitly show the sequence, computational algorithms are used to reconstruct the peptide from the mass data by exploring possible amino acid combinations that explain the observed fragment masses.

---

## 📁 Files in This Repository

- `generate_theoretical_spectrum.py`: Generates the theoretical mass spectrum for a given peptide sequence.
- `branch_and_bound_cyclopeptide_sequencing.py`: Finds all possible candidate peptides for a given mass spectrum.
- `cyclopeptide_sequencing_with_leaderboard_and_convolution.py`: Finds the most likely cadidate peptide by using scoring and mass differences (convolution).
- `linear_spectrum_output.txt`: An example input file of a theoretical mass spectrum for a linear peptide.
- `cyclic_spectrum_output.txt`: An example input file of a theoretical mass spectrum for a cyclic peptide.
- `cyclic_spectrum_with_error_output.txt`: An example input file of a simulated mass spectrum for a cyclic peptide.

---

## ⚙️ How to Use

### 1. Prepare Input

There are two main inputs depending on the program being used:

- Enter a (one-letter code) peptide sequence as a string text, if using `generate_theoretical_spectrum.py`
- Use experimental or theoretical mass spectrum data (with `generate_theoretical_spectrum.py`) as a text file containing values rounded to 1 d.p. (e.g. `cyclic_spectrum_output.txt`), if using `branch_and_bound_cyclopeptide_sequencing.py` or `cyclopeptide_sequencing_with_leaderboard_and_convolution.py`

### 2. Run the Algorithms

Each script reads from the input file and prints:

- Three text files with theoretical **mass spectra** for linear and cyclic peptides (`generate_theoretical_spectrum.py`)
- A list of all candidate peptides, for the `branch_and_bound_cyclopeptide_sequencing.py`
- The peptide sequence and the score of the highest scoring candidate, if using `cyclopeptide_sequencing_with_leaderboard_and_convolution.py`

---

#### Generate Theoretical Mass Spectra

  bash
```generate_theoretical_spectrum.py```

#### Branch-and-bound Sequencing Algorithm

  bash
```branch_and_bound_cyclopeptide_sequencing.py```

#### Cyclopeptide Sequencing with Leaderboard and Convolution

  bash
```cyclopeptide_sequencing_with_leaderboard_and_convolution.py```

The parameters N (leaderboard size), M (allowed masses), T (tolerance threshold), and C (maximum charge) are defined, but can be changed at the bottom when calling the functions.

---

## 🧠 Algorithm Overviews

### Theoretical Spectrum Generation

- Computes the theoretical spectrum of a peptide by calculating all possible subpeptide masses in linear or cyclic form.
- Useful for simulating expected mass spectra from known peptides and introducing noise/errors to mimic real experimental data.
- Time complexity: O(n^2 * log n)

### Branch and Bound Cyclopeptide Sequencing

- Uses a recursive peptide expansion strategy combined with pruning to find peptides consistent with an input experimental spectrum.
- Filters candidate peptides by checking if their linear spectrum matches subsets of the experimental spectrum, reducing search space.
- Guarantees to find all peptides whose cyclic spectrum exactly matches the experimental spectrum.
- Time complexity: O(A^L * L^2 * log L)

### Cyclopeptide Sequencing with Leaderboard and Convolution

- Improves sequencing accuracy by focusing on amino acid masses frequently observed in the spectral convolution of the input spectrum.
- Maintains a leaderboard of top-scoring candidate peptides, expanding and trimming it iteratively based on how well peptides’ spectra match the experimental data.
- Balances exploration and pruning with parameters controlling the leaderboard size and mass tolerance, making it scalable to noisy or large spectra.
- Time complexity: O(n^2 * log n + N * (n * L + L^3 * log L))

---

## 🧪 Example Output

- Theretical Cyclic Spectrum:
  0.0 113.1 114.1 128.1 129.1 227.2 242.2 242.2 257.2 355.3 356.3 370.3 371.3 484.4

- Cyclocpeptide Sequencing:
  Peptide: SPVM | Score: 14

---

## 👤 Author

Heitor Gelain do Nascimento
Email: heitorgelain@outlook.com
GitHub: @heitor-sg5

---

## 📚 References

Bioinformatics Algorithms: An Active Learning Approach (Chapter 4) by
Phillip Compeau & Pavel Pevzner
https://bioinformaticsalgorithms.com
