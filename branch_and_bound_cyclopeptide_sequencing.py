from collections import Counter

def mass(peptide, amino_acid_mass):
    return round(sum(amino_acid_mass[aa] for aa in peptide), 1)

def parent_mass(spectrum):
    return round(max(spectrum), 1)

def expand(peptides, amino_acid_mass):
    expanded = []
    for peptide in peptides:
        for aa in amino_acid_mass.keys():
            expanded.append(peptide + aa)
    return expanded

def cyclic_spectrum(peptide, amino_acid_mass):
    prefix_mass = [0.0]
    for aa in peptide:
        prefix_mass.append(round(prefix_mass[-1] + amino_acid_mass[aa], 1))
    peptide_mass = prefix_mass[-1]

    spectrum = [0.0]
    n = len(peptide)
    for i in range(n):
        for j in range(i + 1, n + 1):
            sub_mass = round(prefix_mass[j] - prefix_mass[i], 1)
            spectrum.append(sub_mass)
            if i > 0 and j < n:
                spectrum.append(round(peptide_mass - sub_mass, 1))
    return sorted(spectrum)

def consistent(peptide, spectrum, amino_acid_mass):
    def linear_spectrum(peptide, amino_acid_mass):
        prefix_mass = [0.0]
        for aa in peptide:
            prefix_mass.append(round(prefix_mass[-1] + amino_acid_mass[aa], 1))
        linear_spec = [0.0]
        n = len(peptide)
        for i in range(n):
            for j in range(i + 1, n + 1):
                linear_spec.append(round(prefix_mass[j] - prefix_mass[i], 1))
        return linear_spec

    lin_spec = linear_spectrum(peptide, amino_acid_mass)
    spec_counts = Counter(round(m, 1) for m in spectrum)
    lin_spec_counts = Counter(lin_spec)
    for mass_, count_ in lin_spec_counts.items():
        if spec_counts[mass_] < count_:
            return False
    return True

def cyclopeptide_sequencing(spectrum_str, amino_acid_mass):
    spectrum = list(map(float, spectrum_str.strip().split()))
    spectrum = sorted(round(m, 1) for m in spectrum)

    peptides = [""]
    final_peptides = []

    parent_mass_val = parent_mass(spectrum)

    while peptides:
        peptides = expand(peptides, amino_acid_mass)
        peptides_copy = peptides.copy()
        for peptide in peptides_copy:
            peptide_mass_val = mass(peptide, amino_acid_mass)
            if abs(peptide_mass_val - parent_mass_val) < 1e-6:
                if cyclic_spectrum(peptide, amino_acid_mass) == spectrum:
                    final_peptides.append(peptide)
                peptides.remove(peptide)
            elif not consistent(peptide, spectrum, amino_acid_mass):
                peptides.remove(peptide)

    return final_peptides


amino_acid_mass = {
    'G': 57.0,  'A': 71.0,  'S': 87.0,  'P': 97.0,  'V': 99.0,
    'T': 101.1, 'C': 103.0, 'I': 113.1, 'L': 113.1, 'N': 114.1,
    'D': 115.1, 'K': 128.1, 'Q': 128.1, 'E': 129.1, 'M': 131.2,
    'O': 132.2, 'H': 137.1, 'F': 147.2, 'R': 156.2, 'Y': 163.2,
    'U': 168.1, 'W': 186.2
}

with open('cyclic_spectrum_output.txt', 'r') as file:
    spectrum = file.read().strip()

print(cyclopeptide_sequencing(spectrum, amino_acid_mass))
