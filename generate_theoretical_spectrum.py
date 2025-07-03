import random

def cyclic_spectrum(peptide, amino_acid_mass):
    prefix_mass = [0.0]
    for aa in peptide:
        prefix_mass.append(prefix_mass[-1] + amino_acid_mass[aa])
    
    peptide_mass = prefix_mass[-1]
    spectrum = [0.0]
    n = len(peptide)
    for i in range(n):
        for j in range(i + 1, n + 1):
            sub_mass = prefix_mass[j] - prefix_mass[i]
            spectrum.append(round(sub_mass, 1))
            if i > 0 and j < n:
                spectrum.append(round(peptide_mass - sub_mass, 1))

    return sorted(spectrum)

def linear_spectrum(peptide, amino_acid_mass):
    prefix_mass = [0.0]
    for aa in peptide:
        prefix_mass.append(prefix_mass[-1] + amino_acid_mass[aa])
    
    spectrum = [0.0]
    n = len(peptide)
    for i in range(n):
        for j in range(i + 1, n + 1):
            sub_mass = prefix_mass[j] - prefix_mass[i]
            spectrum.append(round(sub_mass, 1))

    return sorted(spectrum)

def cyclic_spectrum_with_error(peptide, amino_acid_mass, P):
    full_spectrum = cyclic_spectrum(peptide, amino_acid_mass)
    n_remove = round(len(full_spectrum) * P / 100)
    if n_remove == 0:
        return full_spectrum
    
    indices_to_remove = set(random.sample(range(len(full_spectrum)), n_remove))
    noisy_spectrum = [mass for i, mass in enumerate(full_spectrum) if i not in indices_to_remove]

    return sorted(noisy_spectrum)

amino_acid_mass = {
    'G': 57.0,  'A': 71.0,  'S': 87.0,  'P': 97.0,  'V': 99.0,
    'T': 101.1, 'C': 103.0, 'I': 113.1, 'L': 113.1, 'N': 114.1,
    'D': 115.1, 'K': 128.1, 'Q': 128.1, 'E': 129.1, 'M': 131.2,
    'O': 132.2, 'H': 137.1, 'F': 147.2, 'R': 156.2, 'Y': 163.2,
    'U': 168.1, 'W': 186.2
}

peptide = "MEVPLSPIGT"
P = 10
with open("linear_spectrum_output.txt", "w") as f:
    f.write(' '.join(map(str, linear_spectrum(peptide, amino_acid_mass))))
with open("cyclic_spectrum_output.txt", "w") as f:
    f.write(' '.join(map(str, cyclic_spectrum(peptide, amino_acid_mass))))
with open("cyclic_spectrum_with_error_output.txt", "w") as f:
    f.write(' '.join(map(str, cyclic_spectrum_with_error(peptide, amino_acid_mass, P))))
