from collections import Counter

def mass(peptide, amino_acid_mass):
    return sum(amino_acid_mass[aa] for aa in peptide)

def parent_mass(spectrum):
    return max(spectrum)

def expand(peptides, amino_acid_mass):
    expanded = []
    for peptide in peptides:
        for aa in amino_acid_mass.keys():
            expanded.append(peptide + aa)
    return expanded

def cyclic_spectrum(peptide, amino_acid_mass):
    prefix_mass = [0]
    for aa in peptide:
        prefix_mass.append(prefix_mass[-1] + amino_acid_mass[aa])
    peptide_mass = prefix_mass[-1]

    spectrum = [0]
    n = len(peptide)
    for i in range(n):
        for j in range(i + 1, n + 1):
            sub_mass = prefix_mass[j] - prefix_mass[i]
            spectrum.append(sub_mass)
            if i > 0 and j < n:
                spectrum.append(peptide_mass - sub_mass)
    return sorted(spectrum)

def linear_spectrum(peptide, amino_acid_mass):
    prefix_mass = [0]
    for aa in peptide:
        prefix_mass.append(prefix_mass[-1] + amino_acid_mass[aa])
    spectrum = [0]
    n = len(peptide)
    for i in range(n):
        for j in range(i + 1, n + 1):
            spectrum.append(prefix_mass[j] - prefix_mass[i])
    return sorted(spectrum)

def score(peptide, spectrum, amino_acid_mass, T, cyclic=True):
    peptide_spectrum = cyclic_spectrum(peptide, amino_acid_mass) if cyclic else linear_spectrum(peptide, amino_acid_mass)
    peptide_spectrum = sorted(peptide_spectrum)
    spectrum = sorted(spectrum)

    score, i, j = 0, 0, 0
    while i < len(peptide_spectrum) and j < len(spectrum):
        diff = peptide_spectrum[i] - spectrum[j]
        if abs(diff) <= T:
            score += 1
            i += 1
            j += 1
        elif peptide_spectrum[i] < spectrum[j]:
            i += 1
        else:
            j += 1

    return score

def trim(leaderboard, spectrum, N, amino_acid_mass, T):
    scored = [(peptide, score(peptide, spectrum, amino_acid_mass, T, cyclic=False)) for peptide in leaderboard]
    scored.sort(key=lambda x: x[1], reverse=True)

    if len(scored) <= N:
        return [p for p, s in scored]

    cutoff_score = scored[N - 1][1]
    trimmed = [p for p, s in scored if s >= cutoff_score]
    return trimmed

def spectral_convolution(spectrum, M):
    convolution = []
    for i in range(len(spectrum)):
        for j in range(i+1, len(spectrum)):
            diff = spectrum[j] - spectrum[i]
            if 57 <= diff <= 200:
                convolution.append(round(diff, 1))
    counts = Counter(convolution)
    if not counts:
        return []

    most_common = counts.most_common()
    result = []
    last_count = None
    for mass, count in most_common:
        if len(result) < M:
            result.append(mass)
            last_count = count
        elif count == last_count:
            result.append(mass)
        else:
            break
    return sorted(result)

def filter_amino_acids_by_mass(amino_acid_mass, allowed_masses, T):
    filtered = {}
    for aa, mass in amino_acid_mass.items():
        if any(abs(mass - m) <= T for m in allowed_masses):
            filtered[aa] = mass
    return filtered

def leaderboard_cyclopeptide_sequencing_with_convolution(spectrum_str, N, M, amino_acid_mass, T, C):
    spectrum = list(map(float, spectrum_str.strip().split()))
    spectrum = [s - C for s in spectrum]
    spectrum = sorted(spectrum)
    parent = parent_mass(spectrum)

    allowed_masses = spectral_convolution(spectrum, M)
    filtered_amino_acid_mass = filter_amino_acids_by_mass(amino_acid_mass, allowed_masses, T)
    if not filtered_amino_acid_mass:
        filtered_amino_acid_mass = amino_acid_mass

    leaderboard = [""]
    leader_peptide = ""

    while leaderboard:
        leaderboard = expand(leaderboard, filtered_amino_acid_mass)
        leaderboard_copy = leaderboard.copy()
        for peptide in leaderboard_copy:
            peptide_mass_val = mass(peptide, filtered_amino_acid_mass)
            if abs(peptide_mass_val - parent) <= T:
                if score(peptide, spectrum, filtered_amino_acid_mass, T=T) > score(leader_peptide, spectrum, filtered_amino_acid_mass, T=T):
                    leader_peptide = peptide
            elif peptide_mass_val > parent + T:
                leaderboard.remove(peptide)
        leaderboard = trim(leaderboard, spectrum, N, filtered_amino_acid_mass, T)

    return leader_peptide

amino_acid_mass = {
    'G': 57.0,  'A': 71.0,  'S': 87.0,  'P': 97.0,  'V': 99.0,
    'T': 101.1, 'C': 103.0, 'I': 113.1, 'L': 113.1, 'N': 114.1,
    'D': 115.1, 'K': 128.1, 'Q': 128.1, 'E': 129.1, 'M': 131.2,
    'O': 132.2, 'H': 137.1, 'F': 147.2, 'R': 156.2, 'Y': 163.2,
    'U': 168.1, 'W': 186.2
}

with open('example_spectrum.txt', 'r') as file:
    spectrum = file.read().strip()
N = 25
M = 20
T = 0.5
C = 1.0
print(leaderboard_cyclopeptide_sequencing_with_convolution(spectrum, N, M, amino_acid_mass, T, C))
