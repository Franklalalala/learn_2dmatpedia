import json
import os
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase import Atom, Atoms
from ase.data import chemical_symbols
from ase.db import connect


def analyze_database(db_file, output_dir):
    """Analyze the contents of the ASE database."""
    db = connect(db_file)

    db_name = os.path.basename(db_file).split('.')[0]
    output_file = os.path.join(output_dir, f"{db_name}_analysis.txt")

    with open(output_file, 'w') as f:
        f.write(f"Database contains {len(db)} materials\n")

        # Get some statistics on material properties
        # formation_energies = []
        band_gaps = []
        formulas = []

        for row in db.select():
            # formation_energies.append(row.data.decomposition_energy)
            band_gaps.append(row.data.bandgap)
            formulas.append(row.formula)

        # f.write("\nFormation energy statistics:\n")
        # f.write(f"  Min: {min(formation_energies):.4f} eV/atom\n")
        # f.write(f"  Max: {max(formation_energies):.4f} eV/atom\n")
        # f.write(f"  Mean: {np.mean(formation_energies):.4f} eV/atom\n")

        f.write("\nBand gap statistics:\n")
        f.write(f"  Min: {min(band_gaps):.4f} eV\n")
        f.write(f"  Max: {max(band_gaps):.4f} eV\n")
        f.write(f"  Mean: {np.mean(band_gaps):.4f} eV\n")

        f.write("\nMaterials in the database:\n")
        for i, formula in enumerate(formulas):
            f.write(f"  {i + 1}. {formula}\n")

    print(f"Analysis saved to {output_file}")


def element_analysis(db_file, output_dir):
    """
    Analyze element occurrence in the database structures.
    For each element in the periodic table, count how many structures contain it.
    """
    db = connect(db_file)

    db_name = os.path.basename(db_file).split('.')[0]
    output_file = os.path.join(output_dir, f"{db_name}_element_analysis.txt")

    # Dictionary to store counts of structures containing each element
    element_counts = {symbol: 0 for symbol in chemical_symbols[1:]}  # Skip the first element (X)

    # Iterate through all structures
    for row in db.select():
        atoms = row.toatoms()
        # Get unique elements in this structure
        unique_elements = set(atoms.get_chemical_symbols())
        # Increment count for each element
        for element in unique_elements:
            if element in element_counts:
                element_counts[element] += 1

    # Filter out elements with zero counts
    element_counts = {k: v for k, v in element_counts.items() if v > 0}

    # Sort by count (descending)
    sorted_elements = sorted(element_counts.items(), key=lambda x: x[1], reverse=True)

    with open(output_file, 'w') as f:
        f.write("Element occurrence in structures:\n")
        for element, count in sorted_elements:
            f.write(f"  {element}: {count} structures\n")

    print(f"Element analysis saved to {output_file}")
    return element_counts


def structure_statistics(db_file, output_dir):
    """
    Generate statistics about the number of atoms and elements in each structure.
    Create histograms with thinner boxes and ensure integer values.
    """
    db = connect(db_file)

    db_name = os.path.basename(db_file).split('.')[0]
    stats_file = os.path.join(output_dir, f"{db_name}_structure_stats.txt")
    plot_file = os.path.join(output_dir, f"{db_name}_structure_statistics.png")

    # Lists to store number of atoms and elements per structure
    atoms_per_structure = []
    elements_per_structure = []

    # Iterate through all structures
    for row in db.select():
        atoms = row.toatoms()
        # Count atoms in this structure - this is already an integer
        num_atoms = len(atoms)
        atoms_per_structure.append(num_atoms)

        # Count unique elements in this structure - this is already an integer
        num_elements = len(set(atoms.get_chemical_symbols()))
        elements_per_structure.append(num_elements)

    # Create histograms
    plt.figure(figsize=(12, 5))  # Wider figure for better spacing

    # Histogram for number of atoms - using thinner boxes
    plt.subplot(1, 2, 1)
    max_atoms = max(atoms_per_structure)
    # Create more bins for thinner boxes - bins of size 5 or 2 depending on range
    bin_size = 2 if max_atoms < 100 else 5
    bins = np.arange(0, max_atoms + bin_size + 1, bin_size)
    plt.hist(atoms_per_structure, bins=bins, color='royalblue', edgecolor='black', linewidth=0.5)
    plt.yscale('log')
    plt.xlabel('Number of atoms')
    plt.ylabel('Number of structures')
    plt.title('Distribution of atoms per structure')

    # Set x-ticks to ensure integer values
    max_xtick = (max_atoms // 10 + 1) * 10  # Round up to nearest 10
    plt.xticks(np.arange(0, max_xtick + 1, 10))

    plt.tight_layout()

    # Histogram for number of elements - using thinner boxes (width=1 for integers)
    plt.subplot(1, 2, 2)
    max_elements = max(elements_per_structure)
    # For elements, always use bin width of 1 for integer counts
    bins = np.arange(0.5, max_elements + 1.5, 1)  # +0.5 offset for proper integer binning
    plt.hist(elements_per_structure, bins=bins, color='royalblue', edgecolor='black', linewidth=0.5)
    plt.yscale('log')
    plt.xlabel('Number of elements')
    plt.ylabel('Number of structures')
    plt.title('Distribution of elements per structure')

    # Force x-ticks to be integers
    plt.xticks(np.arange(1, max_elements + 1, 1))

    plt.tight_layout()

    # Save figure with higher resolution
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()  # Close the figure to free memory

    # Write statistics to file
    with open(stats_file, 'w') as f:
        f.write("Structure statistics:\n")
        f.write(f"  Average atoms per structure: {np.mean(atoms_per_structure):.2f}\n")
        f.write(f"  Median atoms per structure: {np.median(atoms_per_structure):.2f}\n")
        f.write(f"  Min atoms: {min(atoms_per_structure)}, Max atoms: {max(atoms_per_structure)}\n")
        f.write(f"  Average elements per structure: {np.mean(elements_per_structure):.2f}\n")
        f.write(f"  Median elements per structure: {np.median(elements_per_structure):.2f}\n")
        f.write(f"  Min elements: {min(elements_per_structure)}, Max elements: {max(elements_per_structure)}\n")

        # Count occurrences of each number of elements
        element_counts = Counter(elements_per_structure)
        sorted_counts = sorted(element_counts.items())

        f.write("\nBreakdown of number of elements:\n")
        for num_elements, count in sorted_counts:
            f.write(f"  {num_elements} element{'s' if num_elements > 1 else ''}: {count} structures\n")

    print(f"Structure statistics saved to {stats_file}")
    print(f"Structure histograms saved to {plot_file}")

    return atoms_per_structure, elements_per_structure


def analyze_property_distribution(db_file, output_dir, property_name, title, xlabel, ylabel='Number of structures',
                                  log_scale=True):
    """
    Analyze the distribution of a specific property in the database and create a histogram.

    Parameters:
    -----------
    db_file : str
        Path to the ASE database file
    output_dir : str
        Directory to save output files
    property_name : str
        Name of the property to analyze (attribute in row.data)
    title : str
        Title for the histogram
    xlabel : str
        X-axis label for the histogram
    ylabel : str
        Y-axis label for the histogram
    log_scale : bool
        Whether to use log scale for y-axis
    """
    db = connect(db_file)

    db_name = os.path.basename(db_file).split('.')[0]
    stats_file = os.path.join(output_dir, f"{db_name}_{property_name}_stats.txt")
    plot_file = os.path.join(output_dir, f"{db_name}_{property_name}_distribution.png")

    # Collect property values
    property_values = []

    # Iterate through all structures
    total_structures = 0
    structures_with_property = 0

    for row in db.select():
        total_structures += 1
        try:
            value = getattr(row.data, property_name)
            property_values.append(value)
            structures_with_property += 1
        except (AttributeError, KeyError):
            # Skip structures without this property
            continue

    if not property_values:
        print(f"No structures found with property '{property_name}'")
        return None

    # Calculate statistics
    min_val = min(property_values)
    max_val = max(property_values)
    mean_val = np.mean(property_values)
    median_val = np.median(property_values)
    std_val = np.std(property_values)

    # Create histogram
    plt.figure(figsize=(10, 6))

    # Determine appropriate bin width based on the range and number of data points
    range_val = max_val - min_val
    bin_count = min(100, max(10, int(np.sqrt(len(property_values)))))

    plt.hist(property_values, bins=bin_count, color='royalblue', edgecolor='black', linewidth=0.5, alpha=0.7)

    if log_scale:
        plt.yscale('log')

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    # Add vertical lines for mean and median
    plt.axvline(mean_val, color='red', linestyle='dashed', linewidth=1, label=f'Mean: {mean_val:.4f}')
    plt.axvline(median_val, color='green', linestyle='dashed', linewidth=1, label=f'Median: {median_val:.4f}')

    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # Save figure with higher resolution
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()

    # Write statistics to file
    with open(stats_file, 'w') as f:
        f.write(f"{property_name} statistics:\n")
        f.write(f"  Total structures in database: {total_structures}\n")
        f.write(
            f"  Structures with {property_name}: {structures_with_property} ({structures_with_property / total_structures * 100:.2f}%)\n")
        f.write(f"  Min: {min_val:.6f}\n")
        f.write(f"  Max: {max_val:.6f}\n")
        f.write(f"  Range: {range_val:.6f}\n")
        f.write(f"  Mean: {mean_val:.6f}\n")
        f.write(f"  Median: {median_val:.6f}\n")
        f.write(f"  Standard deviation: {std_val:.6f}\n")

        # Add histogram intervals
        f.write("\nHistogram distribution:\n")
        hist, bin_edges = np.histogram(property_values, bins=10)
        for i in range(len(hist)):
            f.write(f"  {bin_edges[i]:.4f} to {bin_edges[i + 1]:.4f}: {hist[i]} structures\n")

    print(f"{property_name} statistics saved to {stats_file}")
    print(f"{property_name} histogram saved to {plot_file}")

    return property_values


def analyze_added_properties(db_file, output_dir):
    """
    Analyze decomposition_energy, exfoliation_energy_per_atom, and total_magnetization.
    """
    print("Analyzing decomposition energy...")
    analyze_property_distribution(
        db_file,
        output_dir,
        'decomposition_energy',
        'Distribution of Decomposition Energy',
        'Decomposition Energy (eV)'
    )

    print("Analyzing exfoliation energy per atom...")
    analyze_property_distribution(
        db_file,
        output_dir,
        'exfoliation_energy_per_atom',
        'Distribution of Exfoliation Energy per Atom',
        'Exfoliation Energy per Atom (eV/atom)'
    )

    print("Analyzing total magnetization...")
    analyze_property_distribution(
        db_file,
        output_dir,
        'total_magnetization',
        'Distribution of Total Magnetization',
        'Total Magnetization (μB)'
    )

    print("Analyzing band gap...")
    analyze_property_distribution(
        db_file,
        output_dir,
        'bandgap',
        'Distribution of band gap',
        'Band Gap (eV)'
    )


if __name__ == "__main__":
    output_folder = './analysis_results'
    db_file = 'all_2DMat.db'

    os.makedirs(output_folder, exist_ok=True)
    print(f"Analyzing database: {db_file}")
    analyze_database(db_file, output_folder)

    print(f"Performing element analysis")
    element_analysis(db_file, output_folder)

    print(f"Generating structure statistics")
    structure_statistics(db_file, output_folder)

    print(f"Analyzing added properties (decomposition_energy, exfoliation_energy_per_atom, total_magnetization)")
    analyze_added_properties(db_file, output_folder)

    print(f"Complete.")