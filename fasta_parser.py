# ============================================
# FASTA File Parser
# Author: Padma Shree
# Phase 1 - Project 3
# Features: Parse, GC content, stats,
#            visualization, NCBI download
# ============================================

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from datetime import datetime
import os

# ── Try importing BioPython ──
try:
    from Bio import SeqIO, Entrez
    BIOPYTHON = True
except ImportError:
    BIOPYTHON = False
    print("BioPython not found. NCBI download disabled.")

# ─────────────────────────────────────────────
# CORE FUNCTIONS
# ─────────────────────────────────────────────

def parse_fasta(filepath):
    """Parse a FASTA file and return list of (header, sequence)."""
    sequences = []
    header    = None
    seq_lines = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    sequences.append((header, "".join(seq_lines)))
                header    = line[1:]
                seq_lines = []
            elif line:
                seq_lines.append(line.upper())
        if header is not None:
            sequences.append((header, "".join(seq_lines)))

    return sequences

def gc_content(sequence):
    """Calculate GC content percentage."""
    if len(sequence) == 0:
        return 0.0
    gc = sequence.count("G") + sequence.count("C")
    return round((gc / len(sequence)) * 100, 2)

def analyze_sequences(sequences):
    """Analyze all sequences and return stats."""
    results = []
    for header, seq in sequences:
        results.append({
            "header"    : header,
            "sequence"  : seq,
            "length"    : len(seq),
            "gc"        : gc_content(seq),
            "a_count"   : seq.count("A"),
            "t_count"   : seq.count("T"),
            "g_count"   : seq.count("G"),
            "c_count"   : seq.count("C"),
        })
    return results

def get_stats(results):
    """Get summary statistics."""
    lengths = [r["length"] for r in results]
    gcs     = [r["gc"] for r in results]
    longest  = max(results, key=lambda x: x["length"])
    shortest = min(results, key=lambda x: x["length"])
    return {
        "total"          : len(results),
        "total_bases"    : sum(lengths),
        "avg_length"     : round(sum(lengths) / len(lengths), 2),
        "avg_gc"         : round(sum(gcs) / len(gcs), 2),
        "longest"        : longest,
        "shortest"       : shortest,
        "max_gc"         : max(results, key=lambda x: x["gc"]),
        "min_gc"         : min(results, key=lambda x: x["gc"]),
    }

# ─────────────────────────────────────────────
# VISUALIZATION
# ─────────────────────────────────────────────

def visualize(results, filename_prefix="fasta"):
    """Generate 4 charts in one figure."""
    labels  = [f"Seq{i+1}" for i in range(len(results))]
    lengths = [r["length"] for r in results]
    gcs     = [r["gc"] for r in results]

    fig = plt.figure(figsize=(16, 10))
    fig.suptitle("FASTA File Analysis Report", fontsize=16, fontweight="bold")
    gs  = gridspec.GridSpec(2, 2, figure=fig)

    # Chart 1 — Sequence lengths
    ax1 = fig.add_subplot(gs[0, 0])
    bars = ax1.bar(labels, lengths, color="#4CAF50", edgecolor="black")
    for bar, val in zip(bars, lengths):
        ax1.text(bar.get_x() + bar.get_width()/2,
                 bar.get_height() + 0.5, str(val),
                 ha="center", fontsize=9, fontweight="bold")
    ax1.set_title("Sequence Lengths", fontweight="bold")
    ax1.set_xlabel("Sequence")
    ax1.set_ylabel("Length (bases)")

    # Chart 2 — GC content
    ax2 = fig.add_subplot(gs[0, 1])
    colors = ["#E91E63" if gc > 60 else "#2196F3" if gc < 40 else "#FF9800" for gc in gcs]
    bars2  = ax2.bar(labels, gcs, color=colors, edgecolor="black")
    for bar, val in zip(bars2, gcs):
        ax2.text(bar.get_x() + bar.get_width()/2,
                 bar.get_height() + 0.5, f"{val}%",
                 ha="center", fontsize=9, fontweight="bold")
    ax2.axhline(y=50, color="gray", linestyle="--", alpha=0.7, label="50% line")
    ax2.set_title("GC Content %", fontweight="bold")
    ax2.set_xlabel("Sequence")
    ax2.set_ylabel("GC Content (%)")
    ax2.legend()

    # Chart 3 — Base composition of first sequence
    ax3  = fig.add_subplot(gs[1, 0])
    r    = results[0]
    bases = ["A", "T", "G", "C"]
    counts = [r["a_count"], r["t_count"], r["g_count"], r["c_count"]]
    ax3.pie(counts, labels=bases, autopct="%1.1f%%",
            colors=["#4CAF50", "#2196F3", "#FF9800", "#E91E63"],
            startangle=90)
    ax3.set_title(f"Base Composition — Seq1", fontweight="bold")

    # Chart 4 — GC content scatter
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.scatter(lengths, gcs, color="#9C27B0", s=100, edgecolors="black", zorder=5)
    for i, label in enumerate(labels):
        ax4.annotate(label, (lengths[i], gcs[i]),
                     textcoords="offset points", xytext=(5, 5), fontsize=8)
    ax4.set_title("Length vs GC Content", fontweight="bold")
    ax4.set_xlabel("Sequence Length")
    ax4.set_ylabel("GC Content (%)")

    plt.tight_layout()
    chart_file = f"{filename_prefix}_analysis.png"
    plt.savefig(chart_file, dpi=150)
    print(f"\n   Chart saved as: {chart_file}")
    plt.show()

# ─────────────────────────────────────────────
# NCBI DOWNLOAD
# ─────────────────────────────────────────────

def download_from_ncbi(accession, email, output_file="ncbi_sequence.fasta"):
    """Download a sequence from NCBI using accession number."""
    if not BIOPYTHON:
        print("BioPython required for NCBI download!")
        return None

    print(f"\n   Downloading {accession} from NCBI...")
    Entrez.email = email
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession,
                               rettype="fasta", retmode="text")
        record = handle.read()
        handle.close()
        with open(output_file, "w") as f:
            f.write(record)
        print(f"   Downloaded and saved as: {output_file}")
        return output_file
    except Exception as e:
        print(f"   Download failed: {e}")
        return None

# ─────────────────────────────────────────────
# SAVE RESULTS
# ─────────────────────────────────────────────

def save_results(results, stats):
    """Save analysis results to a text file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename  = f"fasta_results_{timestamp}.txt"

    with open(filename, "w") as f:
        f.write("FASTA FILE PARSER — RESULTS\n")
        f.write(f"Generated: {datetime.now().strftime('%d-%m-%Y %H:%M:%S')}\n")
        f.write("="*60 + "\n\n")
        f.write("SUMMARY STATISTICS\n")
        f.write("-"*60 + "\n")
        f.write(f"Total sequences    : {stats['total']}\n")
        f.write(f"Total bases        : {stats['total_bases']}\n")
        f.write(f"Average length     : {stats['avg_length']}\n")
        f.write(f"Average GC content : {stats['avg_gc']}%\n")
        f.write(f"Longest sequence   : {stats['longest']['header'][:50]} ({stats['longest']['length']} bases)\n")
        f.write(f"Shortest sequence  : {stats['shortest']['header'][:50]} ({stats['shortest']['length']} bases)\n")
        f.write(f"Highest GC         : {stats['max_gc']['header'][:50]} ({stats['max_gc']['gc']}%)\n")
        f.write(f"Lowest GC          : {stats['min_gc']['header'][:50]} ({stats['min_gc']['gc']}%)\n")
        f.write("\n" + "="*60 + "\n\n")
        f.write("INDIVIDUAL SEQUENCE DETAILS\n")
        f.write("-"*60 + "\n")
        for i, r in enumerate(results):
            f.write(f"\nSequence {i+1}: {r['header'][:60]}\n")
            f.write(f"  Length     : {r['length']} bases\n")
            f.write(f"  GC Content : {r['gc']}%\n")
            f.write(f"  Bases      : A={r['a_count']} T={r['t_count']} G={r['g_count']} C={r['c_count']}\n")
            f.write("-"*40 + "\n")

    print(f"\n   Results saved to: {filename}")

# ─────────────────────────────────────────────
# CREATE SAMPLE FASTA
# ─────────────────────────────────────────────

def create_sample_fasta():
    """Create a sample FASTA file for testing."""
    sample = """>Sequence_1 Homo sapiens BRCA1 partial
ATGCTAGCTAGCTAGCATGCATGCATGCATGCATGCTAGC
TAGCTAGCATGCATGCATGCATGCATGCATGCATGCTAGC
>Sequence_2 Escherichia coli 16S rRNA partial
GGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
>Sequence_3 Mus musculus p53 gene partial
ATGATGATGATGATGATGATGATGATGATGATGATGATGA
TGATGATGATGATGATGATGATGATGATGATGATGATGAT
>Sequence_4 Arabidopsis thaliana chloroplast
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
>Sequence_5 Saccharomyces cerevisiae COX1
ATGTTATTAATTAGAATAAGAATAATAAATAAATAAATAA
ATAAATAAATAAATAAATAAATAAATAAATAAATAAATAA
"""
    with open("sample.fasta", "w") as f:
        f.write(sample)
    print("   Sample FASTA file created: sample.fasta")
    return "sample.fasta"

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

if __name__ == "__main__":
    print("\n" + "="*60)
    print("    FASTA File Parser & Analyzer")
    print("    Author: Padma Shree | Phase 1 - Project 3")
    print("="*60)

    print("\nOptions:")
    print("  1. Parse an existing FASTA file")
    print("  2. Use sample FASTA file (auto-created)")
    print("  3. Download from NCBI")

    choice = input("\nEnter choice (1/2/3): ").strip()

    filepath = None

    if choice == "1":
        filepath = input("Enter path to your FASTA file: ").strip()
        if not os.path.exists(filepath):
            print("File not found!")
            exit()

    elif choice == "2":
        filepath = create_sample_fasta()

    elif choice == "3":
        if not BIOPYTHON:
            print("Please install BioPython: pip install biopython")
            exit()
        email     = input("Enter your email (required by NCBI): ").strip()
        accession = input("Enter NCBI accession number (e.g. NM_007294): ").strip()
        filepath  = download_from_ncbi(accession, email)
        if not filepath:
            exit()
    else:
        print("Invalid choice!")
        exit()

    # Parse and analyze
    print(f"\n   Parsing: {filepath}")
    sequences = parse_fasta(filepath)
    print(f"   Found {len(sequences)} sequences!")

    results = analyze_sequences(sequences)
    stats   = get_stats(results)

    # Display results
    print("\n" + "="*60)
    print("  SUMMARY")
    print("="*60)
    print(f"  Total sequences    : {stats['total']}")
    print(f"  Total bases        : {stats['total_bases']}")
    print(f"  Average length     : {stats['avg_length']} bases")
    print(f"  Average GC content : {stats['avg_gc']}%")
    print(f"  Longest sequence   : {stats['longest']['header'][:40]} ({stats['longest']['length']} bases)")
    print(f"  Shortest sequence  : {stats['shortest']['header'][:40]} ({stats['shortest']['length']} bases)")
    print(f"  Highest GC         : {stats['max_gc']['header'][:40]} ({stats['max_gc']['gc']}%)")
    print(f"  Lowest GC          : {stats['min_gc']['header'][:40]} ({stats['min_gc']['gc']}%)")

    print("\n" + "="*60)
    print("  INDIVIDUAL SEQUENCES")
    print("="*60)
    for i, r in enumerate(results):
        print(f"\n  Sequence {i+1}: {r['header'][:50]}")
        print(f"    Length     : {r['length']} bases")
        print(f"    GC Content : {r['gc']}%")
        print(f"    Bases      : A={r['a_count']} T={r['t_count']} G={r['g_count']} C={r['c_count']}")

    # Visualize and save
    visualize(results)
    save_results(results, stats)

    print("\n  All done! Happy parsing, Paddu! 🧬")
    print("="*60)