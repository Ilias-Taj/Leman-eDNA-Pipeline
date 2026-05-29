"""
Shared plotting and analysis utilities for eDNA metabarcoding notebooks.
Extracted from analysis notebooks to eliminate code duplication.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import re
import gzip
from pathlib import Path
from matplotlib import cm
from matplotlib.gridspec import GridSpec

# ─── Global Configuration ────────────────────────────────────────────────────

CONF_THRESHOLD = 0.8
RANKS = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

_CONF_NOTES = {
    'filtered': f'SINTAX confidence >= {CONF_THRESHOLD:.2f} applied - lower-confidence calls treated as Unassigned.',
    'unfiltered': 'No confidence filter applied - includes low-confidence assignments. Interpret taxa cautiously.',
    'blast': 'Based on BLAST vs NCBI identity (not SINTAX confidence).',
    'qc': 'Pre-taxonomy QC plot - confidence filter not applicable.',
}


def setup_style():
    """Set consistent matplotlib/seaborn theme for all notebooks."""
    sns.set_theme(style="whitegrid")
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.size': 12,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
    })


# ─── Helper Functions ────────────────────────────────────────────────────────

def clean_sample_names(columns):
    """Clean barcode names: 'Sample_barcode01' -> '01'."""
    return [c.replace('Sample_barcode', '').replace('Sample_', '') for c in columns]


def get_tax_prefix(df):
    """Auto-detect taxonomy column prefix (PR2, Porter, MIDORI2, SILVA, eKOI)."""
    for col in df.columns:
        for pfx in ['PR2', 'Porter', 'MIDORI2', 'MIDORI', 'SILVA', 'eKOI']:
            if col.startswith(f"{pfx}_"):
                return pfx
    return "PR2"


def get_sample_cols(df):
    """Get sample columns (barcode read counts)."""
    return [c for c in df.columns if c.startswith('Sample_') and 'unclassified' not in c]


def add_conf_note(fig=None, kind='filtered'):
    """Add confidence interpretation caption below figure."""
    if fig is None:
        fig = plt.gcf()
    txt = _CONF_NOTES.get(kind, _CONF_NOTES['filtered'])
    fig.text(0.5, -0.02, txt, ha='center', va='top', fontsize=9,
             style='italic', color='#555555', wrap=True)


def clean_taxonomy_names(df, prefix):
    """Remove rank prefixes and trailing numeric suffixes from taxonomy columns."""
    rank_re = re.compile(r'^(kingdom|phylum|class|order|family|genus|species)_', re.IGNORECASE)
    tax_cols = [f'{prefix}_{r}' for r in RANKS]
    for col in tax_cols:
        if col in df.columns:
            df[col] = df[col].apply(
                lambda x: rank_re.sub('', str(x)) if pd.notna(x) else x)
            df[col] = df[col].apply(
                lambda x: '_'.join(str(x).rsplit('_', 1)[:-1])
                if pd.notna(x) and '_' in str(x) and str(x).rsplit('_', 1)[-1].isdigit()
                else x)
    return df


# ─── Data Loading ────────────────────────────────────────────────────────────

def load_marker_data(base_path, marker, db):
    """Load comprehensive taxonomy CSV for a marker/db combo.
    
    Args:
        base_path: Path to dataset output dir (e.g., out/Water_eDNA_18S_COI_14_01_26)
        marker: '18S', 'COI', or 'JEDI'
        db: 'pr2', 'silva', 'porter', 'midori2', or 'ekoi'
    
    Returns:
        DataFrame or None if file doesn't exist
    """
    csv_path = Path(base_path) / f'taxonomy_summary/{marker}/{db}/comprehensive_taxonomy_{marker}.csv'
    if csv_path.exists():
        return pd.read_csv(csv_path)
    return None


def print_data_summary(df, label, prefix):
    """Print a compact summary of loaded dataframe."""
    print(f"\n{'='*60}")
    print(f"  {label} ({prefix} database)  -  {df.shape[0]} OTUs x {df.shape[1]} columns")
    print(f"{'='*60}")
    sample_cols = get_sample_cols(df)
    taxonomy_cols = [c for c in df.columns if c.startswith(f"{prefix}_") or c.startswith("NCBI_")]
    meta_cols = [c for c in df.columns if c not in sample_cols + taxonomy_cols]
    print(f"\nMetadata columns:  {meta_cols}")
    print(f"Sample columns:    {sample_cols}")
    print(f"Taxonomy columns:  {taxonomy_cols}")


# ─── QC Plots ────────────────────────────────────────────────────────────────

def confidence_dashboard(datasets, conf_threshold=None):
    """SINTAX confidence distribution dashboard.
    
    Args:
        datasets: list of (df, prefix, label) tuples, e.g. [(df_18s, 'PR2', '18S (PR2)')]
        conf_threshold: override global threshold
    """
    if conf_threshold is None:
        conf_threshold = CONF_THRESHOLD
    
    n_markers = len(datasets)
    fig = plt.figure(figsize=(22, 5 * (1 + n_markers)))
    gs = GridSpec(1 + n_markers, 6, figure=fig, hspace=0.45, wspace=0.35)
    fig.suptitle('SINTAX Confidence Distribution Dashboard', fontsize=18, fontweight='bold')

    # Row 0: Mean confidence per rank (bar charts)
    cols_per = 6 // n_markers
    for col_start, (df, prefix, label) in enumerate(datasets):
        ax = fig.add_subplot(gs[0, col_start*cols_per:(col_start+1)*cols_per])
        means, counts = [], []
        for r in RANKS:
            col = f'{prefix}_{r}_Conf'
            if col in df.columns:
                vals = pd.to_numeric(df[col], errors='coerce').dropna()
                means.append(vals.mean() if len(vals) > 0 else 0)
                counts.append(len(vals))
            else:
                means.append(0); counts.append(0)
        colors = ['#2ecc71' if m >= conf_threshold else '#e74c3c' for m in means]
        bars = ax.bar(RANKS, means, color=colors, edgecolor='white')
        ax.axhline(y=conf_threshold, color='red', linestyle='--', alpha=0.7,
                   label=f'Threshold ({conf_threshold:.2f})')
        ax.set_ylim(0, 1.05)
        ax.set_ylabel('Mean Confidence')
        ax.set_title(f'{label}: Mean Confidence per Rank', fontweight='bold')
        ax.legend(fontsize=7)
        ax.tick_params(axis='x', rotation=45)
        for bar, m, n in zip(bars, means, counts):
            if m > 0:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                        f'{m:.2f} (n={n})', ha='center', va='bottom', fontsize=7)

    # Subsequent rows: histograms per marker
    for row_idx, (df, prefix, label) in enumerate(datasets):
        for rank_idx, rank in enumerate(RANKS):
            ax = fig.add_subplot(gs[1 + row_idx, rank_idx])
            col = f'{prefix}_{rank}_Conf'
            if col in df.columns:
                vals = pd.to_numeric(df[col], errors='coerce').dropna()
                if len(vals) > 0:
                    ax.hist(vals, bins=20, color='#3498db', edgecolor='white', alpha=0.8)
                    ax.axvline(x=conf_threshold, color='red', linestyle='--', alpha=0.7)
                    pct_above = 100 * (vals >= conf_threshold).sum() / len(vals)
                    ax.set_title(f'{label} {rank}\n{pct_above:.0f}% >= {conf_threshold}', fontsize=9)
            ax.set_xlim(0, 1)
            if rank_idx == 0:
                ax.set_ylabel('Count')

    plt.tight_layout()
    add_conf_note(fig, 'qc')
    plt.show()


def db_performance_dashboard(datasets, conf_threshold=None):
    """Database performance dashboard showing assignment rates.
    
    Args:
        datasets: list of (df, prefix, label) tuples
    """
    if conf_threshold is None:
        conf_threshold = CONF_THRESHOLD

    def db_stats(df, prefix):
        stats = {'rank': [], 'pct_any': [], 'pct_conf': [],
                 'pct_reads_any': [], 'pct_reads_conf': []}
        total_otus = len(df)
        sample_cols = get_sample_cols(df)
        total_reads = df['Total_Abundance'].sum() if 'Total_Abundance' in df.columns else df[sample_cols].sum().sum()
        for rank in RANKS:
            col = f'{prefix}_{rank}'
            conf_col = f'{prefix}_{rank}_Conf'
            if col not in df.columns:
                stats['rank'].append(rank)
                for k in ['pct_any', 'pct_conf', 'pct_reads_any', 'pct_reads_conf']:
                    stats[k].append(0)
                continue
            has_any = df[col].notna() & (df[col] != '') & (df[col] != 'Unassigned')
            has_conf = has_any.copy()
            if conf_col in df.columns:
                has_conf = has_any & (pd.to_numeric(df[conf_col], errors='coerce').fillna(0) >= conf_threshold)
            reads_col = 'Total_Abundance' if 'Total_Abundance' in df.columns else None
            if reads_col and reads_col in df.columns:
                reads_any = df.loc[has_any, reads_col].sum()
                reads_conf = df.loc[has_conf, reads_col].sum()
            else:
                reads_any = df.loc[has_any, sample_cols].sum().sum()
                reads_conf = df.loc[has_conf, sample_cols].sum().sum()
            stats['rank'].append(rank)
            stats['pct_any'].append(100 * has_any.sum() / total_otus if total_otus else 0)
            stats['pct_conf'].append(100 * has_conf.sum() / total_otus if total_otus else 0)
            stats['pct_reads_any'].append(100 * reads_any / total_reads if total_reads else 0)
            stats['pct_reads_conf'].append(100 * reads_conf / total_reads if total_reads else 0)
        return pd.DataFrame(stats)

    n = len(datasets)
    fig, axes = plt.subplots(n, 3, figsize=(20, 5*n))
    if n == 1:
        axes = [axes]
    fig.suptitle('Database Performance Dashboard', fontsize=16, fontweight='bold', y=1.02)

    for i, (df, prefix, label) in enumerate(datasets):
        st = db_stats(df, prefix)
        # Panel 1: OTU assignment rate
        ax = axes[i][0] if n > 1 else axes[0]
        x = np.arange(len(RANKS))
        ax.bar(x - 0.2, st['pct_any'], 0.4, label='Any assignment', color='#3498db', alpha=0.8)
        ax.bar(x + 0.2, st['pct_conf'], 0.4, label=f'Confident (>={conf_threshold})', color='#2ecc71', alpha=0.8)
        ax.set_xticks(x); ax.set_xticklabels(RANKS, rotation=45)
        ax.set_ylabel('% of OTUs'); ax.set_ylim(0, 105)
        ax.set_title(f'{label}: OTU Assignment Rate', fontweight='bold')
        ax.legend(fontsize=8)
        # Panel 2: Read-weighted coverage
        ax = axes[i][1] if n > 1 else axes[1]
        ax.bar(x - 0.2, st['pct_reads_any'], 0.4, label='Any assignment', color='#9b59b6', alpha=0.8)
        ax.bar(x + 0.2, st['pct_reads_conf'], 0.4, label=f'Confident (>={conf_threshold})', color='#e67e22', alpha=0.8)
        ax.set_xticks(x); ax.set_xticklabels(RANKS, rotation=45)
        ax.set_ylabel('% of Reads'); ax.set_ylim(0, 105)
        ax.set_title(f'{label}: Read-Weighted Coverage', fontweight='bold')
        ax.legend(fontsize=8)
        # Panel 3: Resolution depth
        ax = axes[i][2] if n > 1 else axes[2]
        ax.plot(RANKS, st['pct_conf'], 'o-', color='#e74c3c', linewidth=2, markersize=8)
        ax.fill_between(RANKS, st['pct_conf'], alpha=0.2, color='#e74c3c')
        ax.set_ylabel('% Confidently Resolved'); ax.set_ylim(0, 105)
        ax.set_title(f'{label}: Resolution Depth', fontweight='bold')
        ax.tick_params(axis='x', rotation=45)

    plt.tight_layout()
    add_conf_note(fig, 'qc')
    plt.show()


# ─── Sequence QC Plots ───────────────────────────────────────────────────────

def plot_read_lengths(base_path, markers, marker_colors=None):
    """Plot raw read length distributions per marker from filtered FASTQ files."""
    if marker_colors is None:
        marker_colors = {"18S": "#2196F3", "COI": "#4CAF50", "JEDI": "#FF9800"}
    
    base = Path(base_path)
    barcode_dirs = sorted(base.glob("barcode*"))
    marker_lengths = {m: [] for m in markers}

    for bd in barcode_dirs:
        for marker in markers:
            fq = bd / f"filtered_reads_{marker}.fastq.gz"
            if fq.exists():
                with gzip.open(str(fq), 'rt') as f:
                    for i, line in enumerate(f):
                        if i % 4 == 1:
                            marker_lengths[marker].append(len(line.strip()))

    fig, axes = plt.subplots(1, len(markers), figsize=(14, 5))
    if len(markers) == 1:
        axes = [axes]
    for i, marker in enumerate(markers):
        ax = axes[i]
        lengths = marker_lengths[marker]
        if lengths:
            ax.hist(lengths, bins=100, color=marker_colors.get(marker, '#333'),
                    edgecolor='white', alpha=0.85)
            ax.axvline(np.median(lengths), color='red', ls='--', lw=1.5,
                       label=f'median={np.median(lengths):.0f}bp')
            ax.set_title(f'{marker} \u2014 {len(lengths):,} reads',
                         fontsize=13, fontweight='bold')
            ax.set_xlabel('Read length (bp)')
            ax.set_ylabel('Number of reads')
            ax.legend()
            print(f'\u2713 {marker}: {len(lengths):,} reads, median={np.median(lengths):.0f}bp, '
                  f'range={min(lengths)}-{max(lengths)}bp')
        else:
            ax.text(0.5, 0.5, f'{marker}\nNo reads found', ha='center',
                    va='center', transform=ax.transAxes)
    plt.suptitle('Raw Read Length Distributions', fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.show()


def plot_consensus_lengths(base_path, markers, marker_colors=None):
    """Plot consensus OTU sequence length distributions."""
    if marker_colors is None:
        marker_colors = {"18S": "#2196F3", "COI": "#4CAF50", "JEDI": "#FF9800"}

    from Bio import SeqIO
    base = Path(base_path)

    fig, axes = plt.subplots(1, len(markers), figsize=(14, 5))
    if len(markers) == 1:
        axes = [axes]
    for i, marker in enumerate(markers):
        ax = axes[i]
        fasta = base / f"temp_clustering/consensus_{marker}_clean.fasta"
        if fasta.exists():
            lengths = [len(rec.seq) for rec in SeqIO.parse(str(fasta), "fasta")]
            ax.hist(lengths, bins=50, color=marker_colors.get(marker, '#333'),
                    edgecolor='white', alpha=0.85)
            ax.axvline(np.median(lengths), color='red', ls='--', lw=1.5,
                       label=f'median={np.median(lengths):.0f}bp')
            ax.set_title(f'{marker} — {len(lengths)} OTUs',
                         fontsize=13, fontweight='bold')
            ax.set_xlabel('Sequence length (bp)')
            ax.set_ylabel('Number of OTUs')
            ax.legend()
            print(f'✓ {marker}: {len(lengths)} OTUs, median={np.median(lengths):.0f}bp, '
                  f'range={min(lengths)}-{max(lengths)}bp')
        else:
            ax.text(0.5, 0.5, f'{marker}\nFASTA not found', ha='center',
                    va='center', transform=ax.transAxes)
    plt.suptitle('Consensus OTU Sequence Length Distributions',
                 fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    add_conf_note(kind='qc')
    plt.show()


def plot_reads_per_otu(base_path, markers, marker_colors=None):
    """Plot reads-per-OTU distribution from consensus FASTA headers (log-scale histogram)."""
    if marker_colors is None:
        marker_colors = {"18S": "#2196F3", "COI": "#4CAF50", "JEDI": "#FF9800"}

    base = Path(base_path)

    def parse_cluster_sizes(fasta_path):
        sizes = {}
        try:
            with open(fasta_path) as fh:
                for line in fh:
                    if line.startswith('>'):
                        m = re.search(r'centroid=([^;]+);seqs=(\d+)', line)
                        if m:
                            sizes[m.group(1)] = int(m.group(2))
        except FileNotFoundError:
            pass
        return sizes

    fig, axes = plt.subplots(1, len(markers), figsize=(14, 5))
    if len(markers) == 1:
        axes = [axes]
    for ax, marker in zip(axes, markers):
        fasta = base / f"temp_clustering/consensus_{marker}.fasta"
        sizes = parse_cluster_sizes(str(fasta))
        if not sizes:
            ax.text(0.5, 0.5, f'No consensus FASTA found:\n{fasta}',
                    ha='center', va='center', transform=ax.transAxes, color='red')
            continue

        vals = sorted(sizes.values())
        n_total = len(vals)
        n_spoa = sum(1 for v in vals if v >= 3)
        n_single = sum(1 for v in vals if v < 3)

        bins = np.logspace(0, np.ceil(np.log10(max(vals))), 30)
        ax.hist([v for v in vals if v < 3], bins=bins, color='#F44336', alpha=0.85,
                label=f'Longest read (1–2 reads, n={n_single})')
        ax.hist([v for v in vals if v >= 3], bins=bins, color='#4CAF50', alpha=0.85,
                label=f'spoa POA (≥3 reads, n={n_spoa})')
        ax.axvline(3, ls='--', color='black', lw=1.5, label='spoa threshold (3 reads)')
        ax.set_xscale('log')
        ax.set_xlabel('Reads per OTU (log scale)')
        ax.set_ylabel('Number of OTUs')
        ax.set_title(f'{marker} — Reads per OTU\n'
                     f'({n_spoa}/{n_total} OTUs use spoa POA, {100*n_spoa/n_total:.1f}%)')
        ax.legend(fontsize=9)

        median_spoa = int(np.median([v for v in vals if v >= 3])) if n_spoa else 'N/A'
        ax.text(0.98, 0.95,
                f'Median cluster size: {int(np.median(vals))}\nMedian spoa size: {median_spoa}\nMax: {max(vals)}',
                transform=ax.transAxes, ha='right', va='top',
                fontsize=9, bbox=dict(boxstyle='round', fc='white', alpha=0.8))

    plt.suptitle('Consensus Quality Assessment — Reads per OTU', fontweight='bold')
    plt.tight_layout()
    plt.show()

    # Summary table
    print(f"{'Marker':<8} {'Total OTUs':>12} {'spoa (≥3)':>12} {'Longest rd (<3)':>16} "
          f"{'% spoa':>8} {'Median reads':>14} {'Max reads':>10}")
    print("-" * 80)
    for marker in markers:
        fasta = base / f"temp_clustering/consensus_{marker}.fasta"
        sizes = parse_cluster_sizes(str(fasta))
        if sizes:
            vals = list(sizes.values())
            n_spoa = sum(1 for v in vals if v >= 3)
            print(f"{marker:<8} {len(vals):>12} {n_spoa:>12} {len(vals)-n_spoa:>16} "
                  f"{100*n_spoa/len(vals):>7.1f}% {int(np.median(vals)):>14} {max(vals):>10}")


def read_lengths_fastq_gz(path):
    """Read a FASTQ.gz file and return list of sequence lengths."""
    lengths = []
    try:
        with gzip.open(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    lengths.append(len(line.strip()))
    except Exception:
        pass
    return lengths


def plot_primer_trimming(base_path, markers, marker_colors=None):
    """Plot read length before vs after primer trimming."""
    if marker_colors is None:
        marker_colors = {"18S": "#2196F3", "COI": "#4CAF50", "JEDI": "#FF9800"}
    
    base = Path(base_path)

    fig, axes = plt.subplots(1, len(markers), figsize=(7*len(markers), 5))
    if len(markers) == 1:
        axes = [axes]
    
    for i, marker in enumerate(markers):
        ax = axes[i]
        before, after = [], []
        for bd in sorted(base.glob("barcode*")):
            raw = bd / f"filtered_reads_{marker}.fastq.gz"
            trimmed = bd / f"trimmed_{marker}.fastq.gz"
            if raw.exists():
                before.extend(read_lengths_fastq_gz(str(raw)))
            if trimmed.exists():
                after.extend(read_lengths_fastq_gz(str(trimmed)))
        
        if before or after:
            if before:
                ax.hist(before, bins=60, alpha=0.5, color='gray', label=f'Before (n={len(before):,})')
            if after:
                ax.hist(after, bins=60, alpha=0.7,
                        color=marker_colors.get(marker, '#333'),
                        label=f'After (n={len(after):,})')
            ax.set_title(f'{marker}: Before vs After Primer Trimming', fontweight='bold')
            ax.set_xlabel('Read Length (bp)')
            ax.set_ylabel('Count')
            ax.legend()
        else:
            ax.set_title(f'{marker}: No trimming data')
    plt.tight_layout()
    add_conf_note(kind='qc')
    plt.show()


def plot_barcode_reads(base_path, markers, marker_colors=None):
    """Plot read counts per barcode per marker (from filtered FASTQ files)."""
    if marker_colors is None:
        marker_colors = {"18S": "#2196F3", "COI": "#4CAF50", "JEDI": "#FF9800"}
    
    base = Path(base_path)
    barcode_dirs = sorted(base.glob("barcode*"))
    barcode_labels = [bd.name for bd in barcode_dirs]

    fig, axes = plt.subplots(1, len(markers), figsize=(7*len(markers), 5))
    if len(markers) == 1:
        axes = [axes]
    for ax, marker in zip(axes, markers):
        color = marker_colors.get(marker, '#333')
        counts = []
        for bd in barcode_dirs:
            fq = bd / f"filtered_reads_{marker}.fastq.gz"
            counts.append(len(read_lengths_fastq_gz(str(fq))) if fq.exists() else 0)
        x = np.arange(len(barcode_labels))
        ax.bar(x, counts, color=color, alpha=0.8, edgecolor="black", lw=0.4)
        ax.set_xticks(x)
        ax.set_xticklabels([lbl.replace("barcode", "bc") for lbl in barcode_labels],
                           rotation=45, ha="right", fontsize=8)
        ax.set_xlabel("Barcode")
        ax.set_ylabel("Number of reads")
        ax.set_title(f"{marker} — Read Count by Barcode")
        ax.yaxis.grid(True, linestyle="--", alpha=0.5)
        ax.set_axisbelow(True)
    plt.suptitle("Read Count by Barcode", fontweight="bold")
    plt.tight_layout()
    plt.show()


# ─── Biodiversity Plots ─────────────────────────────────────────────────────

def stacked_bar_compare(df, rank, prefix, marker_label, sample_cols=None,
                        top_n=10, conf_threshold=None):
    """Dual stacked bar: left with Unassigned, right without (renormalized).
    
    Args:
        df: DataFrame with taxonomy and sample columns
        rank: Taxonomic rank (Phylum, Class, etc.)
        prefix: DB prefix (PR2, SILVA, etc.)
        marker_label: Display label (e.g., '18S', 'COI')
        sample_cols: Sample columns (auto-detected if None)
        top_n: Number of top taxa to show
        conf_threshold: Confidence threshold (default: CONF_THRESHOLD)
    """
    if conf_threshold is None:
        conf_threshold = CONF_THRESHOLD
    if sample_cols is None:
        sample_cols = get_sample_cols(df)
    
    col = f'{prefix}_{rank}'
    conf_col = f'{prefix}_{rank}_Conf'
    d = df.copy()
    if conf_col in d.columns:
        d.loc[pd.to_numeric(d[conf_col], errors='coerce').fillna(0) < conf_threshold, col] = ''
    d[col] = d[col].replace('', 'Unassigned').fillna('Unassigned')
    grouped = d.groupby(col)[sample_cols].sum()

    unassigned_row = grouped.loc['Unassigned'] if 'Unassigned' in grouped.index else None
    assigned = grouped.drop('Unassigned', errors='ignore').copy()
    assigned['Total'] = assigned.sum(axis=1)
    assigned = assigned.sort_values('Total', ascending=False)
    top_idx = assigned.head(top_n).index
    top_data = assigned.loc[top_idx].drop(columns='Total')
    others_data = assigned.loc[~assigned.index.isin(top_idx)].drop(columns='Total').sum()

    left = top_data.copy()
    left.loc['Others'] = others_data
    if unassigned_row is not None:
        left.loc['Unassigned'] = unassigned_row
    left_pct = left.div(left.sum(axis=0).replace(0, np.nan), axis=1).fillna(0) * 100

    right = top_data.copy()
    right.loc['Others'] = others_data
    right_pct = right.div(right.sum(axis=0).replace(0, np.nan), axis=1).fillna(0) * 100

    def colors_for(data_):
        palette = cm.tab20(np.linspace(0, 1, max(len(data_) - 2, 1)))
        out, pi = [], 0
        for name in data_.index:
            if name == 'Unassigned': out.append('#D3D3D3')
            elif name == 'Others': out.append('#696969')
            else: out.append(palette[pi % len(palette)]); pi += 1
        return out

    if unassigned_row is not None:
        total_reads = grouped.sum().sum()
        un_share = 100 * unassigned_row.sum() / total_reads if total_reads else 0
        print(f'[{marker_label}] {rank}: Unassigned = {un_share:.1f}% of reads')

    fig, axes = plt.subplots(1, 2, figsize=(20, 7), sharey=True)
    for ax, data, subtitle in [(axes[0], left_pct, 'including Unassigned'),
                                (axes[1], right_pct, 'excluding Unassigned (renormalized)')]:
        data = data.copy()
        data.columns = clean_sample_names(data.columns)
        data.T.plot(kind='bar', stacked=True, ax=ax, width=0.85, color=colors_for(data))
        ax.set_title(f'{rank}-Level Composition ({marker_label}, {prefix})\n{subtitle}',
                     fontweight='bold', fontsize=11)
        ax.set_ylabel('Relative Abundance (%)')
        ax.set_xlabel('Sample')
        h, l = ax.get_legend_handles_labels()
        ax.legend(reversed(h), reversed(l), bbox_to_anchor=(1.02, 1), loc='upper left',
                  title=rank, fontsize=8)
        for _i, _s in enumerate(data.columns):
            if data[_s].sum() == 0:
                ax.text(_i, 2, 'ND', ha='center', va='bottom', fontsize=9,
                        color='gray', fontstyle='italic')
        ax.set_ylim(0, 105)
    plt.tight_layout()
    add_conf_note(fig, 'filtered')
    plt.show()


def plot_top_genera_confident(df, prefix, marker_label, sample_cols=None,
                              top_n=20, conf_threshold=None):
    """Top N genera with confidence >= threshold."""
    if conf_threshold is None:
        conf_threshold = CONF_THRESHOLD
    if sample_cols is None:
        sample_cols = get_sample_cols(df)
    
    genus_col = f'{prefix}_Genus'
    conf_col = f'{prefix}_Genus_Conf'
    d = df.copy()
    d[genus_col] = d[genus_col].fillna('Unassigned')

    if conf_col in d.columns:
        mask_conf = pd.to_numeric(d[conf_col], errors='coerce') >= conf_threshold
        df_conf = d[mask_conf].copy()
    else:
        df_conf = d[d[genus_col] != 'Unassigned'].copy()

    genus_agg = df_conf.groupby(genus_col)[sample_cols].sum()
    genus_agg = genus_agg.drop('Unassigned', errors='ignore')
    genus_agg['Total'] = genus_agg.div(d[sample_cols].sum(axis=0), axis=1).mean(axis=1) * 100
    genus_agg = genus_agg.sort_values('Total', ascending=False)
    top = genus_agg.head(top_n)

    genus_conf = {}
    if conf_col in df_conf.columns:
        for genus in top.index:
            vals = pd.to_numeric(df_conf.loc[df_conf[genus_col] == genus, conf_col], errors='coerce').dropna()
            genus_conf[genus] = vals.mean() if len(vals) > 0 else None

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.barplot(x=top['Total'], y=top.index, palette='viridis', ax=ax)
    ax.set_title(f'Top {top_n} Genera \u2014 Confident Only (\u2265 {conf_threshold:.2f}) ({marker_label}, {prefix})',
                 fontweight='bold')
    ax.set_xlabel('Mean Relative Abundance (% of total reads)')
    ax.set_ylabel('Genus')

    for j, genus in enumerate(top.index):
        conf = genus_conf.get(genus)
        if conf is not None:
            ax.text(top.loc[genus, 'Total'] + 0.001, j,
                    f' ({conf:.2f})', va='center', fontsize=9, color='darkred', fontweight='bold')
    plt.tight_layout()
    add_conf_note(kind='filtered')
    plt.show()


def plot_top_genera_all(df, prefix, marker_label, sample_cols=None,
                        top_n=20, conf_threshold=None):
    """Top N genera by abundance (all confidence levels), with color-coded bars."""
    if conf_threshold is None:
        conf_threshold = CONF_THRESHOLD
    if sample_cols is None:
        sample_cols = get_sample_cols(df)
    
    genus_col = f'{prefix}_Genus'
    conf_col = f'{prefix}_Genus_Conf'
    d = df.copy()
    d[genus_col] = d[genus_col].fillna('Unassigned')

    genus_abund = d.groupby(genus_col)[sample_cols].sum()
    genus_abund = genus_abund.drop('Unassigned', errors='ignore')
    genus_abund['Total'] = genus_abund.div(d[sample_cols].sum(axis=0), axis=1).mean(axis=1) * 100
    genus_abund = genus_abund.sort_values('Total', ascending=False)
    top = genus_abund.head(top_n)

    raw_conf = {}
    if conf_col in d.columns:
        for genus in top.index:
            vals = pd.to_numeric(d.loc[d[genus_col] == genus, conf_col], errors='coerce').dropna()
            raw_conf[genus] = vals.mean() if len(vals) > 0 else None

    fig, ax = plt.subplots(figsize=(11, 8))
    colors = []
    for genus in top.index:
        c = raw_conf.get(genus)
        if c is None: colors.append('#999999')
        elif c >= 0.9: colors.append('#1a6832')
        elif c >= conf_threshold: colors.append('#2d8a4e')
        else: colors.append('#e6a817')

    bars = ax.barh(range(len(top)), top['Total'].values[::-1], color=colors[::-1], edgecolor='white')
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top.index[::-1], fontsize=9)
    ax.set_xlabel('Mean Relative Abundance (%)')
    ax.set_title(f'Top {top_n} Genera by Abundance - All Confidences ({marker_label}, {prefix})',
                 fontweight='bold')

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#1a6832', label=f'High conf (>=0.9)'),
        Patch(facecolor='#2d8a4e', label=f'Confident (>={conf_threshold})'),
        Patch(facecolor='#e6a817', label=f'Low conf (<{conf_threshold})'),
        Patch(facecolor='#999999', label='No confidence'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
    plt.tight_layout()
    add_conf_note(kind='unfiltered')
    plt.show()


def plot_eukaryotic_genera(df, prefix, marker_label, sample_cols=None,
                           top_n=20, conf_threshold=None):
    """Top N eukaryotic genera (domain-filtered, confidence-filtered)."""
    if conf_threshold is None:
        conf_threshold = CONF_THRESHOLD
    if sample_cols is None:
        sample_cols = get_sample_cols(df)
    
    _df = df.copy()
    genus_col = f'{prefix}_Genus'
    conf_col = f'{prefix}_Genus_Conf'
    domain_col = f'{prefix}_Domain'

    if domain_col in _df.columns:
        euk_mask = _df[domain_col].astype(str).str.contains('Eukaryota', case=False, na=False)
        _df = _df[euk_mask].copy()
        print(f"Eukaryota filter: {len(_df)}/{len(df)} OTUs ({100*len(_df)/len(df):.1f}%)")
    else:
        print(f"No Domain column found ({domain_col}), using all OTUs")

    _df[genus_col] = _df[genus_col].fillna('Unassigned')
    if conf_col in _df.columns:
        conf_mask = pd.to_numeric(_df[conf_col], errors='coerce') >= conf_threshold
        _df = _df[conf_mask].copy()

    genus_agg = _df.groupby(genus_col)[sample_cols].sum()
    genus_agg = genus_agg.drop('Unassigned', errors='ignore')
    genus_agg['Total'] = genus_agg.div(df[sample_cols].sum(axis=0), axis=1).mean(axis=1) * 100
    genus_agg = genus_agg.sort_values('Total', ascending=False)
    top = genus_agg.head(top_n)

    if len(top) == 0:
        print("No Eukaryota genera found with confidence >= threshold")
        return

    genus_conf = {}
    if conf_col in _df.columns:
        for g in top.index:
            vals = pd.to_numeric(_df.loc[_df[genus_col] == g, conf_col], errors='coerce').dropna()
            genus_conf[g] = vals.mean() if len(vals) > 0 else None

    fig, ax = plt.subplots(figsize=(11, 8))
    colors = []
    for g in top.index:
        c = genus_conf.get(g)
        if c is None: colors.append('#999999')
        elif c >= 0.9: colors.append('#1a6832')
        elif c >= conf_threshold: colors.append('#2d8a4e')
        else: colors.append('#e6a817')

    bars = ax.barh(range(len(top)), top['Total'].values[::-1], color=colors[::-1], edgecolor='white')
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top.index[::-1], fontsize=9)
    ax.set_xlabel('Mean Relative Abundance (%)')
    ax.set_title(f'Top {top_n} Eukaryotic Genera ({marker_label}, {prefix})\nDomain=Eukaryota, Conf>={conf_threshold}',
                 fontweight='bold')

    for j, g in enumerate(list(top.index)[::-1]):
        c = genus_conf.get(g)
        if c is not None:
            ax.text(top.loc[g, 'Total'] + 0.001, j,
                    f' ({c:.2f})', va='center', fontsize=8, color='darkred')
    plt.tight_layout()
    add_conf_note(kind='filtered')
    plt.show()


# ─── BLAST Validation ────────────────────────────────────────────────────────

def parse_blast_file(filepath):
    """Parse custom BLAST output from 6_blast_top_otus.py."""
    data = []
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        start_reading = False
        for line in lines:
            if line.startswith('---'):
                start_reading = True
                continue
            if not start_reading or not line.strip():
                continue
            parts = line.split('|')
            if len(parts) >= 4:
                otu_id = parts[0].strip()
                species = parts[2].strip()
                identity_str = parts[3].strip().replace('%', '')
                try:
                    reads = float(parts[1].strip())
                    identity = float(identity_str) if identity_str and identity_str != '-' else None
                    data.append({'OTU': otu_id, 'Abundance': reads,
                                 'Species': species, 'Identity': identity})
                except ValueError:
                    pass
    except FileNotFoundError:
        pass
    return pd.DataFrame(data) if data else pd.DataFrame()


def plot_blast_hits(blast_file, marker_label, sort_by='Abundance'):
    """Horizontal bar chart of BLAST top hits."""
    df_blast = parse_blast_file(blast_file)
    if df_blast.empty:
        print(f'No BLAST data for {marker_label} ({sort_by})')
        return
    
    df_blast = df_blast.sort_values('Abundance', ascending=True)
    labels = [f"{row['Species']}  [{row['OTU'].split('_')[-1]}]" for _, row in df_blast.iterrows()]
    df_blast['Label'] = labels
    colors = []
    for _, row in df_blast.iterrows():
        if row['Identity'] is None: colors.append('#999999')
        elif row['Identity'] >= 97: colors.append('#2d8a4e')
        elif row['Identity'] >= 90: colors.append('#e6a817')
        else: colors.append('#c0392b')

    fig, ax = plt.subplots(figsize=(11, max(3, len(df_blast) * 0.55)))
    ax.barh(range(len(df_blast)), df_blast['Abundance'].values,
            color=colors, edgecolor='white', height=0.7)
    ax.set_yticks(range(len(df_blast)))
    ax.set_yticklabels(df_blast['Label'].values, fontsize=9)
    ax.set_xlabel('Relative Abundance', fontsize=11)
    ax.set_title(f'Top {marker_label} OTUs by {sort_by} - BLAST vs NCBI',
                 fontweight='bold', fontsize=13, pad=12)
    xmax = df_blast['Abundance'].max()
    for j, (_, row) in enumerate(df_blast.iterrows()):
        if row['Identity'] is not None:
            ax.text(row['Abundance'] + xmax * 0.015, j,
                    f"{row['Identity']:.1f}%", va='center', fontsize=9,
                    color='#2d8a4e' if row['Identity'] >= 97 else '#c0392b')
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#2d8a4e', label='>= 97% (species)'),
        Patch(facecolor='#e6a817', label='>= 90% (genus)'),
        Patch(facecolor='#c0392b', label='< 90%'),
        Patch(facecolor='#999999', label='No hit'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
    plt.tight_layout()
    add_conf_note(kind='blast')
    plt.show()


def blast_vs_sintax_table(blast_file, df_taxonomy, prefix, marker_label):
    """Styled table comparing BLAST hits with SINTAX taxonomy."""
    blast_tbl = parse_blast_file(blast_file)
    if blast_tbl.empty:
        print(f'No BLAST data for {marker_label} table.')
        return
    
    rows = []
    for _, br in blast_tbl.iterrows():
        otu = br['OTU']
        tr = df_taxonomy[df_taxonomy['OTU_ID'] == otu]
        row = {'OTU': otu, 'Abundance': br['Abundance'],
               'BLAST Top Hit': br['Species'], 'Identity (%)': br['Identity']}
        for rank in RANKS:
            nm_col = f'{prefix}_{rank}'
            cf_col = f'{prefix}_{rank}_Conf'
            if not tr.empty and nm_col in tr.columns:
                nm = tr[nm_col].values[0]
                cf = tr[cf_col].values[0] if cf_col in tr.columns else None
                if pd.notna(cf):
                    row[rank] = f'{nm} ({float(cf):.2f})'
                elif pd.notna(nm):
                    row[rank] = str(nm)
                else:
                    row[rank] = '\u2013'
            else:
                row[rank] = '\u2013'
        rows.append(row)
    
    df_tbl = pd.DataFrame(rows)
    
    def color_id(val):
        try:
            v = float(val)
            if v >= 97: return 'background-color:#d4f0dc;color:#1a6e35;font-weight:bold'
            if v >= 90: return 'background-color:#fef3cd;color:#856404;font-weight:bold'
            return 'background-color:#fde8e8;color:#8b1a1a;font-weight:bold'
        except:
            return ''
    
    styled = (df_tbl.style
        .applymap(color_id, subset=['Identity (%)'])
        .format({'Identity (%)': lambda x: f'{x:.1f}%' if pd.notna(x) else 'N/A',
                 'Abundance': lambda x: f'{x:.4f}'})
        .set_caption(f'<b>{marker_label} - BLAST vs SINTAX (by Abundance)</b>')
        .set_table_styles([
            {'selector': 'th', 'props': [('background-color', '#2c3e50'),
                ('color', 'white'), ('font-weight', 'bold'), ('padding', '6px 10px')]},
            {'selector': 'td', 'props': [('padding', '5px 8px'),
                ('border', '1px solid #dee2e6'), ('white-space', 'nowrap')]},
            {'selector': 'caption', 'props': [('caption-side', 'top'),
                ('font-size', '14px'), ('padding', '8px')]}
        ]).hide(axis='index'))
    display(styled)


# ─── Cross-Marker Comparison ────────────────────────────────────────────────

def cross_marker_genus_comparison(datasets, conf_threshold=None):
    """Compare genus detections across markers (Venn-style summary).
    
    Args:
        datasets: list of (df, prefix, marker_label) tuples
    """
    if conf_threshold is None:
        conf_threshold = CONF_THRESHOLD
    
    genera_sets = {}
    for df, prefix, label in datasets:
        genus_col = f'{prefix}_Genus'
        conf_col = f'{prefix}_Genus_Conf'
        d = df.copy()
        if conf_col in d.columns:
            mask = pd.to_numeric(d[conf_col], errors='coerce').fillna(0) >= conf_threshold
            d = d[mask]
        genera = set(d[genus_col].dropna().unique()) - {'Unassigned', ''}
        genera_sets[label] = genera
        print(f"Genera detected by {label} ({prefix}): {len(genera)}")

    if len(genera_sets) == 2:
        labels = list(genera_sets.keys())
        s1, s2 = genera_sets[labels[0]], genera_sets[labels[1]]
        shared = s1 & s2
        only1, only2 = s1 - s2, s2 - s1
        print(f"\nShared: {len(shared)}")
        print(f"Only {labels[0]}: {len(only1)}")
        print(f"Only {labels[1]}: {len(only2)}")
        if shared:
            print(f"\nShared genera: {sorted(shared)[:20]}{'...' if len(shared) > 20 else ''}")


# ─── Comparative (Multi-DB) Plots ───────────────────────────────────────────

def compare_db_assignment_rates(db_datasets, marker_label, conf_threshold=None):
    """Compare assignment rates across multiple databases for the same marker.
    
    Args:
        db_datasets: list of (df, prefix, db_label) tuples for same marker
        marker_label: e.g., '18S', 'COI', 'JEDI'
    """
    if conf_threshold is None:
        conf_threshold = CONF_THRESHOLD
    
    db_colors = {'eKOI': '#4CAF50', 'SILVA': '#2196F3', 'PR2': '#9C27B0',
                 'Porter': '#FF9800', 'MIDORI2': '#F44336'}
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    x = np.arange(len(RANKS))
    width = 0.8 / len(db_datasets)
    
    for i, (df, prefix, db_label) in enumerate(db_datasets):
        sample_cols = get_sample_cols(df)
        total_otus = len(df)
        pct_conf = []
        for rank in RANKS:
            col = f'{prefix}_{rank}'
            conf_col = f'{prefix}_{rank}_Conf'
            if col not in df.columns:
                pct_conf.append(0)
                continue
            has_any = df[col].notna() & (df[col] != '') & (df[col] != 'Unassigned')
            has_conf = has_any.copy()
            if conf_col in df.columns:
                has_conf = has_any & (pd.to_numeric(df[conf_col], errors='coerce').fillna(0) >= conf_threshold)
            pct_conf.append(100 * has_conf.sum() / total_otus if total_otus else 0)
        
        color = db_colors.get(prefix, db_colors.get(db_label, f'C{i}'))
        offset = (i - len(db_datasets)/2 + 0.5) * width
        axes[0].bar(x + offset, pct_conf, width, label=db_label, color=color, alpha=0.85)
    
    axes[0].set_xticks(x); axes[0].set_xticklabels(RANKS, rotation=45)
    axes[0].set_ylabel('% OTUs Confidently Assigned')
    axes[0].set_title(f'{marker_label}: Confident Assignment Rate by Database', fontweight='bold')
    axes[0].legend(); axes[0].set_ylim(0, 105)
    
    # Panel 2: Species-level comparison
    species_data = []
    for df, prefix, db_label in db_datasets:
        col = f'{prefix}_Species'
        conf_col = f'{prefix}_Species_Conf'
        if col in df.columns:
            has = df[col].notna() & (df[col] != '')
            conf = has & (pd.to_numeric(df[conf_col], errors='coerce').fillna(0) >= conf_threshold) if conf_col in df.columns else has
            species_data.append({'DB': db_label, 'Any': 100*has.sum()/len(df), 'Confident': 100*conf.sum()/len(df)})
    
    if species_data:
        sdf = pd.DataFrame(species_data)
        x2 = np.arange(len(sdf))
        axes[1].bar(x2 - 0.2, sdf['Any'], 0.4, label='Any species', alpha=0.6, color='#3498db')
        axes[1].bar(x2 + 0.2, sdf['Confident'], 0.4, label='Confident species', alpha=0.85, color='#2ecc71')
        axes[1].set_xticks(x2); axes[1].set_xticklabels(sdf['DB'])
        axes[1].set_ylabel('% OTUs')
        axes[1].set_title(f'{marker_label}: Species-Level Resolution', fontweight='bold')
        axes[1].legend()
    
    plt.tight_layout()
    plt.show()


def compare_db_phylum_composition(db_datasets, marker_label, sample_cols=None,
                                  top_n=8, conf_threshold=None):
    """Side-by-side phylum composition across databases.
    
    Args:
        db_datasets: list of (df, prefix, db_label) tuples
    """
    if conf_threshold is None:
        conf_threshold = CONF_THRESHOLD
    
    n = len(db_datasets)
    fig, axes = plt.subplots(1, n, figsize=(7*n, 6), sharey=True)
    if n == 1:
        axes = [axes]
    
    for i, (df, prefix, db_label) in enumerate(db_datasets):
        ax = axes[i]
        s_cols = sample_cols if sample_cols else get_sample_cols(df)
        col = f'{prefix}_Phylum'
        conf_col = f'{prefix}_Phylum_Conf'
        d = df.copy()
        if conf_col in d.columns:
            d.loc[pd.to_numeric(d[conf_col], errors='coerce').fillna(0) < conf_threshold, col] = 'Unassigned'
        d[col] = d[col].fillna('Unassigned').replace('', 'Unassigned')
        
        grouped = d.groupby(col)[s_cols].sum().sum(axis=1)
        total = grouped.sum()
        pct = (grouped / total * 100).sort_values(ascending=False)
        
        top = pct.head(top_n)
        others = pct.iloc[top_n:].sum() if len(pct) > top_n else 0
        plot_data = top.copy()
        if others > 0:
            plot_data['Others'] = others
        
        colors = cm.tab20(np.linspace(0, 1, len(plot_data)))
        wedges, texts, autotexts = ax.pie(
            plot_data.values, labels=None, autopct='%1.1f%%',
            colors=colors, startangle=90, pctdistance=0.85)
        ax.set_title(f'{db_label}', fontweight='bold', fontsize=12)
        ax.legend(plot_data.index, loc='center left', bbox_to_anchor=(-0.3, 0.5), fontsize=8)
    
    fig.suptitle(f'{marker_label}: Phylum Composition by Database', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.show()
