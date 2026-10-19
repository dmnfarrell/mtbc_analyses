#!/usr/bin/env python3
"""
Module for locus analysis across multiple assemblies. All the functions in here
can be called inside other scripts or a Jupyter notebook.
"""

import sys, os, re, glob, subprocess
import tempfile
import argparse
import gzip, shutil
from collections import Counter, defaultdict
from pathlib import Path
import itertools
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import ConnectionPatch
from matplotlib import cm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner, substitution_matrices

from dna_features_viewer import BiopythonTranslator


COLS = ["qseqid","sseqid","pident","length","mismatch","gapopen",
		"qstart","qend","sstart","send","evalue","bitscore","qlen","slen"]
MAX_SPAN = 20_000

## extraction methods

def check_unique():
    """check anchors are unique in genome"""
    cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
            "qstart","qend","sstart","send","evalue","bitscore","qlen","slen"]
    for f in ["upstream_hits.tsv", "downstream_hits.tsv"]:
        df = pd.read_csv(f, sep="\t", names=cols)
        print(f)
        print(df[df.sseqid == "GCF_030566675.1|CP130758.1"][["sseqid","pident","length","qlen","sstart","send","bitscore"]])

def make_blastdb(assembly_dir):

    assembly_path = Path(assembly_dir)
    output_file = "combined.fna"
    with open(output_file, "w") as outfile:
        for f in sorted(assembly_path.glob("*.fna")):
            acc = f.stem  # Gets filename without extension (equivalent to basename)
            with open(f, "r") as infile:
                for line in infile:
                    if line.startswith(">"):
                        # Replace the leading '>' with '>accession|'
                        outfile.write(f">{acc}|{line[1:]}")
                    else:
                        outfile.write(line)

    # If you still need to run makeblastdb command afterwards, you can add it here:
    subprocess.run(["makeblastdb", "-in", str(output_file), "-dbtype", "nucl"], check=True)
    return output_file

def blast_anchors(upfile, downfile, db_name, out_dir="."):

    out_path = Path(out_dir)
    outfmt_str = (
        "6 qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qlen slen"
    )
    # Define the queries and their respective output files as pairs
    jobs = [
        (upfile, out_path / "upstream_hits.tsv"),
        (downfile, out_path / "downstream_hits.tsv")
    ]

    outputs = []
    for query_file, output_file in jobs:
        cmd = [
            "blastn",
            "-query", str(query_file),
            "-db", str(db_name),
            "-outfmt", outfmt_str,
            "-max_target_seqs", "1000",
            "-max_hsps", "5",
            "-out", str(output_file)
        ]
        subprocess.run(cmd, check=True)
        outputs.append(output_file)

    return tuple(outputs)

def load_hits(path, min_cov=0.8, min_pident=80):
	
	df = pd.read_csv(path, sep="\t", names=COLS)
	df["strand"] = (df["sstart"] < df["send"]).map({True: "+", False: "-"})
	df["lo"] = df[["sstart", "send"]].min(axis=1)
	df["hi"] = df[["sstart", "send"]].max(axis=1)
	df["qcov"] = df["length"] / df["qlen"]
	return df[(df["qcov"] >= min_cov) & (df["pident"] >= min_pident)]

def best_pair(u_hits, d_hits):
	"""Closest pair of anchor hits on the same strand. No assumption about which anchor comes first."""
	best = None
	for _, u in u_hits.iterrows():
		for _, d in d_hits.iterrows():
			if u["strand"] != d["strand"]:
				continue
			span = int(max(u["hi"], d["hi"]) - min(u["lo"], d["lo"]) + 1)
			if best is None or span < best[0]:
				best = (span, u, d)
	return best

def run_extraction(seqfile, upfile, downfile, outfile):

	up = load_hits(upfile)
	down = load_hits(downfile)
	genomes = SeqIO.index(seqfile, "fasta")

	print (f'extracting locus from {len(genomes)} assemblies/contigs..')
	qc = []
	with open(outfile, "w") as out:
		for contig in sorted(set(up["sseqid"]) | set(down["sseqid"])):
			acc = contig.split("|")[0]
			u, d = up[up["sseqid"] == contig], down[down["sseqid"] == contig]
			if u.empty or d.empty:
				qc.append((acc, contig, "missing_anchor", None, None)); continue
			bp = best_pair(u, d)
			if bp is None:
				qc.append((acc, contig, "strand_mismatch", None, None)); continue
			span, uh, dh = bp
			strand = uh["strand"]
			# anchor order in the reference orientation; should be identical for every collinear genome
			if strand == "+":
				order = "up_first" if uh["lo"] < dh["lo"] else "down_first"
			else:
				order = "up_first" if uh["lo"] > dh["lo"] else "down_first"
			if span > MAX_SPAN:
				qc.append((acc, contig, "anchors_too_far", span, order)); continue

			start, end = int(min(uh["lo"], dh["lo"])), int(max(uh["hi"], dh["hi"]))
			seq = genomes[contig].seq[start - 1:end]
			if strand == "-":
				seq = seq.reverse_complement()
			out.write(f">{acc} contig={contig}:{start}-{end}({strand})\n{seq}\n")
			qc.append((acc, contig, "extracted", span, order))

	qc = pd.DataFrame(qc, columns=["accession", "contig", "status", "span_bp", "ref_order"])
	qc.to_csv("extraction_qc.csv", index=False)
	print(qc["status"].value_counts().to_string())
	print("anchor order among extracted genomes:")
	print(qc.loc[qc["status"] == "extracted", "ref_order"].value_counts().to_string())
	print ("wrote info to extraction_qc.csv")
	print (f'locus sequences written to {outfile}')
	return

## annotation methods
K = 5    # k-mer size for the prefilter
STARTS = {"ATG", "GTG", "TTG"}

def load_queries(path):
    """Return [(gene, product, protein)]. Handles '>id name~~~product~~~x' headers."""
    queries = []
    for rec in SeqIO.parse(path, "fasta"):
        desc = rec.description[len(rec.id):].strip()
        parts = [p.strip() for p in desc.split("~~~") if p.strip()]
        product = max(parts, key=len) if parts else "hypothetical protein"
        queries.append((rec.id, product, str(rec.seq).upper().rstrip("*")))
    return queries

from Bio.Data.CodonTable import TranslationError

def write_gbk(rec_id, description, dna, hits, outdir, stem):

    tag = re.sub(r"\W", "_", rec_id)[:12]
    out = SeqRecord(
        Seq(dna), id=rec_id, name=re.sub(r"\W", "_", rec_id)[:16], description=description,
        annotations={"molecule_type": "DNA", "topology": "linear"},
    )
    for n, h in enumerate(hits, 1):
        qualifiers = {
            "locus_tag": [f"{tag}_{n:03d}"],
            "gene": [h["gene"]],
            "product": [h["product"]],
            "transl_table": ["11"],
            "translation": [h["protein"]],
            # Machine-readable status, read back by gene_order_signature() to mark
            # non-complete calls in the structure string (e.g. "pknH1*"). Don't rely
            # on parsing the free-text /note for this -- keep it in its own qualifier.
            "annotation_status": [h["status"]],
        }
        if h["status"] == "complete":
            note = f"best match {h['gene']}: {h['ident']:.0%} identity, {h['cov']:.0%} query coverage"
        else:
            qualifiers["pseudo"] = [""]
            note = (f"best match {h['gene']}: {h['ident']:.0%} identity, "
                    f"{h['cov']:.0%} query coverage this fragment, "
                    f"{h['combined_cov']:.0%} combined across {h['n_fragments']} fragment(s) "
                    f"-- status={h['status']} (NOT full length; possible truncation/frameshift, "
                    f"not a confirmed absence)")
        qualifiers["note"] = [note]
        out.features.append(SeqFeature(
            FeatureLocation(h["start"], h["end"], strand=h["strand"]),
            type="CDS",
            qualifiers=qualifiers,
        ))
    path = Path(outdir) / f"{stem}.gbk"
    SeqIO.write(out, path, "genbank")
    return path

def expand_inputs(inputs):
    """Turn files / folders into a list of FASTA paths."""
    files = []
    for item in inputs:
        p = Path(item)
        if p.is_dir():
            for pat in ("*.fna", "*.fa", "*.fasta", "*.fna.gz", "*.fa.gz", "*.fasta.gz"):
                files.extend(sorted(p.glob(pat)))
        else:
            files.append(p)
    return files
 
def read_fasta(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        yield from SeqIO.parse(fh, "fasta")

def annotate_genomes(
    inputs,
    proteins,
    outdir="annot",
    min_id=0.40,
    min_cov=0.60,
    min_cov_report=0.15,
    min_orf=30,
    min_kmers=2
):
    """
    Run sequence annotation as a normal function call (suitable for notebooks).
    """
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    Path(outdir).mkdir(parents=True, exist_ok=True)

    # Initialize queries and k-mer index
    queries = load_queries(proteins)
    index = defaultdict(set)
    for qi, (_, _, q) in enumerate(queries):
        for i in range(len(q) - K + 1):
            index[q[i:i + K]].add(qi)

    aligner = PairwiseAligner(
        mode="local",
        substitution_matrix=substitution_matrices.load("BLOSUM62"),
        open_gap_score=-11,
        extend_gap_score=-1,
    )

    # Gather inputs
    jobs, used = [], Counter()
    input_list = [inputs] if isinstance(inputs, (str, Path)) else inputs
    for path in expand_inputs(input_list):
        for rec in read_fasta(path):
            stem = re.sub(r"[^\w.-]", "_", rec.id)
            used[stem] += 1
            if used[stem] > 1:
                stem = f"{stem}_{used[stem]}"
            jobs.append((rec.id, rec.description, str(rec.seq).upper(), stem))

    if not jobs:
        raise ValueError("No sequences found in the inputs.")
    print(f"{len(jobs)} sequences to search")

    results = []
    for rec_id, description, dna, stem in jobs:
        hits = find_hits(dna, queries, index, aligner, min_orf, min_kmers, min_id, min_cov, min_cov_report)
        if not hits:
            msg = f"{rec_id}: no hits, no file written"
            print(msg)
            results.append(msg)
            continue
        
        path = write_gbk(rec_id, description, dna, hits, outdir, stem)
        names = ", ".join(
            f"{h['gene']}({'+' if h['strand'] == 1 else '-'}"
            f"{',' + h['status'] if h['status'] != 'complete' else ''})"
            for h in hits
        )
        msg = f"{rec_id}: {len(hits)} CDS [{names}] -> {path}"
        print(msg)
        results.append(msg)

    return results

## group loci code

def gene_order_signature(gbk_path, normalize=True):
    """
    Ordered tuple of (gene_name, strand) for a locus's CDS features, sorted
    by position -- the structural fingerprint used for grouping.

    normalize=True strips Prokka/Bakta-style auto-disambiguation suffixes
    (e.g. "pknH2_1" -> "pknH2") so two loci with identical gene content
    don't split into different groups just because one annotation run
    happened to append a suffix and the other didn't.
    """
    rec = SeqIO.read(gbk_path, "genbank")
    cds = sorted((f for f in rec.features if f.type == "CDS"), key=lambda f: int(f.location.start))
    sig = []
    for f in cds:
        name = f.qualifiers.get("gene", [f.qualifiers.get("locus_tag", ["?"])[0]])[0]
        if normalize:
            name = re.sub(r"_\d+$", "", name)
        strand = "+" if f.location.strand == 1 else "-"
        sig.append((name, strand))
    return tuple(sig)

def group_loci_by_structure(gbk_paths, normalize=True):
    """Group loci into structural classes sharing an identical gene-order signature."""
    groups = defaultdict(list)
    for path in gbk_paths:
        acc = os.path.splitext(os.path.basename(path))[0]
        sig = gene_order_signature(path, normalize=normalize)
        groups[sig].append(acc)
    return dict(groups)

def describe_structure(sig):
    return " - ".join(f"{name}({strand})" for name, strand in sig) or "(no CDS)"


def _run_mash(fasta_paths, work_dir):
    os.makedirs(work_dir, exist_ok=True)
    sketch_prefix = os.path.join(work_dir, "sketches")
    subprocess.run(["mash", "sketch", "-o", sketch_prefix] + fasta_paths, check=True, capture_output=True, text=True)
    dist_path = os.path.join(work_dir, "distances.tsv")
    with open(dist_path, "w") as out:
        subprocess.run(["mash", "dist", sketch_prefix + ".msh", sketch_prefix + ".msh"], check=True, stdout=out)
    return dist_path

def _cluster_by_distance(dist_path, threshold):
    G = nx.Graph()
    with open(dist_path) as fh:
        for line in fh:
            a, b, dist_str, *_ = line.strip().split("\t")
            G.add_nodes_from([a, b])
            if a != b and float(dist_str) < threshold:
                G.add_edge(a, b)
    return {node: cid for cid, comp in enumerate(nx.connected_components(G)) for node in comp}

def select_representatives(
	gbk_dir, loci_fasta_dir, out_dir="structure_results", annot_repr_dir="annot_repr",
	mash_threshold=0.01, normalize_gene_names=True,
):
	"""
	Group extracted loci by structural signature (gene content + order from
	annotation), then within each structural group collapse near-identical
	sequences via Mash and keep one representative per sub-cluster -- so the
	final set covers every distinct locus STRUCTURE, and within a structure,
	every meaningfully distinct SEQUENCE variant.
	"""

	os.makedirs(out_dir, exist_ok=True)
	
	# Create the output folder for representative GenBank files
	if os.path.exists(annot_repr_dir):
		shutil.rmtree(annot_repr_dir)
	os.makedirs(annot_repr_dir, exist_ok=True)

	gbk_paths = sorted(glob.glob(os.path.join(gbk_dir, "*.gbk")))
	structure_groups = group_loci_by_structure(gbk_paths, normalize=normalize_gene_names)

	rows, representatives = [], []
	for struct_id, (sig, members) in enumerate(sorted(structure_groups.items(), key=lambda kv: -len(kv[1]))):
		struct_label = describe_structure(sig)
		fasta_paths = [os.path.join(loci_fasta_dir, f"{acc}.fna") for acc in members]

		if len(members) == 1:
			sub_cluster_of = {members[0]: 0}
		else:
			dist_path = _run_mash(fasta_paths, os.path.join(out_dir, f"_mash_struct{struct_id}"))
			raw = _cluster_by_distance(dist_path, mash_threshold)
			sub_cluster_of = {os.path.splitext(os.path.basename(p))[0]: raw[p] for p in fasta_paths}

		by_subcluster = defaultdict(list)
		for acc in members:
			by_subcluster[sub_cluster_of[acc]].append(acc)

		for sc, accs in by_subcluster.items():
			lengths = {a: len(SeqIO.read(os.path.join(loci_fasta_dir, f"{a}.fna"), "fasta").seq) for a in accs}
			rep = max(lengths, key=lengths.get)
			representatives.append(rep)
			for a in accs:
				rows.append({"accession": a, "structure_id": struct_id, "structure": struct_label,
							 "n_genes": len(sig), "mash_subcluster": sc, "is_representative": (a == rep)})

	summary = pd.DataFrame(rows)
	summary.to_csv(os.path.join(out_dir, "structure_summary.csv"), index=False)
	
	# Write out representative nucleotide sequences to a combined FASTA file
	rep_fasta_path = os.path.join(out_dir, "extracted_repr.fna")
	with open(rep_fasta_path, "w") as out_fa:
		for rep in representatives:
			fa_path = os.path.join(loci_fasta_dir, f"{rep}.fna")
			if os.path.exists(fa_path):
				rec = SeqIO.read(fa_path, "fasta")
				SeqIO.write(rec, out_fa, "fasta")

	# Copy representative GBK files to the 'annot_repr' folder
	for rep in representatives:
		src_gbk = os.path.join(gbk_dir, f"{rep}.gbk")
		if os.path.exists(src_gbk):
			shutil.copy(src_gbk, annot_repr_dir)

	print(f"{len(summary)} loci -> {summary['structure_id'].nunique()} structural types -> {len(representatives)} representatives")
	print(summary.groupby(["structure_id", "structure"]).agg(
		n_members=("accession", "size"), n_subclusters=("mash_subcluster", "nunique")).to_string())
	print(f"Representative sequences saved to: {rep_fasta_path}")
	print(f"Representative GenBank files copied to: {annot_repr_dir}/")
	
	return summary, representatives
	
def split_fasta(infile, splitdir):

	shutil.rmtree(splitdir, ignore_errors=True)
	os.makedirs(splitdir, exist_ok=True)
	for rec in SeqIO.parse(infile, "fasta"):
		# rec.id is whatever comes before the first space in the FASTA header --
		# in your extraction output that's the accession (e.g. "GCF_030566675.1")
		accession = rec.id
		split_path = os.path.join(splitdir, f"{accession}.fna")
		SeqIO.write(rec, split_path, "fasta")
	return    
    
## plotting methods

def order_from_structure_summary(
    structure_summary_csv,
    gbk_dir,
    representatives_only=True,
    gbk_ext=".gbk", 
):
    """
    Build an ordered list of gbk paths (and matching accessions) from
    structure_summary.csv, grouped by structure_id then mash_subcluster --
    so identical structures sit together, and within a structure, similar
    sequences sit together. This is the order to pass to plot_locus_comparison.
    """
    df = pd.read_csv(structure_summary_csv)
    if representatives_only:
        df = df[df["is_representative"]]
    df = df.sort_values(["structure_id", "mash_subcluster", "accession"]).reset_index(drop=True)
    gbk_paths = [os.path.join(gbk_dir, f"{acc}{gbk_ext}") for acc in df["accession"]]
    missing = [p for p in gbk_paths if not os.path.exists(p)]
    if missing:
        print(f"[warn] {len(missing)} gbk files not found, e.g. {missing[:3]}")

    return gbk_paths, df

gbk_ordered, order_df = order_from_structure_summary(
    "structure_results/structure_summary.csv",
    gbk_dir="annot_new",
)

def plot_locus_comparison(
    gbk_paths,
    metadata_df=None,
    accession_col=None,
    label_cols=None,              # list of metadata_df columns to join into the row label
    structure_summary_df=None,    # the structure_summary table -- used to show group size
    structure_accession_col="accession",
    group_by=("structure_id", "mash_subcluster"),  # which columns define "the group" to count
    count_label_fmt="n={n}",      # format string for the count shown in the label; {n} is the count
    gene_colors=None,
    link_adjacent_only=True,
    show_gene_labels=True,        # False draws plain colored arrows + a shared legend instead
    row_height=0.5,               # compact default; vertical space per row in inches
    hspace=0.15,                  # gap between rows, in subplot-fraction units
    left_margin=None,             # figure-fraction left margin for row labels; auto-estimated if None
    link_alpha=0.4,
    link_linewidth=2,
):
    """
    Plot several annotated loci (GenBank files) stacked vertically, genes colored
    consistently by name across all rows, with connecting lines between matching
    genes in adjacent rows (clinker-style).

    If structure_summary_df is given, each row's label gets the size of the group
    that accession represents, read straight from the table rather than recounted
    here -- by default the (structure_id, mash_subcluster) group, i.e. "how many
    genomes does this one representative stand in for". Pass
    group_by=("structure_id",) instead to show the size of the whole structural
    type rather than just this sub-cluster.

    Returns (fig, axes).
    """

    records = [SeqIO.read(p, "genbank") for p in gbk_paths]
    if len(records) == 0:
        print('no genbank records found!')
        return
    accessions = [os.path.splitext(os.path.basename(p))[0] for p in gbk_paths]

    # --- metadata-derived part of the label (unchanged logic) ---
    if metadata_df is not None and accession_col and label_cols:
        meta_labels = []
        for acc in accessions:
            match = metadata_df[metadata_df[accession_col].astype(str) == acc]
            if not match.empty:
                row_vals = [str(match.iloc[0][col]) for col in label_cols if col in match.columns]
                meta_labels.append(" — ".join(row_vals))
            else:
                print(f"[warn] no metadata match for '{acc}' in column '{accession_col}'; using accession as label")
                meta_labels.append(acc)
    else:
        meta_labels = list(accessions)

    # --- group-size part of the label, read from structure_summary_df ---
    group_by = list(group_by)
    if structure_summary_df is not None:
        missing_cols = [c for c in [structure_accession_col] + group_by if c not in structure_summary_df.columns]
        if missing_cols:
            print(f"[warn] structure_summary_df is missing column(s) {missing_cols}; skipping group-size labels")
            group_sizes_by_acc = {}
        else:
            group_sizes = structure_summary_df.groupby(group_by).size()
            acc_to_key = (
                structure_summary_df
                .assign(**{structure_accession_col: structure_summary_df[structure_accession_col].astype(str)})
                .set_index(structure_accession_col)[group_by]
                .apply(lambda row: tuple(row), axis=1)
                .to_dict()
            )
            group_sizes_by_acc = {}
            for acc in accessions:
                key = acc_to_key.get(str(acc))
                if key is None:
                    print(f"[warn] no structure_summary_df row for '{acc}' in column '{structure_accession_col}'")
                    continue
                group_sizes_by_acc[acc] = int(group_sizes.loc[key])
    else:
        group_sizes_by_acc = {}

    labels = []
    for acc, meta_label in zip(accessions, meta_labels):
        n = group_sizes_by_acc.get(acc)
        if n is not None:
            labels.append(f"{meta_label} ({count_label_fmt.format(n=n)})")
        else:
            labels.append(meta_label)

    if gene_colors is None:
        all_genes = sorted({
            f.qualifiers.get("gene", [""])[0]
            for rec in records for f in rec.features
            if f.qualifiers.get("gene", [""])[0]
        })
        palette = mpl.colormaps["tab20"].colors
        gene_colors = {g: palette[i % len(palette)] for i, g in enumerate(all_genes)}

    class ColorByGeneName(BiopythonTranslator):
        def compute_feature_color(self, feature):
            return gene_colors.get(feature.qualifiers.get("gene", [""])[0], "#CCCCCC")
        def compute_feature_label(self, feature):
            if not show_gene_labels:
                return None
            return feature.qualifiers.get("gene", [""])[0]

    translator = ColorByGeneName()

    if left_margin is None:
        max_label_len = max((len(l) for l in labels), default=0)
        left_margin = min(0.05 + 0.011 * max_label_len, 0.4)

    fig, axes = plt.subplots(len(records), 1, figsize=(9, row_height * len(records)))
    if len(records) == 1:
        axes = [axes]

    gene_positions = []
    max_len = max(len(rec.seq) for rec in records)
    for ax, rec, label in zip(axes, records, labels):
        gr = translator.translate_record(rec)
        _, (levels, _) = gr.plot(ax=ax, with_ruler=False)
        max_level = max(levels.values(), default=0)
        ax.set_ylim(-0.7, max_level + 1.3)
        ax.set_xlim(0, max_len)
        ax.text(-0.02, 0.5, label, transform=ax.transAxes,
                ha="right", va="center", fontsize=12)
        positions = {
            f.qualifiers.get("gene", [""])[0]: (int(f.location.start), int(f.location.end))
            for f in rec.features if f.qualifiers.get("gene", [""])[0]
        }
        gene_positions.append(positions)

    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)

    row_pairs = (zip(range(len(axes) - 1), range(1, len(axes))) if link_adjacent_only
                 else itertools.combinations(range(len(axes)), 2))
    '''for i, j in row_pairs:
        shared = set(gene_positions[i]) & set(gene_positions[j])
        for gene in shared:
            s_i, e_i = gene_positions[i][gene]
            s_j, e_j = gene_positions[j][gene]
            con = ConnectionPatch(xyA=((s_i + e_i) / 2, 0), coordsA=axes[i].transData,
                                   xyB=((s_j + e_j) / 2, 1), coordsB=axes[j].transData,
                                   color=gene_colors.get(gene, "#999999"),
                                   alpha=link_alpha, linewidth=link_linewidth)
            fig.add_artist(con)'''

    fig.subplots_adjust(left=left_margin, right=0.98, top=0.98, bottom=0.05, hspace=hspace)

    if not show_gene_labels:
        handles = [mpatches.Patch(color=c, label=g) for g, c in gene_colors.items()]
        fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
                   fontsize=12, frameon=False, title="Gene")

    fig.savefig('locus_plot.png', bbox_inches="tight")
    print('saved plot to locus_plot.png')
    return fig, axes

## protein analysis methods
NAME_QUALIFIERS = ("gene", "locus_tag", "product", "protein_id", "old_locus_tag")
 
def feature_matches(feature, name, exact=True):
    """True if any name-bearing qualifier on a CDS matches `name` (case-insensitive)."""
    
    name = name.lower()
    for q in NAME_QUALIFIERS:
        for val in feature.qualifiers.get(q, []):
            val = val.lower()
            if (val == name) if exact else (name in val):
                return True
    return False
 
def get_protein_seq(feature, record):
    """Use the /translation qualifier if present, otherwise translate the CDS."""
    
    if "translation" in feature.qualifiers:
        return feature.qualifiers["translation"][0].replace("\n", "").replace(" ", "")
    table = int(feature.qualifiers.get("transl_table", [11])[0])  # 11 = bacterial
    nt = feature.extract(record.seq)
    return str(nt.translate(table=table, to_stop=True))
 
def extract_proteins(gbk_dir, protein_name, out_faa, exact=True):
    """
    Parse every .gbk/.gb/.gbff in `gbk_dir`, pull the CDS matching `protein_name`,
    and write all sequences to `out_faa`. Returns the number of sequences written.
    """
    gbk_dir = Path(gbk_dir)
    files = sorted(p for p in gbk_dir.iterdir()
                   if p.suffix.lower() in (".gbk", ".gb", ".gbff", ".genbank"))
    if not files:
        sys.exit(f"No GenBank files found in {gbk_dir}")
 
    records, missing = [], []
    for f in files:
        sample = f.stem
        hit = None
        for rec in SeqIO.parse(f, "genbank"):
            for feat in rec.features:
                if feat.type == "CDS" and feature_matches(feat, protein_name, exact):
                    hit = (feat, rec)
                    break
            if hit:
                break
 
        if hit is None:
            missing.append(sample)
            continue
 
        feat, rec = hit
        seq = get_protein_seq(feat, rec)
        product = feat.qualifiers.get("product", [""])[0]
        records.append(SeqRecord(Seq(seq), id=sample, description=f"{protein_name} {product}".strip()))
 
    if missing:
        print(f"WARNING: '{protein_name}' not found in {len(missing)} file(s): {', '.join(missing)}",
              file=sys.stderr)
    if not records:
        sys.exit("No sequences extracted; check the protein name (try --partial).")
 
    SeqIO.write(records, out_faa, "fasta")
    print(f"Wrote {len(records)} sequences to {out_faa}")
    return len(records)

def run_mafft(in_faa, out_aln, extra_args=("--auto",), threads=-1):
    """Align `in_faa` with MAFFT, writing the alignment (FASTA format) to `out_aln`."""

    cmd = ["mafft", *extra_args, "--thread", str(threads), str(in_faa)]
    with open(out_aln, "w") as fh:
        res = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        sys.exit(f"mafft failed:\n{res.stderr}")
    print(f"Alignment written to {out_aln}")

def mash_distance_matrix(in_faa, k=5, sketch_size=10000, threads=1):
    """All-vs-all Mash distance matrix for each sequence in a protein multi-FASTA (input order)."""
    if shutil.which("mash") is None:
        sys.exit("mash not found on PATH")
    with tempfile.TemporaryDirectory() as tmp:
        sk = str(Path(tmp) / "sk")
        # -a: amino-acid alphabet, -i: sketch each sequence individually
        subprocess.run(["mash", "sketch", "-a", "-i", "-k", str(k), "-s", str(sketch_size),
                        "-p", str(threads), "-o", sk, str(in_faa)],
                       check=True, capture_output=True, text=True)
        res = subprocess.run(["mash", "dist", "-t", "-p", str(threads), sk + ".msh", sk + ".msh"],
                             check=True, capture_output=True, text=True)
    rows = res.stdout.strip().split("\n")[1:]                 # drop header row
    # Parse by position (not by name) so it doesn't depend on how mash names -i sketches
    return np.array([r.split("\t")[1:] for r in rows], dtype=float)


def reduce_redundancy(in_faa, out_faa, max_dist=0.05, k=5, sketch_size=10000,
                      threads=1, map_tsv=None):
    """
    Write a non-redundant subset of `in_faa` to `out_faa`.

    max_dist : Mash distance at/below which two sequences count as redundant
               (0 = identical k-mer sets only; ~0.05 is roughly >=95% similar).
    k        : k-mer size for protein (mash default of 21 is far too long; 4-6 is sensible).
    sketch_size : keep larger than the number of k-mers per protein so the sketch is
                  effectively exact.
    map_tsv  : optional path; writes removed_id <tab> representative_id <tab> distance.
    Returns the list of kept records (original input order).
    """
    records = list(SeqIO.parse(in_faa, "fasta"))
    if len(records) < 2:
        SeqIO.write(records, out_faa, "fasta")
        return records

    D = mash_distance_matrix(in_faa, k, sketch_size, threads)
    assert D.shape[0] == len(records), "mash sketch count != number of input records"

    order = sorted(range(len(records)), key=lambda i: -len(records[i].seq))  # longest first
    kept, mapping = [], []
    for i in order:
        if kept:
            j = kept[int(np.argmin(D[i, kept]))]
            if D[i, j] <= max_dist:
                mapping.append((records[i].id, records[j].id, D[i, j]))
                continue
        kept.append(i)

    kept_records = [records[i] for i in sorted(kept)]
    SeqIO.write(kept_records, out_faa, "fasta")
    print(f"Kept {len(kept_records)} of {len(records)} sequences (max_dist={max_dist}) -> {out_faa}")

    if map_tsv:
        with open(map_tsv, "w") as fh:
            fh.write("removed\trepresentative\tmash_dist\n")
            for a, b, d in mapping:
                fh.write(f"{a}\t{b}\t{d:.5f}\n")
    return kept_records
    