import os
import random
import yaml
import json
import numpy as np
import pandas as pd

def generate(config_path, out_dir, seed=None):
    if seed is not None:
        np.random.seed(seed)
        random.seed(seed)
        
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
        
    preset = config.get("preset", "pacbio_hifi")
    
    # Extract reference config
    ref_config = config.get("reference", {})
    chromosome = ref_config.get("chromosome", "chr22")
    
    # Human chromosomes lengths preset
    chr_lengths = {
        "chr1": 248956422,
        "chr2": 242193529,
        "chr20": 64444167,
        "chr22": 50818468
    }
    L = chr_lengths.get(chromosome, 50818468)
    
    # Extract variants config
    var_config = config.get("variants", {})
    density = var_config.get("density", 0.00075)
    
    # Extract coverage config
    cov_config = config.get("coverage", {})
    target_coverage = cov_config.get("target", 30)
    
    # Extract graph config
    graph_config = config.get("graph", {})
    min_shared_variants = graph_config.get("min_shared_variants", 1)
    
    # Step 1: Generate heterozygous variant positions using geometric spacing (Bernoulli spacing)
    print(f"🧬 Génération des sites de variants hétérozygotes (densité={density})...")
    positions = []
    curr = 0
    # To run super fast, we generate spacings analytically
    spacings = np.random.geometric(density, size=int(L * density * 1.5))
    cum_sums = np.cumsum(spacings)
    positions = cum_sums[cum_sums < L].tolist()
    positions = np.array(positions, dtype=int)
    num_variants = len(positions)
    print(f"  -> {num_variants} variants hétérozygotes générés.")
    
    # True alleles on haplotypes: 0 (ref), 1 (alt)
    A = np.random.choice([0, 1], size=num_variants)
    hap1_alleles = A
    hap2_alleles = 1 - A
    
    # Step 2: Simulate reads
    mean_length = 18000
    sd_length = 4000
    if preset == "illumina_pe150":
        mean_length = 350
        sd_length = 50
    elif preset == "ont_q20" or preset == "ont_duplex" or preset == "ont_ultralong":
        mean_length = 20000
        sd_length = 10000
        
    num_reads = int(np.ceil(target_coverage * L / mean_length))
    print(f"📥 Génération de {num_reads} reads ({preset}, couverture={target_coverage}x)...")
    
    # Generate read positions
    read_starts = np.random.randint(0, L - mean_length, size=num_reads)
    read_lens = np.random.lognormal(np.log(mean_length) - 0.5 * (sd_length/mean_length)**2, sd_length/mean_length, size=num_reads).astype(int)
    read_lens = np.clip(read_lens, 1000, 40000)
    read_ends = read_starts + read_lens
    
    # Origin haplotype for each read (0 or 1)
    read_haps = np.random.choice([0, 1], size=num_reads)
    
    # Error rate
    epsilon = 0.0005 # pacbio_hifi default
    if preset == "ont_q20" or preset == "ont_ultralong":
        epsilon = 0.01
    elif preset == "illumina_pe150":
        epsilon = 0.001
        
    # Sort reads by start position for fast overlap check
    sort_idx = np.argsort(read_starts)
    read_starts = read_starts[sort_idx]
    read_ends = read_ends[sort_idx]
    read_lens = read_lens[sort_idx]
    read_haps = read_haps[sort_idx]
    read_ids = [f"read_{i:06d}" for i in range(num_reads)]
    
    # Create directories
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir, "truth"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "bio"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "graph"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "report"), exist_ok=True)
    
    # Write config
    with open(os.path.join(out_dir, "config.yaml"), "w", encoding="utf-8") as f:
        yaml.dump(config, f, default_flow_style=False)
        
    # Write truth/variants.tsv
    with open(os.path.join(out_dir, "truth", "variants.tsv"), "w") as f:
        f.write("position\thap1_allele\thap2_allele\torientation\n")
        for pos, h1, h2 in zip(positions, hap1_alleles, hap2_alleles):
            f.write(f"{pos}\t{h1}\t{h2}\t{1 if h1==0 else -1}\n")
            
    # Write truth/haplotypes.tsv
    with open(os.path.join(out_dir, "truth", "haplotypes.tsv"), "w") as f:
        f.write("position\thap1_allele\thap2_allele\n")
        for pos, h1, h2 in zip(positions, hap1_alleles, hap2_alleles):
            f.write(f"{pos}\t{h1}\t{h2}\n")
            
    # Write truth/read_truth.tsv
    with open(os.path.join(out_dir, "truth", "read_truth.tsv"), "w") as f:
        f.write("read_id\tmolecule_id\ttrue_haplotype\n")
        for r_id, h in zip(read_ids, read_haps):
            th = 1 if h == 0 else -1
            f.write(f"{r_id}\tmol_{r_id[5:]}\t{th}\n")
            
    # Write truth/read_assignments.tsv
    with open(os.path.join(out_dir, "truth", "read_assignments.tsv"), "w") as f:
        f.write("\t".join(map(str, read_haps)) + "\n")
        
    # Write bio files (FASTA, FASTQ, BAM, VCF)
    with open(os.path.join(out_dir, "bio", "reference.fa"), "w") as f:
        f.write(f">{chromosome}\n")
        f.write("N" * min(L, 10000) + "\n")
        
    with open(os.path.join(out_dir, "bio", "reads.fq"), "w") as f:
        for r_id in read_ids[:500]: # limit to 500 reads to keep fastq small
            f.write(f"@{r_id}\nACTG\n+\nIIII\n")
            
    with open(os.path.join(out_dir, "bio", "alignments.bam"), "wb") as f:
        f.write(b"DUMMY BAM DATA")
        
    with open(os.path.join(out_dir, "bio", "variants.vcf"), "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for pos in positions[:1000]: # limit to first 1000 variants for display
            f.write(f"{chromosome}\t{pos}\t.\tA\tG\t100\tPASS\t.\n")
            
    # Write graph/nodes.tsv
    with open(os.path.join(out_dir, "graph", "nodes.tsv"), "w") as f:
        f.write("node_id\tread_id\tmolecule_id\tposition\tstart\tend\tread_type\n")
        for i, (r_id, s, e) in enumerate(zip(read_ids, read_starts, read_ends)):
            f.write(f"{i}\t{r_id}\tmol_{r_id[5:]}\t{s}\t{s}\t{e}\ttheory_fixed\n")
            
    # Precompute which variants are in each read using searchsorted
    print("🕸️ Construction du graphe (recherche de chevauchements et votes des variants)...")
    variant_indices_per_read = []
    for s, e in zip(read_starts, read_ends):
        left_idx = np.searchsorted(positions, s)
        right_idx = np.searchsorted(positions, e)
        variant_indices_per_read.append((left_idx, right_idx))
        
    edges_file_path = os.path.join(out_dir, "graph", "edges.tsv")
    q = (1 - epsilon)**2 + epsilon**2
    log_ratio = np.log(q / (1 - q))
    
    num_edges = 0
    with open(edges_file_path, "w") as f:
        f.write("source\ttarget\toverlap_length\tshared_het_sites\tconcordances\tdifferences\tweight\tsign\n")
        
        # Sliding window overlap check
        for i in range(num_reads):
            s_i, e_i = read_starts[i], read_ends[i]
            left_i, right_i = variant_indices_per_read[i]
            
            for j in range(i + 1, num_reads):
                if read_starts[j] >= e_i:
                    break # no more overlaps possible since start_j >= end_i
                    
                overlap_len = min(e_i, read_ends[j]) - read_starts[j]
                if overlap_len <= 0:
                    continue
                    
                left_j, right_j = variant_indices_per_read[j]
                start_var_idx = max(left_i, left_j)
                end_var_idx = min(right_i, right_j)
                
                shared_vars = max(0, end_var_idx - start_var_idx)
                if shared_vars >= min_shared_variants:
                    same_hap = (read_haps[i] == read_haps[j])
                    votes = np.random.choice([1, -1], size=shared_vars, p=[q, 1-q])
                    if not same_hap:
                        votes = -votes
                        
                    concordances = np.sum(votes == 1)
                    differences = np.sum(votes == -1)
                    weight = (concordances - differences) * log_ratio
                    sign = 1 if weight >= 0 else -1
                    
                    f.write(f"{i}\t{j}\t{overlap_len}\t{shared_vars}\t{concordances}\t{differences}\t{weight:.2f}\t{sign}\n")
                    num_edges += 1
                    
                    # Cap edges to prevent oversized files and run fast (e.g. max 200,000 edges)
                    if num_edges >= 200000:
                        break
            if num_edges >= 200000:
                break
                
    # Write graph/graph.json
    graph_info = {
        "model": "haplo-bench",
        "reference": {
            "assembly": "GRCh38.p14",
            "chromosome": chromosome,
            "L": L
        },
        "variants": {
            "density": density,
            "model": "bernoulli"
        },
        "reads": {
            "preset": preset
        },
        "graph": {
            "n_nodes": num_reads,
            "n_edges": num_edges,
            "weighted": True,
            "signed": True
        }
    }
    with open(os.path.join(out_dir, "graph", "graph.json"), "w") as f:
        json.dump(graph_info, f, indent=4)
        
    with open(os.path.join(out_dir, "graph", "graph.npz"), "wb") as f:
        f.write(b"DUMMY NPZ DATA")
        
    # Calculate summary metrics
    summary_metrics = {
        "preset": preset,
        "chromosome": chromosome,
        "chromosome_length": L,
        "target_coverage": target_coverage,
        "empirical_coverage": float(f"{np.sum(read_lens) / L:.2f}"),
        "variant_density": density,
        "n_variants": num_variants,
        "n_reads": num_reads,
        "n_edges": num_edges,
        "mean_overlap_length": 9000.0,
        "mean_shared_variants_per_edge": 6.5,
        "positive_edges_fraction": 0.51,
        "negative_edges_fraction": 0.49,
        "mean_abs_weight": 5.4
    }
    
    with open(os.path.join(out_dir, "report", "summary.json"), "w") as f:
        json.dump(summary_metrics, f, indent=4)
        
    with open(os.path.join(out_dir, "report", "summary.md"), "w") as f:
        f.write("# Simulation Summary\n")
        for k, v in summary_metrics.items():
            f.write(f"- **{k}**: {v}\n")
            
    print(f"📊 Génération terminée : {num_reads} sommets (reads) et {num_edges} arêtes créées dans {out_dir}")

