import os
import subprocess
from collections import defaultdict
import pysam
import numpy as np

def compute_block_n50(vcf_path, chromosome="chr22"):
    """
    Calcule le N50 de bloc à partir d'un fichier VCF phasé en se basant sur le champ PS.
    """
    blocks = defaultdict(list)
    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf:
            if record.chrom != chromosome:
                continue
            sample = list(record.samples.values())[0]
            if sample.phased and sample.get('PS') is not None:
                blocks[sample['PS']].append(record.pos)

    lengths = []
    for ps, positions in blocks.items():
        if len(positions) >= 2:
            lengths.append(max(positions) - min(positions))
        else:
            lengths.append(1)

    if not lengths:
        return 0

    lengths = sorted(lengths)
    total_len = sum(lengths)
    cum_sum = 0
    for l in lengths:
        cum_sum += l
        if cum_sum >= total_len / 2:
            return l
    return lengths[0]

def parse_whatshap_compare(filepath, vcf_path, chromosome="chr22"):
    """
    Parse le fichier de log de 'whatshap compare' pour extraire les switch errors et la distance de Hamming.
    Calcule le N50 de bloc depuis le VCF.
    """
    metrics = {"switch_errors": 0, "hamming_distance": 0, "block_n50": 0, "phased_variants": 0}

    try:
        metrics["block_n50"] = compute_block_n50(vcf_path, chromosome)
    except Exception as e:
        metrics["block_n50"] = -1

    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
        print(f"⚠️ Warning: File {filepath} not found.")
        return metrics

    with open(filepath, "r") as f:
        in_all_blocks = False
        lines = f.readlines()
        for idx, line in enumerate(lines):
            line_lower = line.lower()

            if "non-singleton blocks in file1:" in line_lower:
                if idx + 1 < len(lines):
                    next_line = lines[idx+1]
                    if "covered variants" in next_line.lower():
                        parts = next_line.split()
                        for p in parts:
                            if p.replace(",", "").isdigit():
                                metrics["phased_variants"] = int(p.replace(",", ""))

            if "all intersection blocks:" in line_lower:
                in_all_blocks = True
            elif "largest intersection block:" in line_lower:
                in_all_blocks = False

            if in_all_blocks:
                if "switch errors" in line_lower:
                    parts = line.split()
                    for p in parts:
                        if p.isdigit():
                            metrics["switch_errors"] = int(p)
                elif "hamming distance" in line_lower:
                    parts = line.split()
                    for p in parts:
                        if p.isdigit():
                            metrics["hamming_distance"] = int(p)
    return metrics

def write_phased_vcf(input_vcf_path, output_vcf_path, active_reads, read_id_to_idx, pred_spins, variants, variant_ps, chromosome="chr22"):
    """
    Exporte le phasing calculé par l'algorithme vers un fichier VCF en utilisant le champ GT et PS.
    """
    with pysam.VariantFile(input_vcf_path) as in_vcf:
        header = in_vcf.header
        with pysam.VariantFile(output_vcf_path, "w", header=header) as out_vcf:
            votes_by_pos = defaultdict(list)
            for rid, profile in active_reads.items():
                if rid not in read_id_to_idx:
                    continue
                u = read_id_to_idx[rid]
                for pos, (allele_val, _) in profile.items():
                    if pos in variants:
                        votes_by_pos[pos].append((u, allele_val))

            phased_count = 0
            for record in in_vcf:
                if record.chrom != chromosome:
                    continue
                pos = record.pos - 1
                if pos in votes_by_pos:
                    votes = votes_by_pos[pos]
                    vote_sum = 0.0
                    for u, allele_val in votes:
                        s_i = 1 - 2 * allele_val
                        vote_sum += s_i * pred_spins[u]

                    sample_name = list(record.samples.keys())[0]
                    if vote_sum >= 0:
                        record.samples[sample_name]['GT'] = (0, 1)
                    else:
                        record.samples[sample_name]['GT'] = (1, 0)
                    record.samples[sample_name].phased = True
                    if 'PS' in in_vcf.header.formats:
                        record.samples[sample_name]['PS'] = variant_ps.get(pos, 1)
                    phased_count += 1
                else:
                    sample_name = list(record.samples.keys())[0]
                    record.samples[sample_name].phased = False
                    if 'PS' in in_vcf.header.formats:
                        record.samples[sample_name]['PS'] = None

                out_vcf.write(record)
            return phased_count
