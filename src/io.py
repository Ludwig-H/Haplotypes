import urllib.request
import pysam
import os
import random
import bisect
from collections import defaultdict

import sys

def download_file(url, local_path):
    """Télécharge un fichier depuis une URL vers un chemin local avec une barre de progression en temps réel."""
    if not os.path.exists(local_path):
        print(f"📥 Démarrage du téléchargement de {url} -> {local_path}...")
        
        def reporthook(block_num, block_size, total_size):
            if total_size <= 0:
                return
            downloaded = block_num * block_size
            percent = min(100.0, downloaded * 100.0 / total_size)
            bar_len = 40
            filled_len = int(round(bar_len * percent / 100.0))
            bar = '█' * filled_len + '-' * (bar_len - filled_len)
            
            downloaded_mb = downloaded / (1024 * 1024)
            total_mb = total_size / (1024 * 1024)
            
            sys.stdout.write(f"\r   [{bar}] {percent:.1f}% ({downloaded_mb:.1f}/{total_mb:.1f} Mo)")
            sys.stdout.flush()
            
        urllib.request.urlretrieve(url, local_path, reporthook=reporthook)
        sys.stdout.write("\n✅ Téléchargement terminé avec succès.\n")
        sys.stdout.flush()
    else:
        print(f"✅ {local_path} déjà présent localement.")

def prepare_cleaned_vcfs(input_phased_path, output_unphased_path, output_phased_chr_path, chromosome):
    """Crée les VCFs nettoyés (unphased et phased) pour le chromosome cible."""
    if os.path.exists(output_unphased_path) and os.path.exists(output_phased_chr_path):
        return
    with pysam.VariantFile(input_phased_path) as in_vcf:
        in_sample_name = list(in_vcf.header.samples)[0]
        header = pysam.VariantHeader()
        header.add_sample("HG001")
        for r in in_vcf.header.records:
            if r.key == 'FORMAT' and r.get('ID') == 'PS':
                header.add_line('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">')
            elif r.key in ['fileformat', 'FILTER', 'FORMAT', 'INFO', 'contig']:
                header.add_record(r)
        
        with pysam.VariantFile(output_unphased_path, "w", header=header) as out_unphased, \
             pysam.VariantFile(output_phased_chr_path, "w", header=header) as out_phased:
            for record in in_vcf.fetch(chromosome):
                # VCF Phased
                rec_phased = header.new_record(
                    contig=record.chrom, start=record.start, stop=record.stop,
                    alleles=record.alleles, id=record.id, qual=record.qual, filter=record.filter.keys()
                )
                in_sample = record.samples[in_sample_name]
                out_sample = rec_phased.samples['HG001']
                for fmt_key in header.formats:
                    if fmt_key == 'PS':
                        ps_val = in_sample.get('PS')
                        if ps_val is not None:
                            try:
                                out_sample['PS'] = int(ps_val)
                            except ValueError:
                                out_sample['PS'] = 1
                    else:
                        if fmt_key in in_sample:
                            out_sample[fmt_key] = in_sample[fmt_key]
                out_sample.phased = in_sample.phased
                out_phased.write(rec_phased)

                # VCF Unphased
                rec_unphased = header.new_record(
                    contig=record.chrom, start=record.start, stop=record.stop,
                    alleles=record.alleles, id=record.id, qual=record.qual, filter=record.filter.keys()
                )
                in_sample = record.samples[in_sample_name]
                out_sample = rec_unphased.samples['HG001']
                for fmt_key in header.formats:
                    if fmt_key != 'PS':
                        if fmt_key in in_sample:
                            out_sample[fmt_key] = in_sample[fmt_key]
                out_sample.phased = False
                out_unphased.write(rec_unphased)

def load_variants_dict(vcf_path, chromosome):
    """Charge les variants hétérozygotes bialléliques du VCF."""
    variants = {}
    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf:
            if record.chrom != chromosome:
                continue
            sample = record.samples[0]
            if sample.allele_indices[0] != sample.allele_indices[1]:
                if len(record.ref) == 1 and len(record.alts[0]) == 1:
                    ref_allele = record.ref.upper()
                    alt_allele = record.alts[0].upper()
                    if ref_allele in ('A', 'C', 'G', 'T') and alt_allele in ('A', 'C', 'G', 'T'):
                        variants[record.pos - 1] = {
                            "ref": record.ref,
                            "alt": record.alts[0],
                            "gt": sample.allele_indices
                        }
    return variants

def build_hic_bam_from_pairs(pairs_path, bam_path, variants, max_pairs=10000000000000, error_rate=0.01):
    """Génère un BAM Hi-C contenant des reads de 1 pb alignés sur les variants."""
    header = {
        'HD': {'VN': '1.0', 'SO': 'unsorted'},
        'SQ': [{'SN': 'chr22', 'LN': 50818468}]
    }
    import gzip
    temp_bam_path = bam_path + ".temp.bam"
    v_positions = sorted(variants.keys())
    V = len(v_positions)
    if V == 0:
        raise ValueError("Aucun variant trouvé dans le VCF pour construire le BAM.")
    
    print("⏳ Génération des reads d'alignement à partir du fichier de paires ENCODE...")
    with pysam.AlignmentFile(temp_bam_path, "wb", header=header) as out_bam:
        pair_idx = 0
        with gzip.open(pairs_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                chrom1, pos1_str, chrom2, pos2_str = parts[1], parts[2], parts[3], parts[4]
                if chrom1 != 'chr22' or chrom2 != 'chr22':
                    continue
                pos1 = int(pos1_str) - 1
                pos2 = int(pos2_str) - 1
                
                hap = random.choice([0, 1])
                
                # Read 1 : Trouver le premier variant couvert
                start_idx1 = bisect.bisect_left(v_positions, pos1)
                var1_pos = -1
                base1 = "N"
                for idx in range(start_idx1, V):
                    v_pos = v_positions[idx]
                    if v_pos >= pos1 + 100:
                        break
                    var1_pos = v_pos
                    v_allele = variants[v_pos]['gt'][hap]
                    if random.random() < error_rate:
                        v_allele = 1 - v_allele
                    base1 = variants[v_pos]['ref'] if v_allele == 0 else variants[v_pos]['alt']
                    break
                
                if var1_pos == -1:
                    var1_pos = pos1
                    base1 = "N"
                    
                # Read 2 : Trouver le premier variant couvert
                start_idx2 = bisect.bisect_left(v_positions, pos2)
                var2_pos = -1
                base2 = "N"
                for idx in range(start_idx2, V):
                    v_pos = v_positions[idx]
                    if v_pos >= pos2 + 100:
                        break
                    var2_pos = v_pos
                    v_allele = variants[v_pos]['gt'][hap]
                    if random.random() < error_rate:
                        v_allele = 1 - v_allele
                    base2 = variants[v_pos]['ref'] if v_allele == 0 else variants[v_pos]['alt']
                    break
                
                if var2_pos == -1:
                    var2_pos = pos2
                    base2 = "N"
                
                if base1 == "N" and base2 == "N":
                    continue
                if base1 == '*' or base2 == '*':
                    continue
                
                read_name = f"hic_pair_{pair_idx}"
                
                # Écriture du Read 1 (1 pb)
                a = pysam.AlignedSegment(out_bam.header)
                a.query_name = read_name
                a.query_sequence = base1
                a.query_qualities = pysam.qualitystring_to_array("I")
                a.reference_id = 0
                a.reference_start = var1_pos
                a.cigar = ((0, 1),)
                a.mapping_quality = 60
                a.is_paired = True
                a.is_read1 = True
                a.is_read2 = False
                a.is_unmapped = False
                a.mate_is_unmapped = False
                a.next_reference_id = 0
                a.mate_is_reverse = False
                a.is_reverse = False
                a.next_reference_start = var2_pos
                a.template_length = var2_pos - var1_pos + 1
                
                # Écriture du Read 2 (1 pb)
                b = pysam.AlignedSegment(out_bam.header)
                b.query_name = read_name
                b.query_sequence = base2
                b.query_qualities = pysam.qualitystring_to_array("I")
                b.reference_id = 0
                b.reference_start = var2_pos
                b.cigar = ((0, 1),)
                b.mapping_quality = 60
                b.is_paired = True
                b.is_read1 = False
                b.is_read2 = True
                b.is_unmapped = False
                b.mate_is_unmapped = False
                b.next_reference_id = 0
                b.mate_is_reverse = False
                b.is_reverse = False
                b.next_reference_start = var1_pos
                b.template_length = -(var2_pos - var1_pos + 1)
                
                out_bam.write(a)
                out_bam.write(b)
                
                pair_idx += 1
                if pair_idx >= max_pairs:
                    break
                    
    pysam.samtools.sort(temp_bam_path, "-o", bam_path)
    pysam.samtools.index(bam_path)
    if os.path.exists(temp_bam_path):
        os.remove(temp_bam_path)
    print(f"✅ BAM Hi-C réel généré avec succès : {bam_path} ({pair_idx} paires d'interactions).")

def extract_read_profiles_fast(bam_path, variants, chromosome, mapq_threshold=20):
    """
    Extrait les profils de reads directement du BAM avec pysam de manière ultra-rapide,
    sans passer par l'alignement local PairHMM lent de WhatsHap.
    """
    print("⏳ [Extraction Rapide] Lecture séquentielle du BAM et extraction des allèles...")
    read_profiles = defaultdict(dict)
    
    with pysam.AlignmentFile(bam_path, "rb") as samfile:
        for read in samfile.fetch(chromosome):
            if read.is_unmapped or read.mapping_quality < mapq_threshold:
                continue
            
            # Pour chaque position alignée
            ref_positions = read.get_reference_positions(full_length=True)
            seq = read.query_sequence
            quals = read.query_qualities
            
            if seq is None:
                continue
                
            for i, ref_pos in enumerate(ref_positions):
                if ref_pos is None or ref_pos not in variants:
                    continue
                
                base = seq[i]
                qual = quals[i] if quals is not None else 30
                # Calcul du taux d'erreur à partir de la qualité Phred
                eps = 10 ** (-qual / 10.0)
                eps = max(1e-6, min(eps, 0.5))
                
                # Allèle numérique : 0 si REF, 1 si ALT, -1 sinon
                if base == variants[ref_pos]["ref"]:
                    allele_val = 0
                elif base == variants[ref_pos]["alt"]:
                    allele_val = 1
                else:
                    allele_val = -1
                    
                if allele_val != -1:
                    read_profiles[read.query_name][ref_pos] = (allele_val, eps)
                    
    return dict(read_profiles)
