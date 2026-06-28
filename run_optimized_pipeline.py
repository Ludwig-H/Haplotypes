import os
import sys
import numpy as np
import time

# S'assurer que le répertoire racine est dans le path
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

import src

def main():
    chromosome = "chr22"
    sequencing_tech = "Hi-C"
    k_hop = 2
    beta_safety = 0.5
    mcmc_steps = 100000
    beta = 1.0
    
    # 1. Préparation des données
    vcf_phased = f"HG001_phased_{chromosome}.vcf"
    vcf_unphased = f"HG001_unphased_{chromosome}.vcf"
    bam_path = f"NA12878_{chromosome}_hic.bam"
    
    print("⏳ [Pipeline Optimisé] Démarrage du pipeline de phasing...")
    
    if not os.path.exists(vcf_phased) or not os.path.exists(vcf_unphased):
        print("🧬 Création des fichiers VCF locaux...")
        src.prepare_cleaned_vcfs("HG001_phased.vcf.gz", vcf_unphased, vcf_phased, chromosome)
        
    print("🧬 Chargement des variants de vérité terrain...")
    variants = src.load_variants_dict(vcf_phased, chromosome)
    print(f"  -> {len(variants)} variants hétérozygotes chargés.")
    
    # Si le BAM Hi-C n'existe pas, on le télécharge/génère
    if not os.path.exists(bam_path):
        PAIRS_URL = "https://www.encodeproject.org/files/ENCFF527IWE/@@download/ENCFF527IWE.pairs.gz"
        pairs_local = "ENCFF527IWE.pairs.gz"
        pairs_filtered = f"ENCFF527IWE_{chromosome}.pairs.gz"
        
        src.download_file(PAIRS_URL, pairs_local)
        if not os.path.exists(pairs_filtered):
            print(f"⏳ Extraction ultra-rapide des contacts {chromosome}...")
            os.system(f"zcat {pairs_local} | grep -E '^#|{chromosome}' | gzip > {pairs_filtered}")
            if os.path.exists(pairs_local):
                os.remove(pairs_local)
                
        # Générer le BAM ultra-rapide (reads de 1 pb)
        src.build_hic_bam_from_pairs(pairs_filtered, bam_path, variants, max_pairs=10000000000000)

    # 2. Extraction des profils de reads
    t0 = time.time()
    read_profiles = src.extract_read_profiles_fast(bam_path, variants, chromosome)
    print(f"  -> {len(read_profiles)} profils de reads extraits en {time.time() - t0:.2f}s.")
    
    # 3. Détermination des spins de vérité terrain pour l'évaluation
    true_spins = {}
    for rid, profile in read_profiles.items():
        votes = []
        for pos, (allele_val, _) in profile.items():
            if pos in variants:
                true_gt = variants[pos]["gt"]
                if allele_val == true_gt[0]:
                    votes.append(1)
                elif allele_val == true_gt[1]:
                    votes.append(-1)
        if len(votes) > 0:
            true_spins[rid] = 1 if np.sum(votes) >= 0 else -1
            
    active_reads = {rid: prof for rid, prof in read_profiles.items() if rid in true_spins}
    read_ids = sorted(active_reads.keys(), key=lambda rid: min(active_reads[rid].keys()))
    R = len(read_ids)
    read_id_to_idx = {rid: i for i, rid in enumerate(read_ids)}
    true_spin_vec = np.array([true_spins[rid] for rid in read_ids])
    
    print(f"  -> {R} reads actifs restants après filtrage de vérité terrain.")
    
    # 4. Construction du graphe signé avec Numba et SciPy sparse
    t0 = time.time()
    W = src.build_signed_graph(active_reads, variants, read_id_to_idx, beta_safety=beta_safety, max_reads_per_variant=150)
    print(f"✅ Graphe construit (CSR sparse) de taille {W.shape[0]}x{W.shape[1]} en {time.time() - t0:.2f}s.")
    
    # 5. Clustering Spectral Signé avec Fallback GPU (CuPy/SciPy)
    t0 = time.time()
    # Nous désactivons l'utilisation du GPU ici dans le test s'il n'y a pas de GPU CUDA dans Codespace
    v_eigen = src.signed_spectral_clustering(W, type="unnormalized", use_gpu=False)
    pred_baseline = np.sign(v_eigen)
    print(f"✅ Clustering Spectral Signé terminé en {time.time() - t0:.2f}s.")
    
    # Aligner le baseline avec la vérité terrain pour mesurer la précision brute
    corr = np.mean(pred_baseline == true_spin_vec)
    if corr < 0.5:
        pred_baseline = -pred_baseline
        corr = 1 - corr
    print(f"  -> Précision du clustering spectral baseline : {corr:.2%}")
    
    # 6. MCMC k-hop par Fenwick et second Clustering Spectral
    t0 = time.time()
    structures = src.build_structures_fast(R, W)
    cross_offsets, cross_left, cross_right, cross_weight = structures[0], structures[1], structures[2], structures[3]
    incident_offsets = structures[4]
    incident_left = structures[5]
    incident_right = structures[6]
    incident_weight = structures[7]
    
    print(f"⏳ Génération des paires {k_hop}-hop...")
    pairs_left, pairs_right = src.build_pairs_sparse(R, W, k_hop, verbose=False)
    P = len(pairs_left)
    print(f"  -> Nombre de paires {k_hop}-hop : {P}")
    
    print("⏳ Initialisation de la MCMC avec les spins de la Baseline spectrale...")
    tree = np.zeros(R, dtype=np.int32)
    tau_init = (pred_baseline[:-1] != pred_baseline[1:]).astype(np.int32)
    for k in range(R - 1):
        if tau_init[k] == 1:
            src.fenwick_update(tree, k, 1)
            
    print(f"⏳ Lancement de la boucle MCMC ({mcmc_steps} étapes)...")
    correlations = src.mcmc_loop(
        mcmc_steps, R, beta, tree,
        cross_offsets, cross_left, cross_right, cross_weight,
        incident_offsets, incident_left, incident_right, incident_weight,
        pairs_left, pairs_right, verbose=False
    )
    
    print("⏳ Résolution spectrale sur la matrice de corrélation...")
    import scipy.sparse as sp
    row_C = np.concatenate([pairs_left, pairs_right])
    col_C = np.concatenate([pairs_right, pairs_left])
    data_C = np.concatenate([correlations, correlations])
    W_C = sp.coo_matrix((data_C, (row_C, col_C)), shape=(R, R)).tocsr()
    
    v_C = src.signed_spectral_clustering(W_C, type="unnormalized", use_gpu=False)
    pred_mcmc = np.sign(v_C)
    pred_mcmc[pred_mcmc == 0] = 1
    print(f"✅ MCMC et second clustering terminés en {time.time() - t0:.2f}s.")
    
    corr_mcmc = np.mean(pred_mcmc == true_spin_vec)
    if corr_mcmc < 0.5:
        pred_mcmc = -pred_mcmc
        corr_mcmc = 1 - corr_mcmc
    print(f"  -> Précision de la MCMC k-hop : {corr_mcmc:.2%}")
    
    # 7. Écriture du VCF phasé final
    # Construire un dictionnaire variant_ps pour les PS blocks (simplifié)
    variant_ps = {}
    for pos in variants.keys():
        variant_ps[pos] = 1 # Un seul bloc pour tout le chromosome
        
    print("⏳ Exportation du VCF phasé...")
    src.write_phased_vcf(vcf_unphased, "phased_mcmc_fast.vcf", active_reads, read_id_to_idx, pred_mcmc, variants, variant_ps, chromosome)
    print("✅ VCF 'phased_mcmc_fast.vcf' exporté avec succès.")
    
    # 8. Évaluation des métriques
    comp = src.parse_whatshap_compare("", "phased_mcmc_fast.vcf", chromosome)
    print("\n📊 RÉSULTATS PIPELINE OPTIMISÉ :")
    print(f"  - Variants phasés : {comp['phased_variants']}")
    print(f"  - N50 de bloc : {comp['block_n50']} pb")

if __name__ == "__main__":
    main()
