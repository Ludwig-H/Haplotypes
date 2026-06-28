import numpy as np
import scipy.sparse as sp
import numba
import random
from collections import defaultdict

@numba.njit
def build_edges_raw_numba(rids_per_variant_flat, variant_offsets, read_values, read_errors, max_reads):
    """
    Construit les arêtes du graphe signé de manière ultra-rapide et vectorisée avec Numba.
    Retourne des tableaux plats d'indices sources, cibles et de poids de contribution.
    """
    V = len(variant_offsets) - 1
    estimated_size = 0
    for v in range(V):
        start = variant_offsets[v]
        end = variant_offsets[v+1]
        n_rids = end - start
        if n_rids > max_reads:
            n_rids = max_reads
        estimated_size += (n_rids * (n_rids - 1)) // 2
        
    sources = np.empty(estimated_size, dtype=np.int32)
    targets = np.empty(estimated_size, dtype=np.int32)
    weights = np.empty(estimated_size, dtype=np.float64)
    
    edge_idx = 0
    for v in range(V):
        start = variant_offsets[v]
        end = variant_offsets[v+1]
        n_rids = end - start
        
        # Capping local des reads pour éviter les OOM quadratiques
        if n_rids > max_reads:
            n_rids = max_reads
            
        for i in range(n_rids):
            u = rids_per_variant_flat[start + i]
            val_i = read_values[start + i]
            eps_i = read_errors[start + i]
            
            for j in range(i + 1, n_rids):
                v_idx = rids_per_variant_flat[start + j]
                val_j = read_values[start + j]
                eps_j = read_errors[start + j]
                
                # Calcul de la probabilité d'accord sous le même haplotype
                q_z = (1.0 - eps_i) * (1.0 - eps_j) + eps_i * eps_j
                if q_z < 1e-6:
                    q_z = 1e-6
                elif q_z > 1.0 - 1e-6:
                    q_z = 1.0 - 1e-6
                log_ratio = np.log(q_z / (1.0 - q_z))
                
                is_concordant = (val_i == val_j)
                w = log_ratio if is_concordant else -log_ratio
                
                if u < v_idx:
                    sources[edge_idx] = u
                    targets[edge_idx] = v_idx
                else:
                    sources[edge_idx] = v_idx
                    targets[edge_idx] = u
                weights[edge_idx] = w
                edge_idx += 1
                
    return sources[:edge_idx], targets[:edge_idx], weights[:edge_idx]

def build_signed_graph(read_profiles, variants, read_id_to_idx, beta_safety=0.5, max_reads_per_variant=150):
    """
    Construit la matrice d'adjacence signée et pondérée W du graphe de reads.
    Gère le capage de couverture et l'accumulation des doublons de manière optimisée.
    """
    R = len(read_id_to_idx)
    
    # Étape 1 : Grouper les reads par variant
    reads_per_variant = defaultdict(list)
    for rid, profile in read_profiles.items():
        if rid not in read_id_to_idx:
            continue
        for pos in profile.keys():
            if pos in variants:
                reads_per_variant[pos].append(rid)
                
    # Étape 2 : Préparer les tableaux plats pour Numba
    V = len(reads_per_variant)
    variant_offsets = np.zeros(V + 1, dtype=np.int32)
    
    # Mélange reproductible et capage
    rids_list = []
    values_list = []
    errors_list = []
    
    offset = 0
    for idx, (pos, rids) in enumerate(reads_per_variant.items()):
        if len(rids) > max_reads_per_variant:
            random.seed(42)
            random.shuffle(rids)
            rids = rids[:max_reads_per_variant]
            
        variant_offsets[idx] = offset
        for rid in rids:
            # On stocke directement l'index numérique du read pour éviter de manipuler des strings dans Numba
            rids_list.append(read_id_to_idx[rid])
            val, eps = read_profiles[rid][pos]
            values_list.append(val)
            errors_list.append(eps)
        offset += len(rids)
    variant_offsets[V] = offset
    
    # Tableaux NumPy pour Numba
    rids_flat = np.array(rids_list, dtype=np.int32)
    values_flat = np.array(values_list, dtype=np.int32)
    errors_flat = np.array(errors_list, dtype=np.float64)
    
    # Étape 3 : Appeler le kernel Numba pour construire les arêtes brutes
    sources, targets, weights = build_edges_raw_numba(
        rids_flat, variant_offsets, values_flat, errors_flat, max_reads_per_variant
    )
    
    # Appliquer le facteur de sécurité beta
    weights = weights * beta_safety
    
    # Étape 4 : Regrouper et sommer les doublons en C avec SciPy sparse
    W_raw = sp.coo_matrix((weights, (sources, targets)), shape=(R, R))
    W_csr = W_raw.tocsr()
    
    return W_csr
