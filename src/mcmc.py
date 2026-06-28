import numpy as np
import numba
import random

@numba.njit
def build_cross_arrays_numba(R, left, right, weight):
    cut_counts = np.zeros(R, dtype=np.int32)
    for idx in range(len(left)):
        u = left[idx]
        v = right[idx]
        for q in range(u, v):
            cut_counts[q] += 1

    cross_offsets = np.zeros(R, dtype=np.int32)
    for q in range(R - 1):
        cross_offsets[q+1] = cross_offsets[q] + cut_counts[q]

    total_crossings = cross_offsets[R-1]
    cross_left = np.empty(total_crossings, dtype=np.int32)
    cross_right = np.empty(total_crossings, dtype=np.int32)
    cross_weight = np.empty(total_crossings, dtype=np.float64)

    current_idx = cross_offsets.copy()
    for idx in range(len(left)):
        u = left[idx]
        v = right[idx]
        w = weight[idx]
        for q in range(u, v):
            write_pos = current_idx[q]
            cross_left[write_pos] = u
            cross_right[write_pos] = v
            cross_weight[write_pos] = w
            current_idx[q] += 1

    return cross_offsets, cross_left, cross_right, cross_weight

@numba.njit
def build_incident_arrays_numba(R, left, right, weight):
    node_counts = np.zeros(R + 1, dtype=np.int32)
    for idx in range(len(left)):
        u = left[idx]
        v = right[idx]
        node_counts[u] += 1
        node_counts[v] += 1

    incident_offsets = np.zeros(R + 1, dtype=np.int32)
    for r in range(R):
        incident_offsets[r+1] = incident_offsets[r] + node_counts[r]

    total_incident = incident_offsets[R]
    incident_left = np.empty(total_incident, dtype=np.int32)
    incident_right = np.empty(total_incident, dtype=np.int32)
    incident_weight = np.empty(total_incident, dtype=np.float64)

    current_idx = incident_offsets.copy()
    for idx in range(len(left)):
        u = left[idx]
        v = right[idx]
        w = weight[idx]

        pos_u = current_idx[u]
        incident_left[pos_u] = u
        incident_right[pos_u] = v
        incident_weight[pos_u] = w
        current_idx[u] += 1

        pos_v = current_idx[v]
        incident_left[pos_v] = u
        incident_right[pos_v] = v
        incident_weight[pos_v] = w
        current_idx[v] += 1

    return incident_offsets, incident_left, incident_right, incident_weight

def build_structures_fast(R, W):
    """
    Construit les tableaux d'adjacence plats requis par la recherche k-hop et le Gibbs MCMC.
    W est la matrice d'adjacence sparse (CSR) issue de build_signed_graph.
    """
    coo = W.tocoo()
    left = np.minimum(coo.row, coo.col).astype(np.int32)
    right = np.maximum(coo.row, coo.col).astype(np.int32)
    weight = coo.data.astype(np.float64)

    cross_offsets, cross_left, cross_right, cross_weight = build_cross_arrays_numba(R, left, right, weight)
    incident_offsets, incident_left, incident_right, incident_weight = build_incident_arrays_numba(R, left, right, weight)

    return (cross_offsets, cross_left, cross_right, cross_weight,
            incident_offsets, incident_left, incident_right, incident_weight)

@numba.njit
def get_pairs_bfs_counts_and_fill(R, incident_offsets, incident_left, incident_right, k_hop, fill_mode, out_left, out_right):
    visited = np.zeros(R, dtype=np.int32)
    queue = np.empty(R, dtype=np.int32)
    depth = np.empty(R, dtype=np.int32)

    total_pairs = 0
    for start_node in range(R):
        head = 0
        tail = 0

        queue[tail] = start_node
        depth[tail] = 0
        tail += 1

        visited[start_node] = start_node + 1

        while head < tail:
            u = queue[head]
            d = depth[head]
            head += 1

            if d > 0:
                if fill_mode == 1:
                    out_left[total_pairs] = start_node
                    out_right[total_pairs] = u
                total_pairs += 1

            if d < k_hop:
                start_offset = incident_offsets[u]
                end_offset = incident_offsets[u+1]
                for idx in range(start_offset, end_offset):
                    v = incident_right[idx] if incident_left[idx] == u else incident_left[idx]
                    if visited[v] != start_node + 1:
                        visited[v] = start_node + 1
                        queue[tail] = v
                        depth[tail] = d + 1
                        tail += 1

    return total_pairs

def get_pairs_bfs(R, incident_offsets, incident_left, incident_right, k_hop):
    dummy = np.empty(0, dtype=np.int32)
    total_pairs = get_pairs_bfs_counts_and_fill(
        R, incident_offsets, incident_left, incident_right, k_hop, 0, dummy, dummy
    )
    
    out_left = np.empty(total_pairs, dtype=np.int32)
    out_right = np.empty(total_pairs, dtype=np.int32)
    
    get_pairs_bfs_counts_and_fill(
        R, incident_offsets, incident_left, incident_right, k_hop, 1, out_left, out_right
    )
    
    return out_left, out_right

@numba.njit
def mcmc_gibbs_sampling_numba(R, incident_offsets, incident_left, incident_right, incident_weight,
                              initial_spins, steps, beta):
    """
    Simule la MCMC avec échantillonnage de Gibbs local sur les spins des reads.
    Retourne la probabilité marginale de spin de chaque read.
    """
    spins = initial_spins.copy()
    spin_sums = np.zeros(R, dtype=np.float64)
    
    for step in range(steps):
        u = random.randint(0, R - 1)
        
        start = incident_offsets[u]
        end = incident_offsets[u+1]
        local_field = 0.0
        
        for idx in range(start, end):
            v = incident_right[idx] if incident_left[idx] == u else incident_left[idx]
            w = incident_weight[idx]
            local_field += w * spins[v]
            
        p = 1.0 / (1.0 + np.exp(-2.0 * beta * local_field))
        
        if random.random() < p:
            spins[u] = 1.0
        else:
            spins[u] = -1.0
            
        if step >= steps // 2:
            spin_sums += spins
            
    n_samples = steps - (steps // 2)
    mean_spins = spin_sums / n_samples
    post_probs = (mean_spins + 1.0) / 2.0
    
    return post_probs
