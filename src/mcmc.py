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
    Construit les structures compressées requises pour la recherche k-hop et la MCMC par Fenwick.
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

def build_pairs_sparse(R, W, k, verbose=False):
    """
    Génère les paires k-hop à partir de la matrice d'adjacence sparse W.
    """
    if verbose:
        print(f"💻 [build_pairs_sparse] Génération des paires {k}-hop via BFS Numba...")
        
    coo = W.tocoo()
    left = np.minimum(coo.row, coo.col).astype(np.int32)
    right = np.maximum(coo.row, coo.col).astype(np.int32)
    weight = coo.data.astype(np.float64)
    
    incident_offsets, incident_left, incident_right, incident_weight = build_incident_arrays_numba(R, left, right, weight)
    out_left, out_right = get_pairs_bfs(R, incident_offsets, incident_left, incident_right, k)
    return out_left, out_right

@numba.njit
def fenwick_update(tree, idx, val):
    i = idx + 1
    n = len(tree)
    while i < n:
        tree[i] ^= val
        i += i & (-i)

@numba.njit
def fenwick_query(tree, idx):
    if idx < 0:
        return 0
    i = idx + 1
    res = 0
    while i > 0:
        res ^= tree[i]
        i -= i & (-i)
    return res

@numba.njit
def get_cached_query(tree, idx, cache_vals, cache_steps, t):
    if idx < 0:
        return 0
    if cache_steps[idx] == t:
        return cache_vals[idx]
    i = idx + 1
    res = 0
    while i > 0:
        res ^= tree[i]
        i -= i & (-i)
    cache_vals[idx] = res
    cache_steps[idx] = t
    return res

@numba.njit
def evaluate_cut(q, tree, cross_offsets, cross_left, cross_right, cross_weight, cache_vals, cache_steps, t):
    start = cross_offsets[q]
    end = cross_offsets[q+1]
    du = 0.0
    for idx in range(start, end):
        i = cross_left[idx]
        j = cross_right[idx]
        w = cross_weight[idx]
        q_j = get_cached_query(tree, j-1, cache_vals, cache_steps, t)
        q_i = get_cached_query(tree, i-1, cache_vals, cache_steps, t)
        xor_val = q_j ^ q_i
        spin_prod = 1.0 - 2.0 * float(xor_val)
        du += w * spin_prod
    return du

@numba.njit
def evaluate_incident(r, tree, incident_offsets, incident_left, incident_right, incident_weight, cache_vals, cache_steps, t):
    start = incident_offsets[r]
    end = incident_offsets[r+1]
    du_L = 0.0
    du_R = 0.0
    for idx in range(start, end):
        i = incident_left[idx]
        j = incident_right[idx]
        w = incident_weight[idx]
        q_j = get_cached_query(tree, j-1, cache_vals, cache_steps, t)
        q_i = get_cached_query(tree, i-1, cache_vals, cache_steps, t)
        xor_val = q_j ^ q_i
        spin_prod = 1.0 - 2.0 * float(xor_val)
        if i == r:
            du_R += w * spin_prod
        else:
            du_L += w * spin_prod
    return du_L, du_R

@numba.njit
def mcmc_loop(steps, R, beta, tree,
              cross_offsets, cross_left, cross_right, cross_weight,
              incident_offsets, incident_left, incident_right, incident_weight,
              pairs_left, pairs_right, verbose=False):
    """
    Exécute la boucle MCMC sur l'arbre de Fenwick et calcule les corrélations temporelles des paires k-hop.
    """
    tree_init = tree.copy()

    max_flips = 2 * steps
    flip_steps = np.empty(max_flips, dtype=np.int32)
    flip_walls = np.empty(max_flips, dtype=np.int32)
    flip_counts = np.zeros(R - 1, dtype=np.int32)
    total_flips = 0

    cache_vals = np.empty(R, dtype=np.int32)
    cache_steps = np.zeros(R, dtype=np.int32)

    for t in range(1, steps + 1):
        r = (t - 1) if t <= R else random.randint(0, R - 1)
        du0 = 0.0

        du_L, du_R = evaluate_incident(r, tree, incident_offsets, incident_left, incident_right, incident_weight, cache_vals, cache_steps, t)
        du1 = du_L + du_R

        if r - 1 >= 0:
            du2 = evaluate_cut(r-1, tree, cross_offsets, cross_left, cross_right, cross_weight, cache_vals, cache_steps, t)
        else:
            du2 = 0.0

        if r < R - 1:
            du3 = du2 - du_L + du_R
        else:
            du3 = 0.0

        min_du = min(du0, du1, du2, du3)
        w0 = np.exp(-beta * (du0 - min_du))
        w1 = np.exp(-beta * (du1 - min_du))
        w2 = np.exp(-beta * (du2 - min_du))
        w3 = np.exp(-beta * (du3 - min_du))
        sum_w = w0 + w1 + w2 + w3
        p0 = w0 / sum_w
        p1 = w1 / sum_w
        p2 = w2 / sum_w
        
        rand_val = random.random()
        chosen_move = 0
        if rand_val < p0:
            chosen_move = 0
        elif rand_val < p0 + p1:
            chosen_move = 1
        elif rand_val < p0 + p1 + p2:
            chosen_move = 2
        else:
            chosen_move = 3

        if chosen_move == 1:
            if r - 1 >= 0:
                fenwick_update(tree, r-1, 1)
                flip_steps[total_flips] = t
                flip_walls[total_flips] = r-1
                total_flips += 1
                flip_counts[r-1] += 1
            if r < R - 1:
                fenwick_update(tree, r, 1)
                flip_steps[total_flips] = t
                flip_walls[total_flips] = r
                total_flips += 1
                flip_counts[r] += 1
        elif chosen_move == 2:
            if r - 1 >= 0:
                fenwick_update(tree, r-1, 1)
                flip_steps[total_flips] = t
                flip_walls[total_flips] = r-1
                total_flips += 1
                flip_counts[r-1] += 1
        elif chosen_move == 3:
            if r < R - 1:
                fenwick_update(tree, r, 1)
                flip_steps[total_flips] = t
                flip_walls[total_flips] = r
                total_flips += 1
                flip_counts[r] += 1

    actual_flip_steps = flip_steps[:total_flips]
    actual_flip_walls = flip_walls[:total_flips]
    offsets = np.zeros(R, dtype=np.int32)
    for k in range(R - 1):
        offsets[k+1] = offsets[k] + flip_counts[k]
        
    grouped_steps = np.empty(total_flips, dtype=np.int32)
    current_idx = offsets.copy()
    for i in range(total_flips):
        k = actual_flip_walls[i]
        t = actual_flip_steps[i]
        pos = current_idx[k]
        grouped_steps[pos] = t
        current_idx[k] += 1

    P = len(pairs_left)
    correlations = np.empty(P, dtype=np.float64)
    for p in range(P):
        u = pairs_left[p]
        v = pairs_right[p]
        sum_counts = 0
        for k in range(u, v):
            sum_counts += flip_counts[k]

        xor_init = fenwick_query(tree_init, v-1) ^ fenwick_query(tree_init, u-1)
        initial_sign = 1.0 - 2.0 * float(xor_init)

        if sum_counts == 0:
            correlations[p] = initial_sign
            continue
            
        temp = np.empty(sum_counts, dtype=np.int32)
        idx_temp = 0
        for k in range(u, v):
            start = offsets[k]
            count = flip_counts[k]
            for idx_grouped in range(start, start + count):
                temp[idx_temp] = grouped_steps[idx_grouped]
                idx_temp += 1
                
        temp = np.sort(temp)
        sum_prod = 0.0
        current_sign = initial_sign
        last_t = 1
        idx = 0
        while idx < sum_counts:
            t_val = temp[idx]
            cnt = 1
            while idx + 1 < sum_counts and temp[idx + 1] == t_val:
                cnt += 1
                idx += 1
            if cnt % 2 == 1:
                sum_prod += current_sign * (t_val - last_t)
                current_sign = -current_sign
                last_t = t_val
            idx += 1
        sum_prod += current_sign * (steps + 1 - last_t)
        correlations[p] = sum_prod / float(steps)
        
    return correlations

@numba.njit
def mcmc_gibbs_sampling_numba(R, incident_offsets, incident_left, incident_right, incident_weight,
                              initial_spins, steps, beta):
    """
    Simule la MCMC avec échantillonnage de Gibbs local sur les spins des reads.
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
