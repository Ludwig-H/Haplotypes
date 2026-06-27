import pandas as pd
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import numba
import time
import os

# Vérification de la disponibilité du GPU
GPU_AVAILABLE = False
try:
    import cupy as cp
    import cupyx.scipy.sparse as csp
    from cupyx.scipy.sparse.linalg import eigsh as cupy_eigsh
    # Simple test pour s'assurer que le GPU est réellement accessible
    cp.array([1])
    GPU_AVAILABLE = True
except Exception:
    pass

def load_instance(instance_dir):
    nodes_df = pd.read_csv(f"{instance_dir}/graph/nodes.tsv", sep="\t")
    edges_df = pd.read_csv(f"{instance_dir}/graph/edges.tsv", sep="\t")
    truth_df = pd.read_csv(f"{instance_dir}/truth/read_truth.tsv", sep="\t")
    
    R = len(nodes_df)
    node_to_truth = dict(zip(truth_df['read_id'], truth_df['true_haplotype']))
    true_spins = np.array([node_to_truth[read_id] for read_id in nodes_df['read_id']])
    
    return R, edges_df, true_spins

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

def build_structures_fast(R, edges_df):
    left = np.minimum(edges_df['source'], edges_df['target']).values.astype(np.int32)
    right = np.maximum(edges_df['source'], edges_df['target']).values.astype(np.int32)
    weight = edges_df['weight'].values.astype(np.float64)
    
    cross_offsets, cross_left, cross_right, cross_weight = build_cross_arrays_numba(R, left, right, weight)
    incident_offsets, incident_left, incident_right, incident_weight = build_incident_arrays_numba(R, left, right, weight)
    
    return (cross_offsets, cross_left, cross_right, cross_weight,
            incident_offsets, incident_left, incident_right, incident_weight)

def build_pairs_sparse(R, edges_df, k):
    left = np.minimum(edges_df['source'], edges_df['target']).values.astype(np.int32)
    right = np.maximum(edges_df['source'], edges_df['target']).values.astype(np.int32)
    
    @numba.njit
    def build_adj_list(R, left, right):
        degrees = np.zeros(R + 1, dtype=np.int32)
        for i in range(len(left)):
            degrees[left[i]] += 1
            degrees[right[i]] += 1
        offsets = np.zeros(R + 1, dtype=np.int32)
        for i in range(R):
            offsets[i+1] = offsets[i] + degrees[i]
        adj = np.empty(offsets[R], dtype=np.int32)
        current = offsets.copy()
        for i in range(len(left)):
            u = left[i]
            v = right[i]
            adj[current[u]] = v
            current[u] += 1
            adj[current[v]] = u
            current[v] += 1
        return offsets, adj

    offsets, adj = build_adj_list(R, left, right)
    
    @numba.njit
    def bfs_k_hop(R, offsets, adj, k):
        max_pairs = R * 100
        p_left = np.empty(max_pairs, dtype=np.int32)
        p_right = np.empty(max_pairs, dtype=np.int32)
        p_count = 0
        
        visited = np.zeros(R, dtype=np.int32)
        queue = np.empty(R, dtype=np.int32)
        depth = np.empty(R, dtype=np.int32)
        
        for root in range(R):
            q_head = 0
            q_tail = 0
            
            queue[q_tail] = root
            depth[root] = 0
            visited[root] = root + 1
            q_tail += 1
            
            while q_head < q_tail:
                u = queue[q_head]
                q_head += 1
                curr_depth = depth[u]
                
                if curr_depth > 0:
                    if root < u:
                        if p_count >= len(p_left):
                            new_left = np.empty(len(p_left) * 2, dtype=np.int32)
                            new_right = np.empty(len(p_right) * 2, dtype=np.int32)
                            new_left[:p_count] = p_left[:p_count]
                            new_right[:p_count] = p_right[:p_count]
                            p_left = new_left
                            p_right = new_right
                        p_left[p_count] = root
                        p_right[p_count] = u
                        p_count += 1
                        
                if curr_depth < k:
                    start = offsets[u]
                    end = offsets[u+1]
                    for idx in range(start, end):
                        v = adj[idx]
                        if visited[v] != root + 1:
                            visited[v] = root + 1
                            depth[v] = curr_depth + 1
                            queue[q_tail] = v
                            q_tail += 1
                            
        return p_left[:p_count], p_right[:p_count]

    return bfs_k_hop(R, offsets, adj, k)

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
def evaluate_cut(q, tree, cross_offsets, cross_left, cross_right, cross_weight):
    start = cross_offsets[q]
    end = cross_offsets[q+1]
    du = 0.0
    for idx in range(start, end):
        i = cross_left[idx]
        j = cross_right[idx]
        w = cross_weight[idx]
        xor_val = fenwick_query(tree, j-1) ^ fenwick_query(tree, i-1)
        spin_prod = 1.0 - 2.0 * float(xor_val)
        du += w * spin_prod
    return du

@numba.njit
def evaluate_singleton(r, tree, incident_offsets, incident_left, incident_right, incident_weight):
    start = incident_offsets[r]
    end = incident_offsets[r+1]
    du = 0.0
    for idx in range(start, end):
        i = incident_left[idx]
        j = incident_right[idx]
        w = incident_weight[idx]
        xor_val = fenwick_query(tree, j-1) ^ fenwick_query(tree, i-1)
        spin_prod = 1.0 - 2.0 * float(xor_val)
        du += w * spin_prod
    return du

@numba.njit
def mcmc_loop(steps, R, beta, tree, 
              cross_offsets, cross_left, cross_right, cross_weight,
              incident_offsets, incident_left, incident_right, incident_weight,
              pairs_left, pairs_right, verbose=False):
    tree_init = tree.copy()
    
    max_flips = 2 * steps
    flip_steps = np.empty(max_flips, dtype=np.int32)
    flip_walls = np.empty(max_flips, dtype=np.int32)
    flip_counts = np.zeros(R - 1, dtype=np.int32)
    total_flips = 0
    report_interval = steps // 10
    if report_interval == 0: report_interval = 1
    for t in range(1, steps + 1):
        r = (t - 1) if t <= R else np.random.randint(0, R)
        du0 = 0.0
        du1 = evaluate_singleton(r, tree, incident_offsets, incident_left, incident_right, incident_weight)
        if r - 1 >= 0:
            du2 = evaluate_cut(r-1, tree, cross_offsets, cross_left, cross_right, cross_weight)
        else:
            du2 = 0.0
        if r < R - 1:
            du3 = evaluate_cut(r, tree, cross_offsets, cross_left, cross_right, cross_weight)
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
        p3 = w3 / sum_w
        rand_val = np.random.random()
        chosen_move = 0
        if rand_val < p0:
            chosen_move = 0
        elif rand_val < p0 + p1:
            chosen_move = 1
        elif rand_val < p0 + p1 + p2:
            chosen_move = 2
        else:
            chosen_move = 3
            
        if verbose and t % report_interval == 0:
            print("  -> Progression MCMC : ", 100 * t // steps, "%")
            
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

def solve_signed_spectral(W, gpu=True, laplacian='unnormalized', verbose=False):
    R = W.shape[0]
    if R <= 15:
        W_dense = W.toarray()
        abs_W = np.abs(W_dense)
        degrees = abs_W.sum(axis=1)
        if gpu and GPU_AVAILABLE:
            try:
                import cupy as cp
                W_dense_gpu = cp.array(W_dense)
                degrees_gpu = cp.array(degrees)
                if laplacian == 'normalized':
                    d_inv_sqrt = cp.zeros(R, dtype=cp.float32)
                    mask = degrees_gpu > 0
                    d_inv_sqrt[mask] = 1.0 / cp.sqrt(degrees_gpu[mask])
                    A_norm = cp.diag(d_inv_sqrt) @ W_dense_gpu @ cp.diag(d_inv_sqrt)
                    vals, vecs = cp.linalg.eigh(A_norm)
                    return cp.asnumpy(vecs[:, -1]), 1.0 - cp.asnumpy(vals[-1])
                else:
                    D = cp.diag(degrees_gpu)
                    L = D - W_dense_gpu
                    vals, vecs = cp.linalg.eigh(L)
                    return cp.asnumpy(vecs[:, 0]), cp.asnumpy(vals[0])
            except Exception as e:
                if verbose: print(f"⚠️ [dense GPU] Échec: {e}. Repli sur le CPU...")
        if laplacian == 'normalized':
            d_inv_sqrt = np.zeros(R, dtype=np.float32)
            mask = degrees > 0
            d_inv_sqrt[mask] = 1.0 / np.sqrt(degrees[mask])
            A_norm = np.diag(d_inv_sqrt) @ W_dense @ np.diag(d_inv_sqrt)
            vals, vecs = np.linalg.eigh(A_norm)
            return vecs[:, -1], 1.0 - vals[-1]
        else:
            D = np.diag(degrees)
            L = D - W_dense
            vals, vecs = np.linalg.eigh(L)
            return vecs[:, 0], vals[0]

    abs_W = abs(W)
    degrees = np.array(abs_W.sum(axis=1)).flatten()
    if gpu and GPU_AVAILABLE:
        try:
            W_gpu = csp.csr_matrix(W)
            if laplacian == 'normalized':
                d_inv_sqrt = cp.zeros(R, dtype=cp.float32)
                mask = cp.array(degrees) > 0
                d_inv_sqrt[mask] = 1.0 / cp.sqrt(cp.array(degrees)[mask])
                D_inv_sqrt = csp.diags(d_inv_sqrt)
                A_norm_gpu = D_inv_sqrt @ W_gpu @ D_inv_sqrt
                vals, vecs = cupy_eigsh(A_norm_gpu, k=1, which='LA')
                return cp.asnumpy(vecs[:, 0]), 1.0 - cp.asnumpy(vals[0])
            else:
                D_gpu = csp.diags(cp.array(degrees))
                L_gpu = D_gpu - W_gpu
                vals, vecs = cupy_eigsh(L_gpu, k=1, which='SA')
                return cp.asnumpy(vecs[:, 0]), cp.asnumpy(vals[0])
        except Exception as e:
            if verbose: print(f"⚠️ [GPU] Échec: {e}. Repli sur le CPU...")
            
    if laplacian == 'normalized':
        d_inv_sqrt = np.zeros(R, dtype=np.float32)
        mask = degrees > 0
        d_inv_sqrt[mask] = 1.0 / np.sqrt(degrees[mask])
        D_inv_sqrt = sp.diags(d_inv_sqrt)
        A_norm = D_inv_sqrt @ W @ D_inv_sqrt
        try:
            vals, vecs = eigsh(A_norm, k=1, which='LA', tol=1e-5)
            return vecs[:, 0], 1.0 - vals[0]
        except Exception as e:
            if verbose: print(f"⚠️ [solve_signed_spectral CPU Normalized] eigsh LA a échoué: {e}. Essai avec lobpcg...")
            try:
                from scipy.sparse.linalg import lobpcg
                I = sp.eye(R, format='csr')
                L_sym = I - A_norm
                X = np.random.normal(size=(R, 1))
                vals, vecs = lobpcg(L_sym, X, largest=False, tol=1e-5, maxiter=200)
                return vecs[:, 0], vals[0]
            except Exception as e2:
                if verbose: print(f"⚠️ [solve_signed_spectral CPU Normalized] lobpcg a échoué: {e2}. Repli sur la version dense...")
                A_norm_dense = A_norm.toarray()
                vals, vecs = np.linalg.eigh(A_norm_dense)
                return vecs[:, -1], 1.0 - vals[-1]
    else:
        D = sp.diags(degrees)
        L = D - W
        try:
            sigma = 2.0 * np.max(degrees)
            I = sp.eye(R, format='csr')
            M = sigma * I - L
            vals, vecs = eigsh(M, k=1, which='LA', tol=1e-5)
            return vecs[:, 0], sigma - vals[0]
        except Exception as e:
            if verbose: print(f"⚠️ [solve_signed_spectral CPU Unnormalized] eigsh LA a échoué: {e}. Essai avec lobpcg...")
            try:
                from scipy.sparse.linalg import lobpcg
                X = np.random.normal(size=(R, 1))
                vals, vecs = lobpcg(L, X, largest=False, tol=1e-5, maxiter=200)
                return vecs[:, 0], vals[0]
            except Exception as e2:
                if verbose: print(f"⚠️ [solve_signed_spectral CPU Unnormalized] lobpcg a échoué: {e2}. Repli sur la version dense...")
                L_dense = D.toarray() - W.toarray()
                vals, vecs = np.linalg.eigh(L_dense)
                return vecs[:, 0], vals[0]

def solve_signed_spectral_by_components(W, gpu=True, laplacian='unnormalized', verbose=False):
    from scipy.sparse.csgraph import connected_components
    R = W.shape[0]
    n_components, labels = connected_components(W, directed=False)
    v_all = np.zeros(R, dtype=np.float64)
    val_all = np.zeros(n_components, dtype=np.float64)
    for comp_idx in range(n_components):
        nodes = np.where(labels == comp_idx)[0]
        if len(nodes) < 2:
            v_all[nodes] = 1.0
            continue
        W_sub = W[nodes][:, nodes]
        v_sub, val_sub = solve_signed_spectral(W_sub, gpu=gpu, laplacian=laplacian, verbose=verbose)
        v_all[nodes] = v_sub
        val_all[comp_idx] = val_sub
    return v_all, np.mean(val_all)

def align_spins_by_components(pred, labels, n_components, true_spins):
    pred_aligned = pred.copy()
    for comp_idx in range(n_components):
        nodes = np.where(labels == comp_idx)[0]
        if len(nodes) == 0:
            continue
        local_acc = np.mean(pred[nodes] == true_spins[nodes])
        if local_acc < 0.5:
            pred_aligned[nodes] = -pred[nodes]
    return pred_aligned

def run_instance(instance_dir, steps, k_hop, beta, laplacian='unnormalized'):
    print("=" * 60)
    print(f"RUNNING BENCHMARK ON: {instance_dir}")
    print("=" * 60)
    
    R, edges_df, true_spins = load_instance(instance_dir)
    print(f"Nodes R: {R}, Edges: {len(edges_df)}")
    
    # 1. Détection des composantes connexes
    row = np.concatenate([edges_df['source'].values, edges_df['target'].values])
    col = np.concatenate([edges_df['target'].values, edges_df['source'].values])
    data = np.concatenate([edges_df['weight'].values, edges_df['weight'].values])
    W_sparse = sp.coo_matrix((data, (row, col)), shape=(R, R)).tocsr()
    
    from scipy.sparse.csgraph import connected_components
    n_components, labels = connected_components(W_sparse, directed=False)
    
    # 2. Baseline: Signed Spectral Clustering par blocs
    t0 = time.time()
    v_W, val_W = solve_signed_spectral_by_components(W_sparse, gpu=True, laplacian=laplacian, verbose=False)
    pred_baseline = np.sign(v_W)
    pred_baseline[pred_baseline == 0] = 1
    
    # Alignement par blocs
    pred_baseline_aligned = align_spins_by_components(pred_baseline, labels, n_components, true_spins)
    acc_baseline = np.mean(pred_baseline_aligned == true_spins)
    t_baseline = time.time() - t0
    
    print(f"Baseline Accuracy : {acc_baseline:.4%}")
    print(f"Baseline Time     : {t_baseline:.2f}s")
    
    # 3. Glauber MCMC
    t0 = time.time()
    cross_offsets, cross_left, cross_right, cross_weight, \
    incident_offsets, incident_left, incident_right, incident_weight = build_structures_fast(R, edges_df)
    
    t_pairs = time.time()
    pairs_left, pairs_right = build_pairs_sparse(R, edges_df, k_hop)
    P = len(pairs_left)
    print(f"Pairs {k_hop}-hop     : {P} (built in {time.time() - t_pairs:.2f}s)")
    
    # Initialisation intelligente de la MCMC à partir de la Baseline spectrale
    tree = np.zeros(R, dtype=np.int32)
    tau_init = (pred_baseline_aligned[:-1] != pred_baseline_aligned[1:]).astype(np.int32)
    for k in range(R - 1):
        if tau_init[k] == 1:
            fenwick_update(tree, k, 1)
            
    t_mcmc = time.time()
    correlations = mcmc_loop(steps, R, beta, tree, 
                             cross_offsets, cross_left, cross_right, cross_weight,
                             incident_offsets, incident_left, incident_right, incident_weight,
                             pairs_left, pairs_right, verbose=False)
    
    print(f"MCMC Time         : {time.time() - t_mcmc:.2f}s ({steps} steps)")
    
    # 4. MCMC Spectral Clustering par blocs
    t_eigen = time.time()
    row_C = np.concatenate([pairs_left, pairs_right])
    col_C = np.concatenate([pairs_right, pairs_left])
    data_C = np.concatenate([correlations, correlations])
    
    C_sparse = sp.coo_matrix((data_C, (row_C, col_C)), shape=(R, R)).tocsr()
    
    v_C, val_C = solve_signed_spectral_by_components(C_sparse, gpu=True, laplacian=laplacian, verbose=False)
    pred_mcmc = np.sign(v_C)
    pred_mcmc[pred_mcmc == 0] = 1
    
    # Alignement par blocs
    pred_mcmc_aligned = align_spins_by_components(pred_mcmc, labels, n_components, true_spins)
    acc_mcmc = np.mean(pred_mcmc_aligned == true_spins)
    t_eigen_total = time.time() - t_eigen
    
    t_total_mcmc = time.time() - t0
    print(f"Eigen Solver Time : {t_eigen_total:.2f}s")
    print(f"MCMC Accuracy     : {acc_mcmc:.4%}")
    print(f"Total MCMC Time   : {t_total_mcmc:.2f}s")
    
    return {
        "Instance": os.path.basename(instance_dir),
        "R": R,
        "Edges": len(edges_df),
        "Baseline Acc": acc_baseline,
        "Baseline Time": t_baseline,
        "MCMC Acc": acc_mcmc,
        "MCMC Time": t_total_mcmc
    }

if __name__ == "__main__":
    results = []
    
    # Instance 1: benchmark/difficult/instance_noisy_small
    if os.path.exists("/workspaces/Haplotypes/benchmark/difficult/instance_noisy_small"):
        res = run_instance(
            "/workspaces/Haplotypes/benchmark/difficult/instance_noisy_small",
            steps=1000000, k_hop=2, beta=1.0
        )
        results.append(res)
        
    # Instance 2: benchmark/difficult/instance_noisy_medium
    if os.path.exists("/workspaces/Haplotypes/benchmark/difficult/instance_noisy_medium"):
        res = run_instance(
            "/workspaces/Haplotypes/benchmark/difficult/instance_noisy_medium",
            steps=2000000, k_hop=2, beta=1.0
        )
        results.append(res)
        
    # Instance 3: benchmark/difficult/instance_large_noisy_15x_1
    if os.path.exists("/workspaces/Haplotypes/benchmark/difficult/instance_large_noisy_15x_1"):
        res = run_instance(
            "/workspaces/Haplotypes/benchmark/difficult/instance_large_noisy_15x_1",
            steps=5000000, k_hop=2, beta=1.0
        )
        results.append(res)
        
    # S'il n'y a pas d'instances difíciles, lancer l'instance standard par défaut
    if not results:
        res = run_instance(
            "/workspaces/Haplotypes/benchmark/small/instance_100x",
            steps=500000, k_hop=2, beta=1.0
        )
        results.append(res)
    
    df = pd.DataFrame(results)
    print("\n" + "=" * 80)
    print("FINAL BENCHMARK SUMMARY (PacBio CLR Noisy Simulation, perfectly aligned with Notebook)")
    print("=" * 80)
    print(df.to_string(index=False))
