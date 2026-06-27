import pandas as pd
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import numba
import time
import os

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
    row = np.concatenate([edges_df['source'].values, edges_df['target'].values])
    col = np.concatenate([edges_df['target'].values, edges_df['source'].values])
    data = np.ones(len(row), dtype=np.bool_)
    
    A = sp.coo_matrix((data, (row, col)), shape=(R, R)).tocsr()
    I = sp.eye(R, dtype=np.bool_, format='csr')
    
    Visited = (A + I)
    Current = Visited
    for _ in range(k - 1):
        Current = Current @ Visited
        
    Current_tri = sp.triu(Current, k=1)
    coo = Current_tri.tocoo()
    return coo.row.astype(np.int32), coo.col.astype(np.int32)

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
              pairs_left, pairs_right, pair_sums, sample_freq):
    
    num_samples = 0
    for t in range(1, steps + 1):
        r = np.random.randint(0, R)
        
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
            
        if chosen_move == 1:
            if r - 1 >= 0:
                fenwick_update(tree, r-1, 1)
            if r < R - 1:
                fenwick_update(tree, r, 1)
        elif chosen_move == 2:
            if r - 1 >= 0:
                fenwick_update(tree, r-1, 1)
        elif chosen_move == 3:
            if r < R - 1:
                fenwick_update(tree, r, 1)
                
        # Periodically sample state
        if t % sample_freq == 0:
            spins = np.zeros(R, dtype=np.float64)
            for i in range(R):
                xor_val = fenwick_query(tree, i-1)
                spins[i] = 1.0 - 2.0 * float(xor_val)
                
            for p in range(len(pairs_left)):
                u = pairs_left[p]
                v = pairs_right[p]
                pair_sums[p] += spins[u] * spins[v]
            num_samples += 1
            
    return num_samples

def evaluate_predictions(pred, true_spins):
    acc = np.mean(pred == true_spins)
    acc = max(acc, 1.0 - acc)
    return acc

def run_instance(instance_dir, steps, k_hop, beta, sample_freq):
    print("=" * 60)
    print(f"RUNNING BENCHMARK ON: {instance_dir}")
    print("=" * 60)
    
    R, edges_df, true_spins = load_instance(instance_dir)
    print(f"Nodes R: {R}, Edges: {len(edges_df)}")
    
    # 1. Baseline: Signed Spectral Clustering
    t0 = time.time()
    row = np.concatenate([edges_df['source'].values, edges_df['target'].values])
    col = np.concatenate([edges_df['target'].values, edges_df['source'].values])
    data = np.concatenate([edges_df['weight'].values, edges_df['weight'].values])
    
    W_sparse = sp.coo_matrix((data, (row, col)), shape=(R, R)).tocsr()
    abs_W = sp.coo_matrix((np.abs(data), (row, col)), shape=(R, R)).tocsr()
    degrees_W = np.array(abs_W.sum(axis=1)).flatten()
    
    D_W = sp.diags(degrees_W)
    L_W = D_W - W_sparse
    
    vals_W, vecs_W = eigsh(L_W, k=1, which='SM')
    v_W = vecs_W[:, 0]
    pred_baseline = np.sign(v_W)
    acc_baseline = evaluate_predictions(pred_baseline, true_spins)
    t_baseline = time.time() - t0
    print(f"Baseline Accuracy : {acc_baseline:.4%}")
    print(f"Baseline Time     : {t_baseline:.2f}s")
    
    # 2. Glauber MCMC
    t0 = time.time()
    cross_offsets, cross_left, cross_right, cross_weight, \
    incident_offsets, incident_left, incident_right, incident_weight = build_structures_fast(R, edges_df)
    
    t_pairs = time.time()
    pairs_left, pairs_right = build_pairs_sparse(R, edges_df, k_hop)
    P = len(pairs_left)
    print(f"Pairs {k_hop}-hop     : {P} (built in {time.time() - t_pairs:.2f}s)")
    
    tree = np.zeros(R, dtype=np.int32)
    pair_sums = np.zeros(P, dtype=np.float64)
    
    t_mcmc = time.time()
    num_samples = mcmc_loop(steps, R, beta, tree, 
                            cross_offsets, cross_left, cross_right, cross_weight,
                            incident_offsets, incident_left, incident_right, incident_weight,
                            pairs_left, pairs_right, pair_sums, sample_freq)
    
    correlations = pair_sums / float(num_samples)
    print(f"MCMC Time         : {time.time() - t_mcmc:.2f}s ({steps} steps, {num_samples} samples)")
    
    # 3. MCMC Spectral Clustering
    t_eigen = time.time()
    row_C = np.concatenate([pairs_left, pairs_right])
    col_C = np.concatenate([pairs_right, pairs_left])
    data_C = np.concatenate([correlations, correlations])
    
    C_sparse = sp.coo_matrix((data_C, (row_C, col_C)), shape=(R, R)).tocsr()
    abs_C = sp.coo_matrix((np.abs(data_C), (row_C, col_C)), shape=(R, R)).tocsr()
    degrees_C = np.array(abs_C.sum(axis=1)).flatten()
    
    D_C = sp.diags(degrees_C)
    L_C = D_C - C_sparse
    
    vals_C, vecs_C = eigsh(L_C, k=1, which='SM')
    v_C = vecs_C[:, 0]
    pred_mcmc = np.sign(v_C)
    acc_mcmc = evaluate_predictions(pred_mcmc, true_spins)
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
            steps=100000, k_hop=2, beta=1.5, sample_freq=200
        )
        results.append(res)
        
    # Instance 2: benchmark/difficult/instance_noisy_medium
    if os.path.exists("/workspaces/Haplotypes/benchmark/difficult/instance_noisy_medium"):
        res = run_instance(
            "/workspaces/Haplotypes/benchmark/difficult/instance_noisy_medium",
            steps=150000, k_hop=2, beta=1.2, sample_freq=300
        )
        results.append(res)
        
    # Instance 3: benchmark/difficult/instance_large_noisy_15x_1
    if os.path.exists("/workspaces/Haplotypes/benchmark/difficult/instance_large_noisy_15x_1"):
        res = run_instance(
            "/workspaces/Haplotypes/benchmark/difficult/instance_large_noisy_15x_1",
            steps=200000, k_hop=2, beta=1.0, sample_freq=400
        )
        results.append(res)
        
    # S'il n'y a pas d'instances difíciles, lancer les instances standard par défaut
    if not results:
        res = run_instance(
            "/workspaces/Haplotypes/benchmark/small/instance_100x",
            steps=50000, k_hop=2, beta=2.0, sample_freq=100
        )
        results.append(res)
    
    df = pd.DataFrame(results)
    print("\n" + "=" * 80)
    print("FINAL BENCHMARK SUMMARY (PacBio CLR Noisy Simulation)")
    print("=" * 80)
    print(df.to_string(index=False))
