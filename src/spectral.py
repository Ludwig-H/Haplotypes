import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse.linalg import lobpcg

# Détection et import conditionnel de CuPy pour accélération GPU
GPU_AVAILABLE = False
try:
    import cupy as cp
    import cupyx.scipy.sparse as cps
    import cupyx.scipy.sparse.linalg as cpsla
    GPU_AVAILABLE = True
except ImportError:
    pass

def compute_signed_laplacian(W):
    r"""
    Calcule le Laplacien signé d'une matrice d'adjacence signée W.
    L = D - W_sym, où D_ii = \sum_j |W_ij| et W_sym = W + W^T.
    """
    W_abs = abs(W)
    deg_cols = np.array(W_abs.sum(axis=0)).flatten()
    deg_rows = np.array(W_abs.sum(axis=1)).flatten()
    degrees = deg_cols + deg_rows
    
    D = sp.diags(degrees)
    W_sym = W + W.T
    L = D - W_sym
    return L, degrees

def signed_spectral_clustering(W, type="unnormalized", use_gpu=True):
    """
    Effectue le clustering spectral signé sur la matrice W.
    Exclut automatiquement les nœuds isolés (degré 0) pour stabiliser la diagonalisation,
    puis ré-injecte les spins par défaut pour les nœuds isolés.
    """
    R = W.shape[0]
    
    # Étape 1 : Identifier et exclure les nœuds de degré 0 (isolés)
    W_abs = abs(W)
    deg_cols = np.array(W_abs.sum(axis=0)).flatten()
    deg_rows = np.array(W_abs.sum(axis=1)).flatten()
    degrees = deg_cols + deg_rows
    active_nodes = np.where(degrees > 0)[0]
    R_active = len(active_nodes)
    
    if R_active < 2:
        return np.ones(R, dtype=np.float64)
        
    # Extraire la sous-matrice pour les nœuds connectés uniquement
    W_active = W[active_nodes, :][:, active_nodes]
    
    # Étape 2 : Calculer le Laplacien sur la sous-matrice active
    L_active, degrees_active = compute_signed_laplacian(W_active)
    
    if type == "normalized":
        degrees_clip = np.clip(degrees_active, 1e-6, None)
        d_inv_sqrt = 1.0 / np.sqrt(degrees_clip)
        D_inv_sqrt = sp.diags(d_inv_sqrt)
        L_active = D_inv_sqrt.dot(L_active).dot(D_inv_sqrt)
        
    run_on_gpu = GPU_AVAILABLE and use_gpu
    
    v_active = None
    if run_on_gpu:
        print(f"🚀 [GPU] Diagonalisation (CuPy) sur {R_active} nœuds actifs...")
        try:
            L_gpu = cps.csr_matrix(L_active)
            if type == "normalized":
                I_gpu = cps.eye(R_active, format='csr')
                M_gpu = 2.0 * I_gpu - L_gpu
                eigenvalues, eigenvectors = cpsla.eigsh(M_gpu, k=1, which='LA')
                v_active = cp.asnumpy(eigenvectors[:, 0])
            else:
                sigma_gpu = float(2.0 * cp.max(cp.array(degrees_active)))
                I_gpu = cps.eye(R_active, format='csr')
                M_gpu = sigma_gpu * I_gpu - L_gpu
                eigenvalues, eigenvectors = cpsla.eigsh(M_gpu, k=1, which='LA')
                v_active = cp.asnumpy(eigenvectors[:, 0])
            print("✅ Diagonalisation GPU terminée.")
        except Exception as e:
            print(f"⚠️ Échec GPU : {e}. Fallback CPU...")
            run_on_gpu = False
            
    if not run_on_gpu:
        print(f"💻 [CPU] Diagonalisation (LOBPCG) sur {R_active} nœuds actifs...")
        np.random.seed(42)
        X = np.random.normal(size=(R_active, 1))
        try:
            vals, vecs = lobpcg(L_active, X, largest=False, tol=1e-5, maxiter=200)
            v_active = vecs[:, 0]
            print("✅ Diagonalisation CPU (LOBPCG) terminée.")
        except Exception as e:
            print(f"⚠️ [solve_signed_spectral CPU] LOBPCG a échoué: {e}. Essai avec shift-invert ARPACK...")
            try:
                I = sp.eye(R_active, format='csr')
                if type == "normalized":
                    M = 2.0 * I - L_active
                    eigenvalues, eigenvectors = spla.eigsh(M, k=1, which='LA', tol=1e-5)
                    v_active = eigenvectors[:, 0]
                else:
                    sigma = float(2.0 * np.max(degrees_active))
                    M = sigma * I - L_active
                    eigenvalues, eigenvectors = spla.eigsh(M, k=1, which='LA', tol=1e-5)
                    v_active = eigenvectors[:, 0]
                print("✅ Diagonalisation CPU (ARPACK) terminée.")
            except Exception as e2:
                print(f"⚠️ Échec ARPACK : {e2}.")
                if R_active < 5000:
                    print("⚠️ Repli sur solveur dense (graphe de petite taille)...")
                    L_dense = L_active.toarray()
                    vals, vecs = np.linalg.eigh(L_dense)
                    v_active = vecs[:, 0]
                else:
                    # Graphe trop grand pour le mode dense, on retourne un vecteur uniforme
                    print("❌ Graphe trop grand pour le mode dense. Utilisation d'un vecteur uniforme.")
                    v_active = np.ones(R_active, dtype=np.float64)
                    
    # Étape 3 : Reconstruire le vecteur de spins global
    v_global = np.ones(R, dtype=np.float64)
    if v_active is not None:
        v_global[active_nodes] = v_active
        
    return v_global
