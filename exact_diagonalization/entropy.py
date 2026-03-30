import numpy as np
import scipy.linalg

# Helper functions to compute the entropy
def mapping(full_basis, basis_A):
    split = len(basis_A[0]) # pos of first site corresponding to system B
    mapping_dict = {}
    for i, a in enumerate(full_basis):
        for j, b in enumerate(full_basis):
            a1 = a[:split]
            a2 = a[split:]
            b1 = b[:split]
            b2 = b[split:]
            if a2 == b2:
      
                idx_a = basis_A.index(a1)
                idx_b = basis_A.index(b1)
                mapping_dict[(i, j)] = (idx_a, idx_b)
    return mapping_dict

def reduced_density_matrix(psi, mapping_dict, dim_A):
    rho_A = np.zeros((dim_A, dim_A), dtype=complex)
    for (i, j) in mapping_dict.keys():
        rho_A[mapping_dict[(i, j)][0], mapping_dict[(i, j)][1]] += psi[i] * np.conj(psi[j])
    return rho_A

def entropy(psi, mapping_dict, dim_A):
    rho_A = reduced_density_matrix(psi, mapping_dict, dim_A)
    eigvals, _ = scipy.linalg.eigh(rho_A)
    eigvals = eigvals[eigvals > 1e-12]
    S = -np.sum(eigvals * np.log(eigvals))
    return S




    