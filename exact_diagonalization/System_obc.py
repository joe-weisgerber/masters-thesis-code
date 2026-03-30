from state_space_obc import *
import numpy as np
import scipy.linalg
from scipy.sparse import csr_matrix
import numpy as np
from entropy import *

def construct_hamiltonian(basis, components, factors):
    N = len(basis)
    H = np.zeros((N, N), dtype = complex)
    for i in range(N):
        state = basis[i]
        for j in range(len(components)):
            new_state = components[j](state)
            if new_state is not None:
                if new_state not in basis:
                    print('Error: new state not in basis!')
                    continue
                k = basis.index(new_state)
                H[k, i] += factors[j]
    return H

class System:
    def __init__(self, s, L, perturbation_terms='', factors='', charge_conjugation_breaking_strength=0, charge_conjugation_breaking_pos='everywhere', potential=False, hermitian=True, height_potential=False, constraint=False):
    
        # Perturbation terms is of the form [['types'], ['sites']]

        self.s = s
        self.L = L
        self.charge_conjugation_breaking_strength = charge_conjugation_breaking_strength
        self.charge_conjugation_breaking_pos = charge_conjugation_breaking_pos
        self.potential = potential
        self.height_potential = height_potential
        self.components = []
        for site in range(L):
            self.components.append(component(s, L, site, 'plaquette'))
            self.components.append(component(s, L, site, 'plaquette_dagger'))
        self.factors = [1] * (2*L) # 1 for each U and Ud
        if perturbation_terms != '':
            perturbation_terms = list(perturbation_terms)
            perturbation_terms[0] = list(perturbation_terms[0])
            perturbation_terms[1] = list(perturbation_terms[1])
            factors = list(factors)
            if hermitian:
                perturbation_terms[0] += [term + '_dagger' for term in perturbation_terms[0]]
                perturbation_terms[1] += perturbation_terms[1]
                factors = factors + [np.conj(factor) for factor in factors]
            for i in range(len(factors)):
                self.components.append(component(s, L, perturbation_terms[1][i], perturbation_terms[0][i]))
                self.factors.append(factors[i])
        self.construct_basis()
        self.N = len(self.basis)
        self.construct_hamiltonian()
        self.gauge_breaking_sites = []
        if perturbation_terms != '':
            self.gauge_breaking_sites = perturbation_terms[1]
        self.symmetry_breaking_sites = [i for i in range(L+1) if i in self.gauge_breaking_sites or i-1 in self.gauge_breaking_sites]
        
    def construct_basis(self):
        self.basis = construct_dynamically_connected_states(physical_rep_to_tight([0, (0, 0)] * self.L + [0], self.s),
                                                            self.s, self.L, self.components)
    def construct_hamiltonian(self):
        self.H = construct_hamiltonian(self.basis, self.components, self.factors)
        if self.potential != False:
            for i in range(self.N):
                potential_energy = 0
                state = tight_rep_to_physical(self.basis[i], self.s)
                for j in range(self.L*2+1):
                    if j%2 == 0:
                        potential_energy += state[j]**2
                    else:
                        potential_energy += state[j][0]**2
                        potential_energy += state[j][1]**2
                potential_energy *= self.potential
                self.H[i, i] += potential_energy
        if self.height_potential != False:
            for i in range(self.N):
                potential_energy = 0
                state = tight_rep_to_physical(self.basis[i], self.s)
                for j in range(self.L*2+1):
                    if j%2 == 1:
                        potential_energy += state[j][0]
                potential_energy *= self.height_potential
                self.H[i, i] += potential_energy
                
        if self.charge_conjugation_breaking_strength!=0:
            for i in range(self.N):
                breaking = 0
                if self.charge_conjugation_breaking_pos == 'everywhere':
                    for site in range(self.L):
                        breaking += charge_conjugation_breaking_term(self.basis[i], site, 'vertical', self.s)
                        breaking += charge_conjugation_breaking_term(self.basis[i], site, 'top', self.s)
                        breaking += charge_conjugation_breaking_term(self.basis[i], site, 'bottom', self.s)
                    breaking += charge_conjugation_breaking_term(self.basis[i], self.L, 'vertical', self.s)
                else:
                    for pos in self.charge_conjugation_breaking_pos:
                        breaking += charge_conjugation_breaking_term(self.basis[i], pos[0], pos[1], self.s)
                    
                self.H[i, i] += breaking*self.charge_conjugation_breaking_strength

    def construct_full_system(self):
        full_basis_physical = construct_full_basis(self.s, self.L)
        full_basis = []
        for b in full_basis_physical:
            full_basis.append(physical_rep_to_tight(b, self.s))
        full_H = construct_hamiltonian(full_basis, self.components, self.factors)
        if self.charge_conjugation_breaking!=0:
            for i in range(self.N):
                breaking = 0
                for site in range(self.L):
                    breaking += charge_conjugation_breaking_term(full_basis[i], site, 'vertical', self.s)
                    breaking += charge_conjugation_breaking_term(full_basis[i], site, 'top', self.s)
                    breaking += charge_conjugation_breaking_term(full_basis[i], site, 'bottom', self.s)
                breaking += charge_conjugation_breaking_term(full_basis[i], self.L+1, 'vertical', self.s)
                full_H[i, i] += breaking*self.charge_conjugation_breaking
        return full_basis, full_H
    
    def construct_potential(self):
        self.V = np.zeros((self.N, self.N), dtype=complex)
        for i in range(self.N):
            potential_energy = 0
            state = tight_rep_to_physical(self.basis[i], self.s)
            for j in range(self.L*2+1):
                if j%2 == 0:
                    potential_energy += state[j]**2
                else:
                    potential_energy += state[j][0]**2
                    potential_energy += state[j][1]**2
            potential_energy *= self.potential
            self.V[i, i] += potential_energy


    def get_basis_states(self, physical=True):
        if physical:
            return [tight_rep_to_physical(state, self.s) for state in self.basis]
        else:
            return self.basis
    
    def get_hamiltonian(self):
        return self.H
    def get_hamiltonian_element(self, i, j):
        return self.H[i, j]
    
    def diagonalize(self):
        self.eigenvalues, self.eigenvectors = np.linalg.eigh(self.H)
        return self.eigenvalues, self.eigenvectors
    
    def get_eigenvalues(self):
        return self.eigenvalues
    def get_eigenvectors(self):
        return self.eigenvectors
    
    def time_evolve_all_fidelities(self, T, steps, return_us=False, verbose=False):
        # Computes the fidelities of computational basis states
        fidelities = np.zeros((steps+1, self.N), dtype=complex)
        fidelities[0, ...] = np.ones(self.N)
        dt = T / steps
        
        H = self.H
        U = scipy.linalg.expm(-1j * H * dt)
    
        
        Um = np.eye(self.N, dtype = complex)
        if return_us:
            Us = np.zeros((steps+1, self.N, self.N), dtype=complex)
            Us[0, ...] = Um.copy()
        else:
            Us = None

        for i in range(steps):
            Um = Um @ U
            diagonal = np.diag(Um)
            fidelities[i+1, ...] = np.abs(diagonal) ** 2
            if verbose and np.isclose(dt * (i + 1), round(dt * (i + 1))):
                print(f"t = {dt * (i + 1)}")
            if return_us:
                Us[i+1, ...] = Um.copy()
        self.fidelities = fidelities
        self.Us = Us
        self.time = np.linspace(0, T, steps)
        return self.fidelities, self.Us
    
    def find_scars(self, T, steps, number_of_scars=5, return_us=False):
        # Extract revival scars
        self.time_evolve_all_fidelities(T, steps, return_us=return_us)
        integrated_fidelities = np.trapz(self.fidelities, dx=T/steps, axis=0)
        sorted_indeces = np.argsort(integrated_fidelities)[::-1]
        if return_us:
            return sorted_indeces[:number_of_scars], self.fidelities, self.Us
        else:
            return sorted_indeces[:number_of_scars], self.fidelities
    
    def time_evolve_state_fidelity(self, initial_state, T, steps):
        # Computes the fidelity of any state
        if not hasattr(self, 'fidelities') or not hasattr(self, 'Us'):
            print('Time evolving all fidelities...')
            self.time_evolve_all_fidelities(T, steps, return_us=True)
        elif len(self.fidelities) != steps:
            self.time_evolve_all_fidelities(T, steps, return_us=True)
        
        
        state = physical_rep_to_tight(initial_state, self.s)
        if state not in self.basis:
            print('Error : initial state not in basis!')
            return None
        index = self.basis.index(state)
        return self.fidelities[:, index]
        
    
    def construct_rho(self):
        rho = np.zeros((self.N, self.N))
        for i in range(self.N):
            state =  tight_rep_to_physical(np.array(self.basis[i]), self.s)
            rho[i, i] = np.sum([np.sum(np.abs(x)) if isinstance(x, tuple) else abs(x) for x in state])
        rho /= self.L
        self.rho = rho
        return rho
    
    def construct_scar01(self):
        scar01 = np.zeros((self.N, self.N))
        for i in range(self.N):
            state =  tight_rep_to_physical(np.array(self.basis[i]), self.s)
            reduced_state = [state[2*j+1][0] for j in range(self.L)]
            scar01[i, i] += np.sum([1 if reduced_state[2*j]==0 and reduced_state[2*j+1]==1 else 0 for j in range(self.L//2)])
            scar01[i, i] += np.sum([1 if reduced_state[2*j]==1 and reduced_state[2*j+1]==0 else 0 for j in range(self.L//2)])
            scar01[i, i] += np.sum([1 if reduced_state[2*j]==0 and reduced_state[2*j+1]==-1 else 0 for j in range(self.L//2)])
            scar01[i, i] += np.sum([1 if reduced_state[2*j]==-1 and reduced_state[2*j+1]==0 else 0 for j in range(self.L//2)])
        scar01 *= 2/self.L
        self.scar01 = scar01
        return scar01
    
    def time_evolve_rho(self, T, steps):
        rho = self.construct_rho()
        rho_t = np.zeros((steps, self.N, self.N), dtype=complex)
        if not hasattr(self, 'Us'):
            print('Constructing time evolution operators...')
            self.time_evolve_all_fidelities(T, steps, return_us=True)
        
        for i in range(steps):
            rho_t[i, ...] = np.diag(self.Us[i, ...] @ rho @ self.Us[i, ...].conj().T)
        self.rho_t = rho_t
        return rho_t
    
    def time_evolve_scar01(self, T, steps):
        scar01 = self.construct_scar01()
        scar01_t = np.zeros((steps, self.N, self.N), dtype=complex)
        if not hasattr(self, 'Us'):
            print('Constructing time evolution operators...')
            self.time_evolve_all_fidelities(T, steps, return_us=True)
        
        for i in range(steps):
            scar01_t[i, ...] = np.diag(self.Us[i, ...] @ scar01 @ self.Us[i, ...].conj().T)
        self.scar01_t = scar01_t
        return scar01_t
    
    def magnetization(self, i):
        if not hasattr(self, 'rho'):
            self.construct_rho()
        return self.rho[i, i]
    
    def time_evolve_rho_state(self, initial_state, T, steps):
        if not hasattr(self, 'rho_t'):
            self.time_evolve_rho(T, steps)
        state = physical_rep_to_tight(initial_state, self.s)
        if state not in self.basis:
            print('Error : initial state not in basis!')
            return None
        index = self.basis.index(state)
        return self.rho_t[:, index, index]
    
    def time_evolve_scar01_state(self, initial_state, T, steps):
        if not hasattr(self, 'scar01_t'):
            self.time_evolve_scar01(T, steps)
        state = physical_rep_to_tight(initial_state, self.s)
        if state not in self.basis:
            print('Error : initial state not in basis!')
            return None
        index = self.basis.index(state)
        return self.scar01_t[:, index, index]
    
    def get_vector(self, states, factors):
        factors = np.array(factors)
        factors = factors / np.linalg.norm(factors)
        vec = np.zeros(self.N, dtype=complex)
        for i in range(len(states)):
            state = physical_rep_to_tight(states[i], self.s)
            if state not in self.basis:
                raise ValueError('State not in basis!')
            index = self.basis.index(state)
            vec[index] = factors[i]
        return vec
    
    def time_evolve_vector(self, T, steps, states='', factors='', scar='', verbose=False):
        dt = T / steps
        if not hasattr(self, 'Us'):
            print('Constructing time evolution operators...')
            self.time_evolve_all_fidelities(T, steps, return_us=True, verbose=verbose)
        if scar != '':
            states, factors = scar.get_states_factors()

        vec = self.get_vector(states, factors)
        vec_t = np.zeros((steps+1, self.N), dtype=complex)
        print('Evolving vector...')
        for i in range(steps+1):
            if verbose and np.isclose(dt * (i + 1), round(dt * (i + 1))):
                print(f"t = {self.time[i]}")
            vec_t[i, ...] = self.Us[i, ...] @ vec
        self.vec_t = vec_t
        return vec_t
    
    def time_evolve_vector_fidelity(self, T, steps, states='', factors='', scar=''):
        # Time evolves any state, most relevant for zero-mode scars
        vec_t = self.time_evolve_vector(T, steps, states, factors, scar)
        if scar != '':
            states, factors = scar.get_states_factors()
        vec = self.get_vector(states, factors)
        fidelity = np.abs(vec_t @ vec)**2
        return fidelity
    
    def time_evolve_vector_rho(self, T, steps, states='', factors='', scar='', verbose=False):
        #rho is what we call M in the thesis, i.e. the projector onto the +-1 links
        vec_t = self.time_evolve_vector(T, steps, states, factors, scar, verbose=verbose)
        rho = self.construct_rho()
        rho_t = np.zeros(steps+1)
        for t in range(steps+1):
            rho_t[t] = vec_t[t, ...].conj().T @ rho @ vec_t[t, ...]
        return rho_t
    
    def time_evolve_vector_scar01(self, T, steps, states='', factors='', scar='', verbose=False):
        # scar01 is the observable O1 in the thesis that projects locally on the 
        # space spanned by 01, 10, -10, 0-1
        vec_t = self.time_evolve_vector(T, steps, states, factors, scar, verbose=verbose)
        scar01 = self.construct_scar01()
        scar01_t = np.zeros(steps+1)
        for t in range(steps+1):
            scar01_t[t] = vec_t[t, ...].conj().T @ scar01 @ vec_t[t, ...]
        return scar01_t

    def get_basis_reduced(self):
        basis = [tight_rep_to_physical(state, self.s) for state in self.basis]
        reduced_basis = []
        for b in basis:
            rb = []
            for i in range(len(b)):
                if i%2 == 0:
                    rb.append(b[i])
                else:
                    rb.append(b[i][0])
            reduced_basis.append(rb)
        return reduced_basis
    def get_shortened_basis(self): 
        basis = self.get_basis_reduced()
        optimal_basis = []
        for b in basis:
            ob = []
            for i in range(len(b)):
                if i%2 == 1:
                    ob.append(b[i])
                elif i%2 == 0 and i//2 in self.symmetry_breaking_sites:
                    ob.append(b[i])
            if self.L+1 in self.symmetry_breaking_sites:
                ob.append(b[-1])
            optimal_basis.append(ob)
        return optimal_basis
    
    def get_partial_trace_mapping(self):
        basis_a = []
        for b in self.basis:
            br = b[:self.L]
            if br not in basis_a:
                basis_a.append(br)
        mapping_dict = mapping(self.basis, basis_a)
        return mapping_dict, len(basis_a)
    
    def get_partial_trace_mapping_shortened(self):
        basis = self.get_shortened_basis()
        split = int(self.L/2 + len([site for site in self.symmetry_breaking_sites if site < self.L/2]))
        print("Split : ", split)
        basis_a = []
        for b in basis:
            br = b[:split]
            if br not in basis_a:
                basis_a.append(br)
        mapping_dict = mapping(basis, basis_a)

        print(len(basis_a[0]))
        return mapping_dict, len(basis_a)

    def time_evolve_vector_entropy(self, T, steps, states='', factors='', scar='', shortened=True):
        vec_t = self.time_evolve_vector(T, steps, states, factors, scar)
        if shortened:
            mapping_dict, dim_A = self.get_partial_trace_mapping_shortened()
        else:
            mapping_dict, dim_A = self.get_partial_trace_mapping()
        entropies = np.zeros(steps+1)

        for t in range(steps+1):
            entropies[t] = entropy(vec_t[t, ...], mapping_dict, dim_A)
        return entropies
    
    def spectrum_entropy(self, shortened=True):
        if shortened:
            mapping_dict, dim_A = self.get_partial_trace_mapping_shortened()
        else:
            mapping_dict, dim_A = self.get_partial_trace_mapping()
        eigenvalues, eigenvectors = self.diagonalize()
        
        entropies = np.zeros(len(eigenvalues))
        
        for i in range(len(eigenvalues)):
            entropies[i] = entropy(eigenvectors.T[i], mapping_dict, dim_A)

        return entropies, eigenvalues
    
    def scar_entropy(self, scar, shortened=True):
        states, factors = scar.get_states_factors()
        vec = self.get_vector(states, factors)
        if shortened:
            mapping_dict, dim_A = self.get_partial_trace_mapping_shortened()
        else:
            mapping_dict, dim_A = self.get_partial_trace_mapping()
        S = entropy(vec, mapping_dict, dim_A)
        return S
    
    def vector_entropy(self, vec, shortened=True):
        if not shortened:
            print("Not supported yet!")
        else:
            mapping_dict, dim_A = self.get_partial_trace_mapping_shortened()
        S = entropy(vec, mapping_dict, dim_A)
        return S
    
    def reduced_density_matrix(self, vec):
        mapping_dict, dim_A = self.get_partial_trace_mapping()
        rho_A = reduced_density_matrix(vec, mapping_dict, dim_A)
        return rho_A
    
    def microcanonical_magnetization(self, dE, E0=0, scar='', states='', factors=''):
        if not hasattr(self, 'eigenvalues') or not hasattr(self, 'eigenvectors'):
            self.diagonalize()
        if scar != '':
            states, factors = scar.get_states_factors()
            vec_scar = self.get_vector(states, factors)
            E0 = vec_scar.conj().T @ self.H @ vec_scar
        if states != '' and factors != '':
            vec_states = self.get_vector(states, factors)
            E0 = vec_states.conj().T @ self.H @ vec_states
        indices = [i for i in range(len(self.eigenvalues)) if E0 - dE <= self.eigenvalues[i] <= E0 + dE]
        if len(indices) == 0:
            print('No eigenstates in the given energy window!')
            return None
        else:
            print(f'Number of eigenstates in the energy window: {len(indices)}')
        M = self.construct_rho()
        avg_magnetization = 0
        for i in indices:
            vec = self.eigenvectors.T[i]
            avg_magnetization += vec.conj().T @ M @ vec
        avg_magnetization /= len(indices)
        return np.real(avg_magnetization)
    
    def microcanonical_scar01(self, dE, E0=0, scar='', states='', factors=''):
        if not hasattr(self, 'eigenvalues') or not hasattr(self, 'eigenvectors'):
            self.diagonalize()
        if scar != '':
            states, factors = scar.get_states_factors()
            vec_scar = self.get_vector(states, factors)
            E0 = vec_scar.conj().T @ self.H @ vec_scar
        if states != '' and factors != '':
            vec_states = self.get_vector(states, factors)
            E0 = vec_states.conj().T @ self.H @ vec_states
        indices = [i for i in range(len(self.eigenvalues)) if E0 - dE <= self.eigenvalues[i] <= E0 + dE]
        if len(indices) == 0:
            print('No eigenstates in the given energy window!')
            return None
        else:
            print(f'Number of eigenstates in the energy window: {len(indices)}')
        scar01 = self.construct_scar01()
        avg_scar01 = 0
        for i in indices:
            vec = self.eigenvectors.T[i]
            avg_scar01 += vec.conj().T @ scar01 @ vec
        avg_scar01 /= len(indices)
        return np.real(avg_scar01)

    def prethermal_scar01(self, dE, E0=0, scar='', states='', factors=''):
        self.construct_potential()
        if scar != '':
            states, factors = scar.get_states_factors()
            vec_scar = self.get_vector(states, factors)
            E0 = vec_scar.conj().T @ self.V @ vec_scar
        if states != '' and factors != '':
            vec_states = self.get_vector(states, factors)
            E0 = vec_states.conj().T @ self.V @ vec_states
        print("energy : ", E0)
        indices = [i for i in range(self.N) if E0 - dE <= self.V[i, i] <= E0 + dE]
        if len(indices) == 0:
            print('No eigenstates in the given energy window!')
            return None
        else:
            print(f'Number of eigenstates in the energy window: {len(indices)}')
        scar01 = self.construct_scar01()
        avg_scar01 = 0
        for i in indices:
            avg_scar01 += scar01[i, i]
        avg_scar01 /= len(indices)
        return np.real(avg_scar01)
    
    def get_cc_sector(self, sector=1):
        U = np.zeros((self.N, (self.N+sector)//2), dtype=complex)
        k = 0
        basis = self.get_basis_states()
        for i in range(self.N):
            el1 = basis[i]
            el2 = []

            for b in el1:    
                if isinstance(b, tuple):
                    el2.append((-b[0], -b[1]))
                else:
                    el2.append(-b)
    
            j = basis.index(el2)
            if j == i and sector==1:
                U[i, k] = 1
                k += 1
            elif j > i:
                U[i, k] = 1/np.sqrt(2)
                U[j, k] = sector/np.sqrt(2)
                k += 1
        return U.conj().T @ self.H @ U

    import numpy as np

    def get_reflection_sector(self, sector=1):
        columns = []
        basis = self.get_basis_states()
        
        for i in range(self.N):
            el1 = basis[i]
            el2 = el1[::-1]
            
            for l in range(len(el2)):
                if isinstance(el2[l], tuple):
                    b = el2[l]
                    el2[l] = (-b[0], -b[1])
                    
            j = basis.index(el2)

            if j == i and sector == 1:
                col = np.zeros(self.N, dtype=complex)
                col[i] = 1
                columns.append(col)
            elif j > i:
                col = np.zeros(self.N, dtype=complex)
                col[i] = 1 / np.sqrt(2)
                col[j] = sector / np.sqrt(2)
                columns.append(col)

        U = np.column_stack(columns)
        
        return U.conj().T @ self.H @ U
    
    def scar_overlaps(self):

        scar = Scar(self.s, self.L, 2, 'even')
        vec = scar.get_vector_scar()
        self.diagonalize()
        overlaps = np.abs(self.eigenvectors.conj().T @ vec)**2
        return overlaps

    def scar_projector_overlaps(self, tol=1e-10, i=2):
        # Returns the overlap of the zero-mode scar
        #  with eigenspaces of the perturbed Hamiltonian 
        scar = Scar(self.s, self.L, i, 'even')
        vec = scar.get_vector_scar()
        
        self.diagonalize()
        evals = self.eigenvalues
        evecs = self.eigenvectors
        
        individual_overlaps = np.abs(evecs.conj().T @ vec)**2
        
        unique_evals = []
        weights = []
        
        idx = np.argsort(evals)
        sorted_evals = evals[idx]
        sorted_overlaps = individual_overlaps[idx]

        current_val = sorted_evals[0]
        current_weight = sorted_overlaps[0]
        
        for i in range(1, len(sorted_evals)):
            if np.abs(sorted_evals[i] - current_val) < tol:
                current_weight += sorted_overlaps[i]
            else:
                unique_evals.append(current_val)
                weights.append(current_weight)
                current_val = sorted_evals[i]
                current_weight = sorted_overlaps[i]
        
        unique_evals.append(current_val)
        weights.append(current_weight)
        return np.array(unique_evals), np.array(weights)

class Scar:
    # Class to manage zero-mode scars
    # Using i > S gives the superposition of scars
    # in the +1 charge conjugation sector
    def __init__(self, s, L, i, tiling):
        self.s = s
        self.L = L
        self.i = i
        self.tiling = tiling #'even' or 'odd', but OBC requires even
        self.construct()
    def construct(self):
        if self.i > self.s:
            scar0 = Scar(self.s, self.L, 0, tiling=self.tiling)
            scar1 = Scar(self.s, self.L, 1, tiling=self.tiling)
            states0, factors0 = scar0.get_states_factors()
            states1, factors1 = scar1.get_states_factors()
            self.states = states0 + states1
            self.factors = np.array(list(factors0) + list(factors1))
            self.factors = self.factors / np.linalg.norm(self.factors)
            return self.states, self.factors
        else:
            if self.L%2 == 1:
                raise ValueError('L must be even!')
            base_tiles = []
            base_factors = []
            for k in range(self.s+1):
                base_tiles.append([self.i-self.s+k, self.i-k])
                base_factors.append((-1)**k)
            
            states = base_tiles.copy()
            factors = base_factors.copy()

            for l in range(self.L//2 - 1):
                new_states = []
                new_factors = []
                for x in range(len(states)):
                    for y in range(len(base_tiles)):
                        new_states.append(states[x] + base_tiles[y])
                        new_factors.append(factors[x]*base_factors[y])
                states = new_states
                factors = new_factors

            self.states = []
            for s in states:
                if self.tiling == 'even':
                    self.states.append(height_rep_to_physical(s))
                elif self.tiling == 'odd':
                    state = height_rep_to_physical(s)
                    state += state[:2].copy()
                    del state[1]
                    del state[0]
                    self.states.append(state)
            self.factors = np.array(factors)
            self.factors = self.factors / np.linalg.norm(self.factors)
            return self.states, self.factors
        
    def get_states_factors(self):
        return self.states, self.factors
    
    def get_tensor_scar(self):
        states = []
        for s in self.states:
            new_s = []
            for i in range(len(s)):
                if i%2 == 0:
                    new_s.append(s[i])
                else:
                    new_s.append(s[i][0])
            states.append(new_s)
        return states, self.factors
    
    def get_vector_scar(self):
        system = System(self.s, self.L)
        basis = system.get_basis_states()
        vec = np.zeros(len(basis), dtype=complex)
        for i, state in enumerate(self.states):
            if state not in basis:
                raise ValueError('State not in basis!')
            idx = basis.index(state)
            vec[idx] = self.factors[i]
        norm = np.linalg.norm(vec)
        if norm == 0:
            raise ValueError('Resulting scar vector has zero norm!')
        vec = vec / norm
        return vec
    
    def get_shortened_tensor_scar(self, gauge_breaking_sites):
        states = []
        for s in self.states:
            counter = 0
            new_s = []

            for i in range(self.L):

                if i in gauge_breaking_sites or (i-1)%self.L in gauge_breaking_sites:
                    new_s.append(s[counter])    
                counter += 1

                new_s.append(s[counter][0])
                counter += 1

            states.append(new_s)
        return states, self.factors


        
