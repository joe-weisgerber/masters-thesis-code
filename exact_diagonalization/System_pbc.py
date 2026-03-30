from state_space import *
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
    def __init__(self, s, L, perturbation_terms='', factors='', hermitian=True, charge_conjugation_breaking_strength=0, charge_conjugation_breaking_pos='everywhere', potential=False):
        # Perturbation terms is of the form [['types'], ['sites']]
        self.s = s
        self.L = L
        self.charge_conjugation_breaking_strength = charge_conjugation_breaking_strength
        self.charge_conjugation_breaking_pos = charge_conjugation_breaking_pos
        self.potential = potential
        self.components = []
        for site in range(L):
            self.components.append(component(s, L, site, 'plaquette'))
            self.components.append(component(s, L, site, 'plaquette_dagger'))
        self.factors = [1] * (2*L)
        if perturbation_terms != '':
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
        self.basis = construct_dynamically_connected_states(physical_rep_to_tight([0, (0, 0)] * self.L, self.s),
                                                            self.s, self.L, self.components)

    def construct_hamiltonian(self):
        self.H = construct_hamiltonian(self.basis, self.components, self.factors)
        if self.potential!=False:
            for i in range(self.N):
                potential_energy = 0
                state = tight_rep_to_physical(self.basis[i], self.s)
                for j in range(self.L*2):
                    if j%2 == 0:
                        potential_energy += state[j]**2
                    else:
                        potential_energy += state[j][0]**2
                        potential_energy += state[j][1]**2
                potential_energy *= self.potential
                self.H[i, i] += potential_energy
                
        if self.charge_conjugation_breaking_strength!=0:
            for i in range(self.N):
                breaking = 0
                if self.charge_conjugation_breaking_pos == 'everywhere':
                    for site in range(self.L):
                        breaking += charge_conjugation_breaking_term(self.basis[i], site, 'vertical', self.s)
                        breaking += charge_conjugation_breaking_term(self.basis[i], site, 'top', self.s)
                        breaking += charge_conjugation_breaking_term(self.basis[i], site, 'bottom', self.s)
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
        if self.charge_conjugation_breaking_strength!=0:
            for i in range(self.N):
                breaking = 0
                for site in range(self.L):
                    breaking += charge_conjugation_breaking_term(full_basis[i], site, 'vertical', self.s)
                    breaking += charge_conjugation_breaking_term(full_basis[i], site, 'top', self.s)
                    breaking += charge_conjugation_breaking_term(full_basis[i], site, 'bottom', self.s)
                full_H[i, i] += breaking*self.charge_conjugation_breaking
        return full_basis, full_H

    def get_basis_states(self, physical = True):
        if physical:
            return [tight_rep_to_physical(state, self.s) for state in self.basis]
        else:
            return self.basis
    
    def get_hamiltonian(self):
        return self.H
    
    def diagonalize(self):
        self.eigenvalues, self.eigenvectors = np.linalg.eigh(self.H)
        return self.eigenvalues, self.eigenvectors
    
    def get_eigenvalues(self):
        return self.eigenvalues
    def get_eigenvectors(self):
        return self.eigenvectors
    
    def time_evolve_all_fidelities(self, T, steps, return_us=False, verbose=False):
        
        fidelities = np.zeros((steps, self.N), dtype=complex)
        fidelities[0, ...] = np.ones(self.N)
        dt = T / steps
        
        H = self.H
        U = scipy.linalg.expm(-1j * H * dt)
    
        
        Um = np.eye(self.N, dtype = complex)
        if return_us:
            Us = np.zeros((steps, self.N, self.N), dtype=complex)
            Us[0, ...] = Um.copy()
        else:
            Us = None

        for i in range(steps-1):
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
        self.time_evolve_all_fidelities(T, steps, return_us=return_us)
        integrated_fidelities = np.trapz(self.fidelities, dx=T/steps, axis=0)
        sorted_indeces = np.argsort(integrated_fidelities)[::-1]
        if return_us:
            return sorted_indeces[:number_of_scars], self.fidelities, self.Us
        else:
            return sorted_indeces[:number_of_scars], self.fidelities
    
    def time_evolve_state_fidelity(self, initial_state, T, steps):

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
    
    def time_evolve_rho_state(self, initial_state, T, steps):
        if not hasattr(self, 'rho_t'):
            self.time_evolve_rho(T, steps)
        state = physical_rep_to_tight(initial_state, self.s)
        if state not in self.basis:
            print('Error : initial state not in basis!')
            return None
        index = self.basis.index(state)
        return self.rho_t[:, index, index]
    
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
            self.time_evolve_all_fidelities(T, steps, return_us=True)
        if scar != '':
            states, factors = scar.get_states_factors()
            
        vec = self.get_vector(states, factors)
        vec_t = np.zeros((steps, self.N), dtype=complex)
        print('Evolving vector...')
        for i in range(steps):
            if verbose and np.isclose(dt * (i + 1), round(dt * (i + 1))):
                print(f"t = {self.time[i]}")
            vec_t[i, ...] = self.Us[i, ...] @ vec
        self.vec_t = vec_t
        return vec_t
    
    def time_evolve_vector_fidelity(self, T, steps, states='', factors='', scar=''):
        vec_t = self.time_evolve_vector(T, steps, states, factors, scar)
        if scar != '':
            states, factors = scar.get_states_factors()
        vec = self.get_vector(states, factors)
        sparse_vec = csr_matrix(vec)
        sparse_vect = csr_matrix(vec_t)
        fidelity = np.abs(vec_t @ vec)**2
        return fidelity
    
    def time_evolve_vector_rho(self, T, steps, states='', factors='', scar=''):
        vec_t = self.time_evolve_vector(T, steps, states, factors, scar)
        rho = self.construct_rho()
        rho_t = np.zeros(steps)
        for t in range(steps):
            rho_t[t] = vec_t[t, ...].conj().T @ rho @ vec_t[t, ...]
        return rho_t
    
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
    
    def get_partial_trace_mapping(self):
        basis_a = []
        for b in self.basis:
            br = b[:self.L]
            if br not in basis_a:
                basis_a.append(br)
        mapping_dict = mapping(self.basis, basis_a)
        return mapping_dict, len(basis_a)
    
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
            optimal_basis.append(ob)
        return optimal_basis
    
    def get_partial_trace_mapping_shortened(self):
        basis = self.get_shortened_basis()
        split = int(self.L/2 + len([site for site in self.symmetry_breaking_sites if site <= self.L/2]))
        basis_a = []
        for b in basis:
            br = b[:split]
            if br not in basis_a:
                basis_a.append(br)
        mapping_dict = mapping(basis, basis_a)
        return mapping_dict, len(basis_a)
    
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


    

class Scar:
    def __init__(self, s, L, i, tiling):
        self.s = s
        self.L = L
        self.i = i
        self.tiling = tiling #'even' or 'odd' 
        self.construct()
    def construct(self):
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

