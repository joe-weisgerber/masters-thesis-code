import numpy as np

def tight_rep_to_physical(state, s):
    if state is None:
        return None
    physical_state = []
    for i in range(len(state)):
        if i%2 == 0:
            physical_state.append(state[i]-s)
        else:
            bottom = state[i] // int(2*s+1) - s
            top = state[i] % (2*s+1) - s
            physical_state.append((bottom, top))
    return physical_state

def physical_rep_to_tight(state, s):
    if state is None:
        return None
    tight_state = []
    for i in range(len(state)):
        if i%2 == 0:
            tight_state.append(state[i]+s)
        else:
            tight_state.append((state[i][0]+s)*(2*s+1) + (state[i][1]+s))
    return tight_state

def height_rep_to_physical(state, pbc=False):
    if pbc :
        physical_state = [state[0]-state[-1], (-state[0], state[0])]
    else:
        physical_state = [state[0], (-state[0], state[0])]

    for i in range(len(state)-1):
        physical_state.append(state[i+1]- state[i])
        physical_state.append((-state[i+1], state[i+1]))
    physical_state.append(physical_state[-1][0])
    return physical_state

def apply_U(state, site, direction, s):
    # Apply raising operator U on given link
    # State is in tight representation
    if direction == 'vertical':
        new_state = state.copy()
        if new_state[2*site] >= 2*s:
            return None
        new_state[2*site] += 1
        
    elif direction == 'top':
        new_state = state.copy()
        new_state[2*site+1] += 1
        if new_state[2*site+1] % (2*s+1) == 0:
            return None
       
    elif direction == 'bottom':
        new_state = state.copy()
        new_state[2*site+1] += (2*s+1)
        if new_state[2*site+1] >= (2*s+1)*(2*s+1):
            return None
    else:
        raise ValueError("Direction must be 'vertical', 'top', or 'bottom'")
    return new_state

def apply_U_dagger(state, site, direction, s):
    # Apply raising operator U on given link
    # State is in tight representation
    if direction == 'vertical':
        new_state = state.copy()
        if new_state[2*site] <= 0:
            return None
        new_state[2*site] -= 1
        
    elif direction == 'top':
        new_state = state.copy()
        if new_state[2*site+1] % (2*s+1) == 0:
            return None
        new_state[2*site+1] -= 1
        
       
    elif direction == 'bottom':
        new_state = state.copy()
        new_state[2*site+1] -= 2*s+1
        if new_state[2*site+1] < 0:
            return None
    else:
        raise ValueError("Direction must be 'vertical', 'top', or 'bottom'")
    return new_state

def apply_U_plaquette(state, site, s, L, pbc=False):
    # Apply the plaquette term at a given site
    # State is in tight representation
    new_state = apply_U_dagger(state, site, 'bottom', s)
    if new_state is None:
        return None
    if pbc:
        new_state = apply_U_dagger(new_state, (site+1)%L, 'vertical', s)
    else:
        new_state = apply_U_dagger(new_state, site+1, 'vertical', s)
    if new_state is None:
        return None
    new_state = apply_U(new_state, site, 'top', s)
    if new_state is None:
        return None
    new_state = apply_U(new_state, site, 'vertical', s)

    return new_state

def apply_U_plaquette_dagger(state, site, s, L, pbc=False):
    # Apply the plaquette term at a given site
    # State is in tight representation
    new_state = apply_U(state, site, 'bottom', s)
    if new_state is None:
        return None
    if pbc:
        new_state = apply_U(new_state, (site+1)%L, 'vertical', s)
    else:
        new_state = apply_U(new_state, site+1, 'vertical', s)
    if new_state is None:
        return None
    new_state = apply_U_dagger(new_state, site, 'top', s)
    if new_state is None:
        return None
    new_state = apply_U_dagger(new_state, site, 'vertical', s)

    return new_state

def gauge_breaking_term(state, site, s, L, pbc=False):
    # Apply gauge breaking term on given site
    # State is in tight representation
    new_state = apply_U(state, site, 'vertical', s)
    if new_state is None:
        return None
    if pbc:
        new_state = apply_U_dagger(new_state, (site+1)%L, 'vertical', s)
    else:
        new_state = apply_U_dagger(new_state, site+1, 'vertical', s)
    return new_state

def gauge_breaking_term_dagger(state, site, s, L, pbc=False):
    # Apply gauge breaking term on given site
    # State is in tight representation
    new_state = apply_U_dagger(state, site, 'vertical', s)
    if new_state is None:
        return None
    if pbc:
        new_state = apply_U(new_state, (site+1)%L, 'vertical', s)
    else:
        new_state = apply_U(new_state, site+1, 'vertical', s)
    return new_state

def winding_breaking_term(state, site, s, L):
    # Apply winding breaking term on a given height (0 or 1)
    # State is in tight representation
    new_state = state.copy()
    if site not in [0, 1]:
        raise ValueError('Site must be 0 or 1')
    for i in range(L):
        new_state = apply_U(new_state, i, 'top' if site==1 else 'bottom', s)
        if new_state is None:
            return None
    return new_state

def winding_breaking_term_dagger(state, site, s, L):
    # Apply winding breaking term on a given height (0 or 1)
    # State is in tight representation
    new_state = state.copy()
    if site not in [0, 1]:
        raise ValueError('Site must be 0 or 1')
    for i in range(L):
        new_state = apply_U_dagger(new_state, i, 'top' if site==1 else 'bottom', s)
        if new_state is None:
            return None
    return new_state

def charge_conjugation_breaking_term(state, site, direction, s):
    # Apply charge conjugation breaking term on a given link
    # State is in tight representation
    new_state = tight_rep_to_physical(state, s)
    if direction == 'vertical':
        E = new_state[2*site]
    elif direction == 'top':
        E = new_state[2*site+1][1]
    elif direction == 'bottom':
        E = new_state[2*site+1][0]
    else:
        raise ValueError('invalid direction!')
    
    
    return (-1)**site * E

def construct_dynamically_connected_states(initial_state, s, L, H_components):
    # Construct all states dynamically connected to the initial state
    # initial state is in tight representation
    # H_components contains all the terms the summands of the hamiltonian
    states = [initial_state]
    done = False
    while done == False:
        done = True
        new_states = states.copy()

        for state in states:
            for component in H_components:
                next_state = component(state)
                
                if next_state is not None and next_state not in new_states:
                    new_states.append(next_state)
                    done = False
        states = new_states
    return states

def construct_gauge_sector(s, L, pbc=False):
    # Construct the kernel gauge sector 
    if pbc:
        initial_state = physical_rep_to_tight([0, (0, 0)] * L, s)
    else:
        initial_state = physical_rep_to_tight([0, (0, 0)] * L + [0], s)
    H_components = []
    for site in range(L):
        H_components.append(component(s, L, site, 'plaquette'))
        H_components.append(component(s, L, site, 'plaquette_dagger'))
        

    states = construct_dynamically_connected_states(initial_state, s, L, H_components)
    return states

def construct_gauge_broken_sector(s, L, gauge_breaking_sites, pbc=False):
    if pbc:
        initial_state = physical_rep_to_tight([0, (0, 0)] * L, s)
    else:
        initial_state = physical_rep_to_tight([0, (0, 0)] * L + [0], s)
    H_components = []
    for site in range(L):
        H_components.append(component(s, L, site, 'plaquette'))
        H_components.append(component(s, L, site, 'plaquette_dagger'))
    
    for site in gauge_breaking_sites:
        H_components.append(component(s, L, site, 'gauge_breaking'))
        H_components.append(component(s, L, site, 'gauge_breaking_dagger'))
        
    states = construct_dynamically_connected_states(initial_state, s, L, H_components)
    return states

def construct_full_basis(s, L):
    basis = [[i] for i in range(-s, s+1)]
    for l in range(1, 2*L+1):
        new_basis = []
        for b in basis:
            if l%2 == 0:
                for i in range(-s, s+1):
                    new_basis.append(b+[i])
            else:
                for i in range(-s, s+1):
                    for j in range(-s, s+1):
                        new_basis.append(b+[(i, j)])

        basis = new_basis
    return basis
class component:
    def __init__(self, s, L, site, component_type):
        self.s = s
        self.L = L
        self.site = site
        self.type = component_type
        if component_type == 'plaquette':
            self.func = lambda state : apply_U_plaquette(state, self.site, self.s, self.L)
        elif component_type == 'plaquette_dagger':
            self.func = lambda state : apply_U_plaquette_dagger(state, self.site, self.s, self.L)
        elif component_type == 'gauge_breaking':
            self.func = lambda state : gauge_breaking_term(state, self.site, self.s, self.L)
        elif component_type == 'gauge_breaking_dagger':
            self.func = lambda state : gauge_breaking_term_dagger(state, self.site, self.s, self.L)
        elif component_type == 'winding_breaking':
            self.func = lambda state : winding_breaking_term(state, self.site, self.s, self.L)
        elif component_type == 'winding_breaking_dagger':
            self.func = lambda state : winding_breaking_term_dagger(state, self.site, self.s, self.L)
        else:
            raise ValueError('Component does not exist!')
    def __call__(self, state):
        return self.func(state)

