using ITensors, ITensorMPS, HDF5
include("hamiltonian.jl")
include("system.jl")
include("time_evolution.jl")
include("observables.jl")

function measure_all(filepath_states, filepath_metadata)
    h5 = h5open(filepath_states, "r")
    
    state_names = keys(h5)
    L = parse(Int, match(r"/L=(\d+)", filepath_states).captures[1])
    close(h5)

    states = []

    for name in state_names
        state = h5open(filepath_states, "r") do f
            read(f, name, MPS)
        end
        push!(states, state)
    end
    

    times = h5open(filepath_metadata, "r") do f
            read(f, "times")
        end
    gauge_breaking_sites = !occursin("gb=0_", filepath_states)
    p = parse(Float64, match(r"_p=([\d.]+)", filepath_states).captures[1])
    g = parse(Float64, match(r"_g=([\d.]+)", filepath_states).captures[1])
    gb = parse(Float64, match(r"_gb=([\d.]+)", filepath_states).captures[1])

    sites = siteinds(states[1])
    gauge_breaking_strengths = gb

    if gauge_breaking_sites==true
        gauge_breaking_sites = [i for i in 1:L+1]
        k = gauge_breaking_strengths
        gauge_breaking_strengths = [k for _ in gauge_breaking_sites]
    else
        gauge_breaking_sites = []
    end
    if isa(gauge_breaking_strengths, Number)
        gauge_breaking_strengths = []
    end

    symmetry_breaking_sites = [i for i in 1:L+1 if i ∈ gauge_breaking_sites || i+1 ∈ gauge_breaking_sites]

    Ψ = create_scar_state(L, 0, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    H = hamiltonian(1, L, sites; gauge_breaking_sites=gauge_breaking_sites, gauge_breaking_strengths=gauge_breaking_strengths, potential=p, g=g)
    measure_energy(; state) = abs(inner(state', H, state))
    measure_fidelity(; state) = abs(inner(Ψ, state))^2
    measure_entropy(; state) = entropy(deepcopy(state); symmetry_breaking_sites=symmetry_breaking_sites)
    M = magnetization(L, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    measure_magnetization(; state) = real(inner(state', M, state))
    H0 = hamiltonian_0(1, L, sites; gauge_breaking_sites=gauge_breaking_sites)
    measure_bare_hamiltonian(; state) = real(inner(state', H0, state))
    CC = charge_conjugation(Ψ, L, sites; gauge_breaking_sites=gauge_breaking_sites)
    measure_charge_conjugation(; state) = real(inner(state', CC, state))
    V = potential(1, L, sites; gauge_breaking_sites=gauge_breaking_sites)
    measure_potential(; state) = real(inner(state', V, state))
    UU = gauge_breaking(1, L, sites; gauge_breaking_sites=gauge_breaking_sites)
    measure_gauge_breaking(; state) = real(inner(state', UU, state))
    S01 = scar_01(L, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    measure_scar_01(; state) = real(inner(state', S01, state))
    S000 = scar_000(L, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    measure_scar_000(; state) = real(inner(state', S000, state))
    MS = minus_projector(L, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    measure_minus_projector(; state) = real(inner(state', MS, state))
    measure_middle_linkdim(; state) = linkdim(state, div(L, 2))

    energies = Float64[]
    fidelities = Float64[]
    entropies = Float64[]
    magnetizations = Float64[]
    bare_energies = Float64[]
    charge_conjugations = Float64[]
    potentials = Float64[]
    gauge_breakings = Float64[]
    scar_01s = Float64[]
    scar_000s = Float64[]
    minus_projectors = Float64[]
    maxlinkdims = Int[]
    middle_linkdims = Int[]

    for i in 1:length(times)
        push!(energies, measure_energy(state=states[i]))
        push!(fidelities, measure_fidelity(state=states[i]))
        push!(entropies, measure_entropy(state=states[i]))
        push!(magnetizations, measure_magnetization(state=states[i]))
        push!(bare_energies, measure_bare_hamiltonian(state=states[i]))
        push!(charge_conjugations, measure_charge_conjugation(state=states[i]))
        push!(potentials, measure_potential(state=states[i]))
        push!(gauge_breakings, measure_gauge_breaking(state=states[i]))
        push!(scar_01s, measure_scar_01(state=states[i]))
        push!(maxlinkdims, maxlinkdim(states[i]))
        push!(scar_000s, measure_scar_000(state=states[i]))
        push!(minus_projectors, measure_minus_projector(state=states[i]))
        push!(middle_linkdims, measure_middle_linkdim(state=states[i]))
    end

    results = Dict(
        "times" => times,
        "fidelities" => fidelities,
        "magnetizations" => magnetizations,
        "entropies" => entropies,
        "bare_energies" => bare_energies,
        "potentials" => potentials,
        "charge_conjugations" => charge_conjugations,
        "gauge_breakings" => gauge_breakings,
        "scar_01s" => scar_01s,
        "energies" => energies,
        "maxlinkdims" => maxlinkdims,
        "scar_000s" => scar_000s,
        "minus_projectors" => minus_projectors,
        "middle_linkdims" => middle_linkdims
    )

    return results
end


function measure_singular_distribution(filepath_states, filepath_metadata)
    h5 = h5open(filepath_states, "r")
    
    state_names = keys(h5)
    L = parse(Int, match(r"/L=(\d+)", filepath_states).captures[1])
    close(h5)

    states = []

    for name in state_names
        state = h5open(filepath_states, "r") do f
            read(f, name, MPS)
        end
        push!(states, state)
    end
    

    times = h5open(filepath_metadata, "r") do f
            read(f, "times")
        end
    gauge_breaking_sites = !occursin("gb=0_", filepath_states)
    p = parse(Float64, match(r"_p=([\d.]+)", filepath_states).captures[1])
    g = parse(Float64, match(r"_g=([\d.]+)", filepath_states).captures[1])
    gb = parse(Float64, match(r"_gb=([\d.]+)", filepath_states).captures[1])

    sites = siteinds(states[1])
    gauge_breaking_strengths = gb

    if gauge_breaking_sites==true
        gauge_breaking_sites = [i for i in 1:L+1]
        k = gauge_breaking_strengths
        gauge_breaking_strengths = [k for _ in gauge_breaking_sites]
    else
        gauge_breaking_sites = []
    end
    if isa(gauge_breaking_strengths, Number)
        gauge_breaking_strengths = []
    end

    symmetry_breaking_sites = [i for i in 1:L+1 if i ∈ gauge_breaking_sites || i+1 ∈ gauge_breaking_sites]

    Ψ = create_scar_state(L, 0, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    singular_vals = []
    for i in 1:length(times)
        sv = singular_distribution(Ψ; symmetry_breaking_sites=symmetry_breaking_sites)
        push!(singular_vals, sv)
    end
    return singular_vals
end

function measure_state(state; p=0, g=0, gauge_breaking_sites=[], gauge_breaking_strengths=[])
    sites = siteinds(state)
    L = length(state)
    symmetry_breaking_sites = [i for i in 1:L+1 if i ∈ gauge_breaking_sites || i-1 ∈ gauge_breaking_sites]
    gauge_breaking_sites = gauge_breaking_sites
    gauge_breaking_strengths = gauge_breaking_strengths

    Ψ = create_scar_state(L, 0, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    H = hamiltonian(1, L, sites; gauge_breaking_sites=gauge_breaking_sites, gauge_breaking_strengths=gauge_breaking_strengths, potential=p, g=g)
    measure_energy(; state) = abs(inner(state', H, state))
    measure_fidelity(; state) = abs(inner(Ψ, state))^2
    measure_entropy(; state) = entropy(deepcopy(state); symmetry_breaking_sites=symmetry_breaking_sites)
    M = magnetization(L, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    measure_magnetization(; state) = real(inner(state', M, state))
    H0 = hamiltonian_0(1, L, sites; gauge_breaking_sites=gauge_breaking_sites)
    measure_bare_hamiltonian(; state) = real(inner(state', H0, state))
    CC = charge_conjugation(Ψ, L, sites; gauge_breaking_sites=gauge_breaking_sites)
    measure_charge_conjugation(; state) = real(inner(state', CC, state))
    V = potential(1, L, sites; gauge_breaking_sites=gauge_breaking_sites)
    measure_potential(; state) = real(inner(state', V, state))
    UU = gauge_breaking(1, L, sites; gauge_breaking_sites=gauge_breaking_sites)
    measure_gauge_breaking(; state) = real(inner(state', UU, state))
    S01 = scar_01(L, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    measure_scar_01(; state) = real(inner(state', S01, state))
    S000 = scar_000(L, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    measure_scar_000(; state) = real(inner(state', S000, state))
    MS = minus_projector(L, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    measure_minus_projector(; state) = real(inner(state', MS, state))
    CCS = CC_sector(sites)
    measure_CC_sector(; state) = real(inner(state', CCS, state))
    measure_middle_linkdim(; state) = linkdim(state, div(L, 2))
    CP = constraint_projector(L, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    measure_constraint_projector(; state) = real(inner(state', CP, state))

    results = Dict(
        "fidelity" => measure_fidelity(state=state),
        "magnetization" => measure_magnetization(state=state),
        "entropy" => measure_entropy(state=state),
        "bare_energies" => measure_bare_hamiltonian(state=state),
        "potentials" => measure_potential(state=state),
        "charge_conjugations" => measure_charge_conjugation(state=state),
        "gauge_breakings" => measure_gauge_breaking(state=state),
        "scar_01s" => measure_scar_01(state=state),
        "energies" => measure_energy(state=state),
        "maxlinkdims" => maxlinkdim(state),
        "scar_000s" => measure_scar_000(state=state),
        "minus_projectors" => measure_minus_projector(state=state),
        "CC_sector" => measure_CC_sector(state=state),
        "middle_linkdim" => measure_middle_linkdim(state=state),
        "constraint_projector" => measure_constraint_projector(state=state)
    )

    return results

end

function measure_rest(state; p=0, g=0, gauge_breaking_sites=[], gauge_breaking_strengths=[])
    sites = siteinds(state)
    L = length(state)
    symmetry_breaking_sites = [i for i in 1:L+1 if i ∈ gauge_breaking_sites || i-1 ∈ gauge_breaking_sites]
    gauge_breaking_sites = gauge_breaking_sites
    gauge_breaking_strengths = gauge_breaking_strengths

    Ψ = create_scar_state(L, 0, sites; symmetry_breaking_sites=symmetry_breaking_sites)

    CCS = CC_sector(sites)
    measure_CC_sector(; state) = real(inner(state', CCS, state))

    CP = constraint_projector(L, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    measure_constraint_projector(; state) = real(inner(state', CP, state))

    results = Dict(
        "constraint_projector" => measure_constraint_projector(state=state)
    )

    return results

end