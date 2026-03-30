using ITensors, ITensorMPS
using Observers: observer
using HDF5
include("observables.jl")

"""Implements time evolution with TDVP, so far only for the fidelity"""


function tdvp_step(Ψ, H_eff, tau, cutoff, maxdim)
    Ψ = tdvp(H_eff, -im * tau, Ψ;
        nsteps=1,
        cutoff = cutoff,
        maxdim = maxdim,
        normalize = false,
        reverse_step = true,
        outputlevel = 0,
    )

    return Ψ
end


function tdvp_projector(Ψ, H, P, T, steps; observable="fidelity", cutoff=1e-10, maxdim=50)
    tau = T/steps
    PHP = apply(P, H, P)
    if observable == "fidelity"
        fidelities = Float64[]
        Ψ0 = deepcopy(Ψ)
        for t in 0.0:tau:T
            F = abs(inner(Ψ0, Ψ))^2
            push!(fidelities, F)
            t ≈ T && break
            if isapprox(t, round(t); atol=1e-8)
                println("t = $(round(t, digits=3)), F = $(round(F, digits=6)), max bond dim = $(maxlinkdim(Ψ))")
            end
                
            Ψ = tdvp(H_eff, -im * tau, Ψ;
                nsteps=1,
                cutoff = cutoff,
                maxdim = maxdim,
                normalize = false,
                reverse_step = true,
                outputlevel = 0,
            )
        end
        return fidelities
    else
        println("Observable not supported")
    end
end

function tdvp_new(Ψ, H, T, steps; observable="fidelity", cutoff=1e-10, maxdim=50)
    tau = T/steps
    if observable == "fidelity"
        fidelities = Float64[]
        Ψ0 = deepcopy(Ψ)
        for t in 0.0:tau:T
            F = abs(inner(Ψ0, Ψ))^2
            push!(fidelities, F)
            t ≈ T && break
            if isapprox(t, round(t); atol=1e-8)
                println("t = $(round(t, digits=3)), F = $(round(F, digits=6)), max bond dim = $(maxlinkdim(Ψ))")
            end
                
            Ψ = tdvp(H_eff, -im * tau, Ψ;
                nsteps=1,
                cutoff = cutoff,
                maxdim = maxdim,
                normalize = false,
                reverse_step = true,
                outputlevel = 0,
            )
        end
        return fidelities
    else
        println("Observable not supported")
    end
end



function tdvp_observer(Ψ, H, T, steps; cutoff=1e-12, maxdim=50, observable="fidelity", symmetry_breaking_sites=[], L=0, sites="")

    tau = T / steps
    
    Ψ₀ = deepcopy(Ψ)

    

    if observable == "fidelity"
        measure_fidelity(; state) = abs(inner(Ψ₀, state))^2
        obs = observer("fidelities" => measure_fidelity)
    elseif observable == "entropy"
        measure_entropy(; state) = entropy(deepcopy(state); symmetry_breaking_sites=symmetry_breaking_sites)
        obs = observer("entropies" => measure_entropy)
    elseif observable == "magnetization"
        M = magnetization(L, sites; symmetry_breaking_sites=symmetry_breaking_sites)
        measure_magnetization(; state) = real(inner(state', M, state))
        obs = observer("magnetizations" => measure_magnetization)
    elseif observable == "state"
        measure_state(; state) = deepcopy(state)
        obs = observer("states" => measure_state)
    else
        error("Observable not supported")
    end

    tdvp(H, -im * T, Ψ; 
        time_step = -im * tau,
        cutoff = cutoff, 
        maxdim = maxdim,
        normalize = true,
        step_observer! = obs, 
        outputlevel = 1
    )
    
    times = collect(0:steps) * tau
    
    if observable == "fidelity"
        fidelities = [1.0; obs.fidelities]
        return times, fidelities
    elseif observable == "entropy"
        entropies = [entropy(Ψ₀; symmetry_breaking_sites=symmetry_breaking_sites); obs.entropies]
        return times, entropies
    elseif observable == "magnetization"
        magnetizations = [real(inner(Ψ₀, M, Ψ₀')); obs.magnetizations]
        return times, magnetizations
    elseif observable == "state"
        states = [deepcopy(Ψ₀); obs.states]
        return times, states
    end
    
end

function thermalize(Ψ, H, tau, T, gauge_breaking_sites, L, sites; cutoff=1e-12, maxdim=50, stop=true)
    symmetry_breaking_sites = [i for i in 1:L+1 if i ∈ gauge_breaking_sites || i-1 ∈ gauge_breaking_sites]
    Ψ0 = deepcopy(Ψ)
    measure_fidelity(; state) = abs(inner(Ψ0, state))^2
    measure_state(; state) = deepcopy(state)
    measure_entropy(; state) = entropy(deepcopy(state); symmetry_breaking_sites=symmetry_breaking_sites)
    measure_maxlinkdim(; state) = maxlinkdim(state)
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


    fidelities = [1.0]
    states = [deepcopy(Ψ0)]
    magnetizations = [measure_magnetization(state=Ψ0)]
    entropies = [measure_entropy(state=Ψ0)]
    maxlinkdims = [maxlinkdim(Ψ0)]
    bare_hamiltonians = [measure_bare_hamiltonian(state=Ψ0)]
    charge_conjugations = [measure_charge_conjugation(state=Ψ0)]
    potentials = [measure_potential(state=Ψ0)]
    gauge_breakings = [measure_gauge_breaking(state=Ψ0)]
    Scar_01s = [measure_scar_01(state=Ψ0)]

    times = [0.0]
    t = 0.0
    tdvp_times = Float64[]

    while t < T && ((maxlinkdims[end] < maxdim && stop==true) || stop==false)
        t += tau
        if maxlinkdims[end] < 4
            Ψ = expand(Ψ, H; alg="global_krylov")
            normalize!(Ψ)
        end

        elapsed = @elapsed begin
            Ψ = tdvp(H, -im * tau, Ψ; 
                nsteps = 1,
                cutoff = cutoff, 
                maxdim = maxdim,
                normalize = true,
                outputlevel = 0
            )
        end
        push!(tdvp_times, elapsed)
        push!(fidelities, measure_fidelity(state=Ψ))
        push!(states, measure_state(state=Ψ))
        push!(magnetizations, measure_magnetization(state=Ψ))
        push!(entropies, measure_entropy(state=Ψ))
        push!(maxlinkdims, measure_maxlinkdim(state=Ψ))
        push!(bare_hamiltonians, measure_bare_hamiltonian(state=Ψ))
        push!(charge_conjugations, measure_charge_conjugation(state=Ψ))
        push!(potentials, measure_potential(state=Ψ))
        push!(gauge_breakings, measure_gauge_breaking(state=Ψ))
        push!(Scar_01s, measure_scar_01(state=Ψ))
        push!(times, round(t, digits=1))
        println("t = $(round(t, digits=3)), F = $(round(fidelities[end], digits=6)), maxbonddim = $(maxlinkdim(Ψ)), sweep duration = $(round(elapsed, digits=2)) seconds")
    end

    results = Dict(
        "times" => times,
        "fidelities" => fidelities,
        "magnetizations" => magnetizations,
        "entropies" => entropies,
        "maxlinkdims" => maxlinkdims,
        "bare_hamiltonians" => bare_hamiltonians,
        "charge_conjugations" => charge_conjugations,
        "potentials" => potentials,
        "gauge_breakings" => gauge_breakings,
        "Scar_01s" => Scar_01s,
        "tdvp_times" => tdvp_times
    )
    return results, states;
end


function thermalize_states(Ψ, H, tau, T, gauge_breaking_sites, L, sites; cutoff=1e-12, maxdim=50, stop=true)
    symmetry_breaking_sites = [i for i in 1:L+1 if i ∈ gauge_breaking_sites || i-1 ∈ gauge_breaking_sites]
    Ψ0 = deepcopy(Ψ)
    measure_maxlinkdim(; state) = maxlinkdim(state)
    
    maxlinkdims = [maxlinkdim(Ψ0)]
    states = [deepcopy(Ψ0)]

    times = [0.0]
    t = 0.0
    tdvp_times = Float64[]

    while t < T && ((maxlinkdims[end] < maxdim && stop==true) || stop==false)
        t += tau
        if maxlinkdims[end] < 4
            Ψ = expand(Ψ, H; alg="global_krylov")
            normalize!(Ψ)
        end

        elapsed = @elapsed begin
            Ψ = tdvp(H, -im * tau, Ψ; 
                nsteps = 1,
                cutoff = cutoff, 
                maxdim = maxdim,
                normalize = true,
                outputlevel = 0
            )
        end
        push!(tdvp_times, elapsed)
        push!(maxlinkdims, measure_maxlinkdim(state=Ψ))
        push!(times, round(t, digits=1))
        push!(states, deepcopy(Ψ))
        println("t = $(round(t, digits=3)), maxbonddim = $(maxlinkdim(Ψ)), sweep duration = $(round(elapsed, digits=2)) seconds")
    end

    results = Dict(
        "times" => times,
        "maxlinkdims" => maxlinkdims,
        "tdvp_times" => tdvp_times,
    )
    return results, states;
end

function thermalize_and_measure(Ψ, H, tau, T, gauge_breaking_sites, L, sites, p, g, gb, cb, folder; cutoff=1e-12, maxdim=50, stop=true)
    symmetry_breaking_sites = [i for i in 1:L+1 if i ∈ gauge_breaking_sites || i-1 ∈ gauge_breaking_sites]
    Ψ0 = deepcopy(Ψ)

    meas = measure_state(Ψ0; p=p, g=g, gauge_breaking_sites=gauge_breaking_sites, gauge_breaking_strengths=[])
    results = Dict()
    
    for k in keys(meas)
        results[k] = [meas[k]]
    end
    results["times"] = [0.0]

    t=0.0
    maxlinkdims = [maxlinkdim(Ψ0)]
    while t < T && ((maxlinkdims[end] < maxdim && stop==true) || stop==false)
        t += tau
        if maxlinkdims[end] < 4
            Ψ = expand(Ψ, H; alg="global_krylov")
            normalize!(Ψ)
        end

        elapsed = @elapsed begin
            Ψ = tdvp(H, -im * tau, Ψ; 
                nsteps = 1,
                cutoff = cutoff, 
                maxdim = maxdim,
                normalize = true,
                outputlevel = 0
            )
        end
        push!(maxlinkdims, maxlinkdim(Ψ))
        meas = measure_state(Ψ)
        meas["times"] = round(t, digits=1)
        for k in keys(results)
            push!(results[k], meas[k])
        end

        h5open("$(folder)/L=$(L)_p=$(p)_g=$(g)_gb=$(gb)_cb=$(cb)_cutoff=$(cutoff)_maxdim=$(maxdim)_stop=$(stop)_states_t=$(round(t, digits=2)).h5", "w") do file
            write(file, "state", Ψ)
        end
        
        println("t = $(round(t, digits=3)), maxbonddim = $(maxlinkdim(Ψ)), sweep duration = $(round(elapsed, digits=2)) seconds")
    end

    return results
end

function measure_all_sites(state, L, sites)
    results = Dict()
    
    for i in 1:L÷2
        obs = local_scar(L, sites, 2*i-1)
        results["$(i)"] = real(inner(state', obs, state))
    end
    return results
end

function  thermalize_and_measure_velocity(Ψ, H, tau, T, p, g, gb, gauge_breaking_sites, L, sites, folder; cutoff=1e-12, maxdim=50, stop=true)
    symmetry_breaking_sites = [i for i in 1:L+1 if i ∈ gauge_breaking_sites || i-1 ∈ gauge_breaking_sites]
    Ψ0 = deepcopy(Ψ)

    results = Dict()
    meas = measure_all_sites(Ψ0, L, sites)
    for k in keys(meas)
        results[k] = [meas[k]]
    end
    
    results["times"] = [0.0]
    maxlinkdims = [maxlinkdim(Ψ0)]
    t = 0.0
    while t < T && ((maxlinkdims[end] < maxdim && stop==true) || stop==false)
        t += tau
        if maxlinkdims[end] < 4
            Ψ = expand(Ψ, H; alg="global_krylov")
            normalize!(Ψ)
        end

        elapsed = @elapsed begin
            Ψ = tdvp(H, -im * tau, Ψ; 
                nsteps = 1,
                cutoff = cutoff, 
                maxdim = maxdim,
                normalize = true,
                outputlevel = 0
            )
        end

        meas = measure_all_sites(Ψ, L, sites)
        for k in keys(meas)
            push!(results[k], meas[k])
        end
        push!(results["times"], round(t, digits=1))
        push!(maxlinkdims, maxlinkdim(Ψ))

        h5open("$(folder)/L=$(L)_p=$(p)_g=$(g)_gb=$(gb)_cutoff=$(cutoff)_maxdim=$(maxdim)_stop=$(stop)_states_t=$(round(t, digits=2)).h5", "w") do file
            write(file, "state", Ψ)
        end
        
        println("t = $(round(t, digits=3)), maxbonddim = $(maxlinkdim(Ψ)), sweep duration = $(round(elapsed, digits=2)) seconds")
    end

    return results
end
