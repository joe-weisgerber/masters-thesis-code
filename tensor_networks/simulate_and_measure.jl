using ITensors, ITensorMPS, HDF5, Printf
include("hamiltonian.jl")
include("system.jl")
include("time_evolution.jl")
include("Measure_all.jl")

function generate_data(L, T, tau, g, p, gauge_breaking_sites, gauge_breaking_strengths, cb, maxdim, cutoff, stop, folder)
    gb = gauge_breaking_strengths
    if isa(gauge_breaking_sites, Bool)
            
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
    end

    symmetry_breaking_sites = [i for i in 1:L+1 if i in gauge_breaking_sites || i-1 in gauge_breaking_sites]
    sites = get_sites(1, L; symmetry_breaking_sites=symmetry_breaking_sites)
    Ψ = create_scar_state(L, 0, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    H = hamiltonian(1, L, sites; gauge_breaking_sites=gauge_breaking_sites, gauge_breaking_strengths=gauge_breaking_strengths, potential=p, g=g, constraint_breaking=cb)
    Ψ = expand(Ψ, H; alg="global_krylov")
    normalize!(Ψ)

    results = thermalize_and_measure(Ψ, H, tau, T, symmetry_breaking_sites, L, sites, p, g, gb, cb, folder; maxdim=maxdim, cutoff=cutoff, stop=stop)

    h5open("$(folder)/results$(L)_p=$(p)_g=$(g)_gb=$(gb)_cb=$(cb)_T=$(T)_tau=$(tau)_cutoff=$(cutoff)_maxdim=$(maxdim)_stop=$(stop)_measurements.h5","w") do f
        for (name, series) in results
            write(f, name, series)
        end
    end
end
    

# Prompt user for parameters
println("Enter L (integer): ")
L = parse(Int, readline())

println("Enter T (number): ")
T = parse(Float64, readline())

println("Enter tau (number): ")
tau = parse(Float64, readline())

println("Enter cutoffs (number, e.g. 1e-20): ")
cutoff_input = readline()
cutoffs = parse.(Float64, split(cutoff_input, ","))

println("Enter maxdims (integer): ")
maxdim_input = readline()
maxdims = parse.(Int64, split(maxdim_input, ","))

println("Enter stop (true/false): ")
stop = lowercase(strip(readline())) == "true"

println("Enter type of perturbation (p, g, gb or cb)")
perturbation_type = strip(readline())

println("Enter perturbation strengths (number, e.g. 0.1, or 0.0 for no perturbation): ")
strengths_input = readline()
strengths = parse.(Float64, split(strengths_input, ","))

println("Enter directory to save results (e.g. results/): ")
folder = readline()

if perturbation_type == "p"
    for p in strengths
        for (cutoff, maxdim) in zip(cutoffs, maxdims)
            println("Running simulation for potential p = ", p, ", cutoff = ", cutoff, ", maxdim = ", maxdim)
            generate_data(L, T, tau, 0, p, false, 0, 0, maxdim, cutoff, stop, folder)
            println("------------------------------------------------------")
        end
    end
end

if perturbation_type == "g"
    for g in strengths
        for (cutoff, maxdim) in zip(cutoffs, maxdims)
            println("Running simulation for charge conjugation breaking g = ", g, ", cutoff = ", cutoff, ", maxdim = ", maxdim)
            generate_data(L, T, tau, g, 0, false, 0, 0,maxdim, cutoff, stop, folder)
            println("------------------------------------------------------")
        end
    end
end


if perturbation_type == "gb"
    for b in strengths
        for (cutoff, maxdim) in zip(cutoffs, maxdims)
            println("Running simulation for gauge breaking b = ", b, ", cutoff = ", cutoff, ", maxdim = ", maxdim)
            generate_data(L, T, tau, 0, 0, true, b, 0, maxdim, cutoff, stop, folder)
            println("------------------------------------------------------")
        end
    end
end

if perturbation_type == "cb"
    for b in strengths
        for (cutoff, maxdim) in zip(cutoffs, maxdims)
            println("Running simulation for constraint breaking cb = ", b, ", cutoff = ", cutoff, ", maxdim = ", maxdim)
            generate_data(L, T, tau, 0, 0, false, 0, b, maxdim, cutoff, stop, folder)
            println("------------------------------------------------------")
        end
    end
end