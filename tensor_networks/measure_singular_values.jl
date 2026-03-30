using HDF5, ITensors, ITensorMPS

singular_values = []

println("Enter L :")
L = parse(Int, readline())
println("Enter p : ")
p = parse(Float64, readline())
println("Enter g : ")
g = parse(Float64, readline())
println("Enter gb : ")
gb = parse(Float64, readline())
println("Enter cb : ")
cb = parse(Float64, readline())
println("Enter t1 : ")
t1 = parse(Float64, readline())
println("Enter t2 : ")
t2 = parse(Float64, readline())
println("Enter folder : ")
folder = readline()

if p == 0.0
    p = Int(p)
end
if g == 0.0
    g = Int(g)
end
if gb == 0.0
    gb = Int(gb)
end
if cb == 0.0
    cb = Int(cb)
end

for t in t1:0.1:t2
    println("Processing time t = ", t)
    filename = "$(folder)/L=$(L)_p=$(p)_g=$(g)_gb=$(gb)_cb=$(cb)_cutoff=1.0e-12_maxdim=2000_stop=true_states_t=$(t).h5"
    
    f = h5open(filename, "r")
    Ψ = read(f, "state", MPS)
    close(f)
    
    b = Int(L/2)
    orthogonalize!(Ψ, b)
    U, S, V = svd(Ψ[b], (linkind(Ψ, b-1), siteind(Ψ, b)))
    sv = [S[i, i] for i in 1:dim(S, 1)]
    
    push!(singular_values, sv)
end

println("Singular values collected for all times")

# Save singular_values (vector of vectors) to an HDF5 file
out_file = "$(folder)/singular_values_L=$(L)_p=$(p)_g=$(g)_gb=$(gb)_cb=$(cb)_cutoff=1.0e-12_maxdim=2000_stop=true.h5"

times = collect(t1:0.1:t2)
@assert length(times) == length(singular_values) "times and singular_values length mismatch"

h5open(out_file, "w") do f
    # Store metadata
    f["times"] = times
    f["num_times"] = length(singular_values)

    # Store each singular-value vector as its own dataset
    g = create_group(f, "singular_values")
    for (i, sv) in enumerate(singular_values)
        g["t$(i)"] = sv
    end
end

println("Saved singular values to: ", out_file)