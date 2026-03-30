using ITensors, ITensorMPS, ITensorMPOConstruction

ITensors.op(::OpName"M", ::SiteType"S=1") = 
    [1 0 0
     0 0 0
     0 0 1]
ITensors.op(::OpName"00", ::SiteType"S=1") = 
    [0 0 0
     0 1 0
     0 0 0]
ITensors.op(::OpName"11", ::SiteType"S=1") = 
    [1 0 0
     0 0 0
     0 0 0]
ITensors.op(::OpName"-1-1", ::SiteType"S=1") = 
    [0 0 0
     0 0 0
     0 0 1]
ITensors.op(::OpName"E", ::SiteType"S=1") = 
    [1 0 0
     0 0 0
     0 0 -1]
ITensors.op(::OpName"01", ::SiteType"S=1") = 
    [0 1 0
     0 0 0
     0 0 0]
ITensors.op(::OpName"10", ::SiteType"S=1") = 
    [0 0 0
     1 0 0
     0 0 0]
ITensors.op(::OpName"0-1", ::SiteType"S=1") = 
    [0 0 0
     0 0 0
     0 1 0]
ITensors.op(::OpName"-10", ::SiteType"S=1") = 
    [0 0 0
     0 0 1
     0 0 0]
ITensors.op(::OpName"C", ::SiteType"S=1") = 
    [0 0 1
     0 1 0
     1 0 0]


function entropy(Ψ; symmetry_breaking_sites=[])
    L = length(Ψ) - length(symmetry_breaking_sites) 
    b = Int(L/2) + length(symmetry_breaking_sites[symmetry_breaking_sites .<= L/2])
    orthogonalize!(Ψ, b)
    U,S,V = svd(Ψ[b], (linkind(Ψ, b-1), siteind(Ψ,b)))
    SvN = 0.0
    for n=1:dim(S, 1)
        p = S[n,n]^2
        SvN -= p * log(p)
    end
    return SvN
end

function magnetization(L, sites; symmetry_breaking_sites=[])
    os = OpSum()
    
    counter = 1
    if 1 ∈ symmetry_breaking_sites
        os += 1/L, "M", counter
        counter += 1
        os += 2/L, "M", counter
        counter += 1
    else
        os += 3/L, "M", counter
        counter += 1
    end

    for i in 2:L-1
        if i ∉ symmetry_breaking_sites
            os += 2/L, "M", counter
            os += 1/L, "00", counter-1, "M", counter
            os += 1/L, "M", counter-1, "00", counter
            counter += 1
            
        else          
            os += 1/L, "M", counter
            counter += 1
            os += 2/L, "M", counter
            counter += 1
            
        end
    end

    if L ∉ symmetry_breaking_sites && L+1 ∉ symmetry_breaking_sites
        os += 3/L, "M", counter
        os += 1/L, "00", counter-1, "M", counter
        os += 1/L, "M", counter-1, "00", counter

        counter += 1
    elseif L ∈ symmetry_breaking_sites && L+1 ∉ symmetry_breaking_sites
        
        os += 1/L, "M", counter
        counter += 1
        os += 3/L, "M", counter
        counter += 1
    elseif L ∉ symmetry_breaking_sites && L+1 ∈ symmetry_breaking_sites
        println("Case not supported yet")
    else
        os += 1/L, "M", counter
        counter += 1
        os += 2/L, "M", counter
        counter += 1
        os += 1/L, "M", counter
    end

    
    return MPO_new(os, sites)
end

function compute_entropies(states; symmetry_breaking_sites=[])
    entropies = []
    for state in states
        S = entropy(deepcopy(state); symmetry_breaking_sites=symmetry_breaking_sites)
        push!(entropies, S)
    end
    return entropies
end

function compute_fidelities(states, Ψ0)
    fidelities = []
    for state in states
        F = abs(inner(Ψ0, state))^2
        push!(fidelities, F)
    end
    return fidelities
end

function compute_magnetizations(states, sites, L; symmetry_breaking_sites=[])
    magnetizations = []
    M = magnetization(L, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    for state in states
        m = real(inner(state', M, state))
        push!(magnetizations, m)
    end
    return magnetizations
end

function scar_01(L, sites; symmetry_breaking_sites=[])
    os = OpSum()
    counter = 1
    i = 1
    while i <= L
        if i ∈ symmetry_breaking_sites 
            counter += 1
        end
        if i+1 ∈ symmetry_breaking_sites
            os += 2/L, "00", counter, "11", counter + 2
            os += 2/L, "11", counter, "00", counter + 2
            os += 2/L, "00", counter, "-1-1", counter + 2
            os += 2/L, "-1-1", counter, "00", counter + 2
            counter += 3
        else
            os += 2/L, "00", counter, "11", counter + 1
            os += 2/L, "11", counter, "00", counter + 1
            os += 2/L, "00", counter, "-1-1", counter + 1
            os += 2/L, "-1-1", counter, "00", counter + 1
            counter += 2
        end
        i += 2
    end
    return MPO_new(os, sites)
end

function minus_projector(L, sites; symmetry_breaking_sites=[])
    os = OpSum()
    counter = 1
    i = 1
    while i <= L
        if i ∈ symmetry_breaking_sites 
            counter += 1
        end
        if i+1 ∈ symmetry_breaking_sites
            os += 1/L, "00", counter, "11", counter + 2
            os += 1/L, "11", counter, "00", counter + 2
            os += 1/L, "00", counter, "-1-1", counter + 2
            os += 1/L, "-1-1", counter, "00", counter + 2
            os += -1/L, "01", counter, "10", counter + 2
            os += -1/L, "10", counter, "01", counter + 2
            os += -1/L, "0-1", counter, "-10", counter + 2
            os += -1/L, "-10", counter, "0-1", counter + 2
            counter += 3
        else
            os += 1/L, "00", counter, "11", counter + 1
            os += 1/L, "11", counter, "00", counter + 1
            os += 1/L, "00", counter, "-1-1", counter + 1
            os += 1/L, "-1-1", counter, "00", counter + 1
            os += -1/L, "01", counter, "10", counter + 1
            os += -1/L, "10", counter, "01", counter + 1
            os += -1/L, "0-1", counter, "-10", counter + 1
            os += -1/L, "-10", counter, "0-1", counter + 1
            counter += 2
        end
        i += 2
    end
    return MPO_new(os, sites)
end

function scar_000(L, sites; symmetry_breaking_sites=[])
    os = OpSum()
    counter = 1
    for i in 1:L-2
        if i ∈ symmetry_breaking_sites
            counter += 1
        end
        cp1 = counter + 1
        cp2 = counter + 2
        if i+1 ∈ symmetry_breaking_sites
            cp1 += 1
            cp2 += 1
        end
        if i+2 ∈ symmetry_breaking_sites
            cp2 += 1
        end
        os += 1, "00", counter, "00", cp1, "00", cp2
        os += 1, "11", counter, "11", cp1, "11", cp2
        os += 1, "-1-1", counter, "-1-1", cp1, "-1-1", cp2
        counter += 1
    end
    return MPO_new(os, sites)
end

function singular_distribution(Ψ; symmetry_breaking_sites=[])
    L = length(Ψ) - length(symmetry_breaking_sites) 
    b = Int(L/2) + length(symmetry_breaking_sites[symmetry_breaking_sites .<= L/2])
    orthogonalize!(Ψ, b)
    U,S,V = svd(Ψ[b], (linkind(Ψ, b-1), siteind(Ψ,b)))
    singular_vals = []
    for n=1:dim(S, 1)
        push!(singular_vals, S[n, n])
    end
    return singular_vals
end

function local_scar(L, sites, p)
    os = OpSum()
    os += 1, "00", p, "11", p + 1
    os += 1, "11", p, "00", p + 1
    os += 1, "00", p, "-1-1", p + 1
    os += 1, "-1-1", p, "00", p + 1

    return MPO_new(os, sites)
end

function CC_sector(sites)
    return MPO(sites, "C")
end

function get_A(type)
    if type == "normal"
        s = Index(3)
        a = Index(3)
        b = Index(3)
        A = ITensor(s, a, b)
        A[s=>1, a=>1, b=>1] = 1.0
        A[s=>1, a=>1, b=>2] = 1.0
        A[s=>2, a=>2, b=>1] = 1.0
        A[s=>2, a=>2, b=>2] = 1.0
        A[s=>2, a=>2, b=>3] = 1.0
        A[s=>3, a=>3, b=>2] = 1.0
        A[s=>3, a=>3, b=>3] = 1.0
        return (t::Index{Int64}, m::Index{Int64}, n::Index{Int64}) -> replaceinds(A, s=>t, a=>m, b=>n)
    elseif type == "left"
        s = Index(3)
        b = Index(3)
        A = ITensor(s,b)
        A[s=>1, b=>1] = 1.0
        A[s=>1, b=>2] = 1.0
        A[s=>2, b=>1] = 1.0
        A[s=>2, b=>2] = 1.0
        A[s=>2, b=>3] = 1.0
        A[s=>3, b=>2] = 1.0
        A[s=>3, b=>3] = 1.0
        return (t::Index{Int64}, n::Index{Int64}) -> replaceinds(A, s=>t, b=>n)
    elseif type == "right"
        s = Index(3)
        a = Index(3)
        A = ITensor(s,a)
        A[s=>1, a=>1] = 1.0
        A[s=>2, a=>2] = 1.0
        A[s=>3, a=>3] = 1.0
        return (t::Index{Int64}, m::Index{Int64}) -> replaceinds(A, s=>t, a=>m)
    else
        error("Unknown type $type for A tensor")
    end
end

function tensor_types(L; symmetry_breaking_sites=Int[])
    types = String[]
    counter = 1
    for i in 1:L
        if i == 1 && i+1 ∈ symmetry_breaking_sites
            #   |
            # _ .
            push!(types, "delta")
        elseif i == 1 && i+1 ∉ symmetry_breaking_sites
            # _ .
            push!(types, "left")
        elseif i == L && i ∈ symmetry_breaking_sites
            # |
            # . _
            push!(types, "delta")
            push!(types, "delta")
        elseif i == L && i ∉ symmetry_breaking_sites
            # . _
            push!(types, "right")
        elseif i ∈ symmetry_breaking_sites && i+1 ∈ symmetry_breaking_sites
            # |   |
            # . _ .
            push!(types, "delta")
            push!(types, "delta")
        elseif i ∈ symmetry_breaking_sites && i+1 ∉ symmetry_breaking_sites
            #  |
            #  . _ .
            push!(types, "delta")
            push!(types, "left")
        elseif i ∉ symmetry_breaking_sites && i+1 ∈ symmetry_breaking_sites
            #      |
            #  . _ .
            push!(types, "right")
        else

            #     
            # . _ .
            push!(types, "normal")
        end
    end

    return types
end



function constraint_projector(L, sites; symmetry_breaking_sites=Int[])
    if length(symmetry_breaking_sites) > 0
        println("Gauge breaking not supported yet ! ")
        return
    end
    N = length(sites)

    left = get_A("left")
    normal = get_A("normal")
    right = get_A("right")

    sites_primed = sites'

    links = [Index(3) for _ in 1:N-1]
    auxiliary_inds = [Index(3, "Auxiliary,n=$i") for i = 1:N]
    mpo = MPO(N)
    types = tensor_types(L; symmetry_breaking_sites=symmetry_breaking_sites)
    for i in 1:N
        if types[i] == "normal"
            mpo[i] = normal(auxiliary_inds[i], links[i-1], links[i]) * δ(auxiliary_inds[i], sites_primed[i], sites[i])
        elseif types[i] == "left"
            mpo[i] = left(auxiliary_inds[i], links[i]) * δ(auxiliary_inds[i], sites_primed[i], sites[i])
        elseif types[i] == "right"
            mpo[i] = right(auxiliary_inds[i], links[i-1]) * δ(auxiliary_inds[i], sites_primed[i], sites[i])
        elseif types[i] == "delta"
            mpo[i] = δ(sites_primed[i], sites[i])
        end
    end
            
    # for m = 2:N-1
    #     mpo[m] = normal(auxiliary_inds[m], links[m-1], links[m]) * δ(auxiliary_inds[m], sites_primed[m], sites[m])
    # end
    # mpo[N] = right(auxiliary_inds[N], links[N-1]) * δ(auxiliary_inds[N], sites_primed[N], sites[N])

    return mpo
end          