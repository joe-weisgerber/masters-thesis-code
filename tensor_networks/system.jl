using ITensors, ITensorMPS

"""Implements the creation of the sites and of the MPS of the scar state"""

function get_sites(s, L; symmetry_breaking_sites=[])
    sites = siteinds("S=1", L+length(symmetry_breaking_sites))
    return sites
end 

function entry_to_state(x)
    state = String[]
    for i in x
        if i == -1
            push!(state, "Dn")
        elseif i == 0
            push!(state, "Z0")
        elseif i == 1
            push!(state, "Up")
        else
            error("Invalid entry $i in basis state")
        end
    end
    return state
end

function create_scar_state_pos(L, i, sites; symmetry_breaking_sites=[])
    A = Vector{ITensor}(undef, L+length(symmetry_breaking_sites))
    #println(length(A))
    # For now we assume that site 1 is not symmetry broken.
    # This means that the boundaries are always of reduced form
    links = []
    counter = 1
    #println("Sites : ", length(sites))
    for k in 1:L-1
        if k%2 == 1
            if k ∈ symmetry_breaking_sites
                # Vertical site in between bricks
                push!(links, Index(2, "link,$counter"))
                counter += 1
                #println("Added link with bond dimension 2")
                # Start of brick -> horizontal site itself 
                # is 2x2 if former site is symmetry breaking
                push!(links, Index(2, "link,$counter"))
                counter += 1
                #println("Added link with bond dimension 2")
            else
                # No symmetry breaking within bricks
                # -> horizontal site is 1x2
                push!(links, Index(2, "link,$counter"))
                counter += 1
                #println("Added link with bond dimension 2")
            end
        else
            if k ∈ symmetry_breaking_sites
                # Vertical site within brick
                push!(links, Index(2, "link,$counter"))
                counter += 1
                #println("Added link with bond dimension 2")
            end
            if k+1 ∈ symmetry_breaking_sites
                # End of brick -> horizontal site is 2x2 if former site is symmetry breaking
                push!(links, Index(2, "link,$counter"))
                counter += 1
                #println("Added link with bond dimension 2")
            else
                push!(links, Index(1, "link,$counter"))
                counter += 1
                #println("Added link with bond dimension 1")
            end
        end
    end

    if L ∈ symmetry_breaking_sites
        push!(links, Index(2, "link,$counter"))
        #println("Added link with bond dimension 2")
    end

    if L+1 ∈ symmetry_breaking_sites
        push!(links, Index(2, "link,$counter"))
        #println("Added link with bond dimension 2")
    end

    #println(length(links))
    
    if 1 ∉ symmetry_breaking_sites
        A[1] = ITensor(sites[1], links[1])
        A[1][sites[1]=>1, links[1]=>2] =  1.0/sqrt(2)
        A[1][sites[1]=>2, links[1]=>1] = 1.0/sqrt(2)
        #println("Counter : ", 1, " Site : ", 1, " Adding initial tensor")
        counter = 2
    else
        A[1] = ITensor(sites[1], links[1])
        A[1][sites[1]=>3, links[1]=>2] =  1.0
        A[1][sites[1]=>2, links[1]=>1] = 1.0

        A[2] = ITensor(sites[2], links[1], links[2])
        A[2][sites[2]=>1, links[1]=>2, links[2]=>2] = 1.0/sqrt(2)
        A[2][sites[2]=>2, links[1]=>1, links[2]=>1] = 1.0/sqrt(2)
        counter = 3
    end

    for i in 2:L-1
        if i%2 == 0

            if i ∈ symmetry_breaking_sites
                # Vertical site within brick
                #println("Counter : ", counter, " Site : ", i, " Adding vertical link within brick")
                A[counter] = ITensor(sites[counter], links[counter-1], links[counter])
                A[counter][sites[counter]=>1, links[counter-1]=>2, links[counter]=>2] = 1.0
                A[counter][sites[counter]=>3, links[counter-1]=>1, links[counter]=>1] = 1.0
                counter += 1
            end
            if i+1 ∈ symmetry_breaking_sites
                #println("Counter : ", counter, " Site : ", i, " Adding extended tensor at end of brick")
                # End of brick -> Need extended tensors
                A[counter] = ITensor(sites[counter], links[counter-1], links[counter])
                A[counter][sites[counter]=>1, links[counter-1]=>1, links[counter]=>2] = -1.0
                A[counter][sites[counter]=>2, links[counter-1]=>2, links[counter]=>1] = 1.0
                counter += 1
            else
                # Use the standard tensor
                #println("Counter : ", counter, " Site : ", i, " Adding standard tensor at end of brick")
                A[counter] = ITensor(sites[counter], links[counter-1], links[counter])
                A[counter][sites[counter]=>1, links[counter-1]=>1, links[counter]=>1] = -1.0
                A[counter][sites[counter]=>2, links[counter-1]=>2, links[counter]=>1] = 1.0
                counter += 1
            end
        else
            if i ∈ symmetry_breaking_sites
                # Vertical site in between bricks
                #println("Counter : ", counter, " Site : ", i, " Adding vertical link in between bricks")
                A[counter] = ITensor(sites[counter], links[counter-1], links[counter])
                A[counter][sites[counter]=>2, links[counter-1]=>1, links[counter]=>1] = 1.0
                A[counter][sites[counter]=>2, links[counter-1]=>2, links[counter]=>2] = 1.0
                A[counter][sites[counter]=>1, links[counter-1]=>2, links[counter]=>1] = 1.0
                A[counter][sites[counter]=>3, links[counter-1]=>1, links[counter]=>2] = 1.0
                counter += 1
                # Next brick is to the right of symmetry breaking + begin of new block
                # Choose extended site tensor
                #println("Counter : ", counter, " Site : ", i, " Adding extended tensor at start of new brick")
                A[counter] = ITensor(sites[counter], links[counter-1], links[counter])
                A[counter][sites[counter]=>1, links[counter-1]=>2, links[counter]=>2] = 1.0/sqrt(2)
                A[counter][sites[counter]=>2, links[counter-1]=>1, links[counter]=>1] = 1.0/sqrt(2)
                counter += 1
            else
                # Standard tensor for new block
                #println("Counter : ", counter, " Site : ", i, " Adding standard tensor at start of new brick")
                A[counter] = ITensor(sites[counter], links[counter-1], links[counter])
                A[counter][sites[counter]=>1, links[counter-1]=>1, links[counter]=>2] = 1.0/sqrt(2)
                A[counter][sites[counter]=>2, links[counter-1]=>1, links[counter]=>1] = 1.0/sqrt(2)
                counter += 1
            end
        end
    end

  
    if L ∈ symmetry_breaking_sites
        #println("Counter : ", counter, " Site : ", L, " Adding vertical link within brick")
        A[counter] = ITensor(sites[counter], links[counter-1], links[counter])
        A[counter][sites[counter]=>1, links[counter-1]=>2, links[counter]=>2] = 1.0
        A[counter][sites[counter]=>3, links[counter-1]=>1, links[counter]=>1] = 1.0
        counter += 1
    end

    if L+1 ∉ symmetry_breaking_sites
    #println("Counter : ", counter, " Site : ", L, " Adding final tensor")
        A[counter] = ITensor(sites[counter], links[counter-1])
        A[counter][sites[counter]=>1, links[counter-1]=>1] = -1.0
        A[counter][sites[counter]=>2, links[counter-1]=>2] = 1.0
    else
        A[counter] = ITensor(sites[counter], links[counter-1], links[counter])
        A[counter][sites[counter]=>1, links[counter-1]=>1, links[counter]=>2] = -1.0
        A[counter][sites[counter]=>2, links[counter-1]=>2, links[counter]=>1] = 1.0

        counter += 1
        A[counter] = ITensor(sites[counter], links[counter-1])
        A[counter][sites[counter]=>1, links[counter-1]=>2] = 1.0
        A[counter][sites[counter]=>2, links[counter-1]=>1] = 1.0
    end

    # A[1] = ITensor(sites[1], links[1])
    # A[1][sites[1]=>1, links[1]=>2] =  1.0/sqrt(2)
    # A[1][sites[1]=>2, links[1]=>1] = 1.0/sqrt(2)
    # for n in 2:L-1
    #     if n%2 == 1
    #         A[n] = ITensor(sites[n], links[n-1], links[n])
    #         A[n][sites[n]=>1, links[n-1]=>1, links[n]=>2] = 1.0/sqrt(2)
    #         A[n][sites[n]=>2, links[n-1]=>1, links[n]=>1] = 1.0/sqrt(2)
    #     else
    #         A[n] = ITensor(sites[n], links[n-1], links[n])
    #         A[n][sites[n]=>1, links[n-1]=>1, links[n]=>1] = -1.0
    #         A[n][sites[n]=>2, links[n-1]=>2, links[n]=>1] = 1.0
    #     end
    # end

    # A[L] = ITensor(sites[L], links[L-1])
    # A[L][sites[L]=>1, links[L-1]=>1] = -1.0
    # A[L][sites[L]=>2, links[L-1]=>2] = 1.0

    psi = MPS(A)
    return psi
end


function create_scar_state_neg(L, i, sites; symmetry_breaking_sites=[])
    A = Vector{ITensor}(undef, L+length(symmetry_breaking_sites))
    #println(length(A))
    # For now we assume that site 1 is not symmetry broken.
    # This means that the boundaries are always of reduced form
    links = []
    counter = 1
    #println("Sites : ", length(sites))
    for k in 1:L-1
        if k%2 == 1
            if k ∈ symmetry_breaking_sites
                # Vertical site in between bricks
                push!(links, Index(2, "link,$counter"))
                counter += 1
                #println("Added link with bond dimension 2")
                # Start of brick -> horizontal site itself 
                # is 2x2 if former site is symmetry breaking
                push!(links, Index(2, "link,$counter"))
                counter += 1
                #println("Added link with bond dimension 2")
            else
                # No symmetry breaking within bricks
                # -> horizontal site is 1x2
                push!(links, Index(2, "link,$counter"))
                counter += 1
                #println("Added link with bond dimension 2")
            end
        else
            if k ∈ symmetry_breaking_sites
                # Vertical site within brick
                push!(links, Index(2, "link,$counter"))
                counter += 1
                #println("Added link with bond dimension 2")
            end
            if k+1 ∈ symmetry_breaking_sites
                # End of brick -> horizontal site is 2x2 if former site is symmetry breaking
                push!(links, Index(2, "link,$counter"))
                counter += 1
                #println("Added link with bond dimension 2")
            else
                push!(links, Index(1, "link,$counter"))
                counter += 1
                #println("Added link with bond dimension 1")
            end
        end
    end

    if L ∈ symmetry_breaking_sites
        push!(links, Index(2, "link,$counter"))
        #println("Added link with bond dimension 2")
    end

    if L+1 ∈ symmetry_breaking_sites
        push!(links, Index(2, "link,$counter"))
        #println("Added link with bond dimension 2")
    end

    #println(length(links))
    
    if 1 ∉ symmetry_breaking_sites
        A[1] = ITensor(sites[1], links[1])
        A[1][sites[1]=>3, links[1]=>2] =  1.0/sqrt(2)
        A[1][sites[1]=>2, links[1]=>1] = 1.0/sqrt(2)
        #println("Counter : ", 1, " Site : ", 1, " Adding initial tensor")
        counter = 2
    else
        A[1] = ITensor(sites[1], links[1])
        A[1][sites[1]=>1, links[1]=>2] =  1.0
        A[1][sites[1]=>2, links[1]=>1] = 1.0

        A[2] = ITensor(sites[2], links[1], links[2])
        A[2][sites[2]=>3, links[1]=>2, links[2]=>2] = 1.0/sqrt(2)
        A[2][sites[2]=>2, links[1]=>1, links[2]=>1] = 1.0/sqrt(2)
        counter = 3
    end

    for i in 2:L-1
        if i%2 == 0

            if i ∈ symmetry_breaking_sites
                # Vertical site within brick
                #println("Counter : ", counter, " Site : ", i, " Adding vertical link within brick")
                A[counter] = ITensor(sites[counter], links[counter-1], links[counter])
                A[counter][sites[counter]=>3, links[counter-1]=>2, links[counter]=>2] = 1.0
                A[counter][sites[counter]=>1, links[counter-1]=>1, links[counter]=>1] = 1.0
                counter += 1
            end
            if i+1 ∈ symmetry_breaking_sites
                #println("Counter : ", counter, " Site : ", i, " Adding extended tensor at end of brick")
                # End of brick -> Need extended tensors
                A[counter] = ITensor(sites[counter], links[counter-1], links[counter])
                A[counter][sites[counter]=>3, links[counter-1]=>1, links[counter]=>2] = -1.0
                A[counter][sites[counter]=>2, links[counter-1]=>2, links[counter]=>1] = 1.0
                counter += 1
            else
                # Use the standard tensor
                #println("Counter : ", counter, " Site : ", i, " Adding standard tensor at end of brick")
                A[counter] = ITensor(sites[counter], links[counter-1], links[counter])
                A[counter][sites[counter]=>3, links[counter-1]=>1, links[counter]=>1] = -1.0
                A[counter][sites[counter]=>2, links[counter-1]=>2, links[counter]=>1] = 1.0
                counter += 1
            end
        else
            if i ∈ symmetry_breaking_sites
                # Vertical site in between bricks
                #println("Counter : ", counter, " Site : ", i, " Adding vertical link in between bricks")
                A[counter] = ITensor(sites[counter], links[counter-1], links[counter])
                A[counter][sites[counter]=>2, links[counter-1]=>1, links[counter]=>1] = 1.0
                A[counter][sites[counter]=>2, links[counter-1]=>2, links[counter]=>2] = 1.0
                A[counter][sites[counter]=>3, links[counter-1]=>2, links[counter]=>1] = 1.0
                A[counter][sites[counter]=>1, links[counter-1]=>1, links[counter]=>2] = 1.0
                counter += 1
                # Next brick is to the right of symmetry breaking + begin of new block
                # Choose extended site tensor
                #println("Counter : ", counter, " Site : ", i, " Adding extended tensor at start of new brick")
                A[counter] = ITensor(sites[counter], links[counter-1], links[counter])
                A[counter][sites[counter]=>3, links[counter-1]=>2, links[counter]=>2] = 1.0/sqrt(2)
                A[counter][sites[counter]=>2, links[counter-1]=>1, links[counter]=>1] = 1.0/sqrt(2)
                counter += 1
            else
                # Standard tensor for new block
                #println("Counter : ", counter, " Site : ", i, " Adding standard tensor at start of new brick")
                A[counter] = ITensor(sites[counter], links[counter-1], links[counter])
                A[counter][sites[counter]=>3, links[counter-1]=>1, links[counter]=>2] = 1.0/sqrt(2)
                A[counter][sites[counter]=>2, links[counter-1]=>1, links[counter]=>1] = 1.0/sqrt(2)
                counter += 1
            end
        end
    end

  
    if L ∈ symmetry_breaking_sites
        #println("Counter : ", counter, " Site : ", L, " Adding vertical link within brick")
        A[counter] = ITensor(sites[counter], links[counter-1], links[counter])
        A[counter][sites[counter]=>3, links[counter-1]=>2, links[counter]=>2] = 1.0
        A[counter][sites[counter]=>1, links[counter-1]=>1, links[counter]=>1] = 1.0
        counter += 1
    end

    if L+1 ∉ symmetry_breaking_sites
    #println("Counter : ", counter, " Site : ", L, " Adding final tensor")
        A[counter] = ITensor(sites[counter], links[counter-1])
        A[counter][sites[counter]=>3, links[counter-1]=>1] = -1.0
        A[counter][sites[counter]=>2, links[counter-1]=>2] = 1.0
    else
        A[counter] = ITensor(sites[counter], links[counter-1], links[counter])
        A[counter][sites[counter]=>3, links[counter-1]=>1, links[counter]=>2] = -1.0
        A[counter][sites[counter]=>2, links[counter-1]=>2, links[counter]=>1] = 1.0

        counter += 1
        A[counter] = ITensor(sites[counter], links[counter-1])
        A[counter][sites[counter]=>3, links[counter-1]=>2] = 1.0
        A[counter][sites[counter]=>2, links[counter-1]=>1] = 1.0
    end

    # A[1] = ITensor(sites[1], links[1])
    # A[1][sites[1]=>1, links[1]=>2] =  1.0/sqrt(2)
    # A[1][sites[1]=>2, links[1]=>1] = 1.0/sqrt(2)
    # for n in 2:L-1
    #     if n%2 == 1
    #         A[n] = ITensor(sites[n], links[n-1], links[n])
    #         A[n][sites[n]=>1, links[n-1]=>1, links[n]=>2] = 1.0/sqrt(2)
    #         A[n][sites[n]=>2, links[n-1]=>1, links[n]=>1] = 1.0/sqrt(2)
    #     else
    #         A[n] = ITensor(sites[n], links[n-1], links[n])
    #         A[n][sites[n]=>1, links[n-1]=>1, links[n]=>1] = -1.0
    #         A[n][sites[n]=>2, links[n-1]=>2, links[n]=>1] = 1.0
    #     end
    # end

    # A[L] = ITensor(sites[L], links[L-1])
    # A[L][sites[L]=>1, links[L-1]=>1] = -1.0
    # A[L][sites[L]=>2, links[L-1]=>2] = 1.0

    psi = MPS(A)
    return psi
end

function create_scar_state(L, i, sites; symmetry_breaking_sites=[])
    # i is ignored here and would only matter for higher spins
    scar_pos = create_scar_state_pos(L, i, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    scar_neg = create_scar_state_neg(L, i, sites; symmetry_breaking_sites=symmetry_breaking_sites)
    scar_state = scar_pos + scar_neg
    normalize!(scar_state)
    return scar_state
end

function create_revival_state(L, sites; gauge_breaking_sites=[])
    entries = []
    symmetry_breaking_sites = [i for i in 1:L if i in gauge_breaking_sites || mod1(i-1, L) in gauge_breaking_sites]
    for i in 1:L
        if i ∈ symmetry_breaking_sites
            push!(entries, "Z0")
        end
        push!(entries, "Up")
    end
    push!(entries, "Z0")  
    println(entries)
    return MPS(sites, entries)
end

function get_symmetry_breaking_sites(L, gauge_breaking_sites)
    symmetry_breaking_sites = [i for i in 1:L if i in gauge_breaking_sites || i-1 in gauge_breaking_sites]
    return symmetry_breaking_sites
end