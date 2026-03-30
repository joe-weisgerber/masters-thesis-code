using ITensors, ITensorMPS, ITensorMPOConstruction

"""Implements the hamiltonian"""

ITensors.op(::OpName"P01",::SiteType"S=1") =
 [1 0 0
  0 1 0
  0 0 0]
ITensors.op(::OpName"P-10",::SiteType"S=1") =
 [0 0 0
  0 1 0
  0 0 1]
ITensors.op(::OpName"U01",::SiteType"S=1") =
 [0 1 0
  1 0 0
  0 0 0]
ITensors.op(::OpName"U0-1",::SiteType"S=1") =
 [0 0 0
  0 0 1
  0 1 0]
ITensors.op(::OpName"U",::SiteType"S=1") =
 [0 1 0
  0 0 1
  0 0 0]
ITensors.op(::OpName"Ud",::SiteType"S=1") =
 [0 0 0
  1 0 0
  0 1 0]
ITensors.op(::OpName"01",::SiteType"S=1") =
 [0 0 0
  1 0 0
  0 0 0]
ITensors.op(::OpName"10",::SiteType"S=1") =
 [0 1 0
  0 0 0
  0 0 0]
ITensors.op(::OpName"0-1",::SiteType"S=1") =
 [0 0 0
  0 0 1
  0 0 0]
ITensors.op(::OpName"-10",::SiteType"S=1") =
 [0 0 0
  0 0 0
  0 1 0]
ITensors.op(::OpName"E2",::SiteType"S=1") =
 [1 0 0
  0 0 0
  0 0 1]
ITensors.op(::OpName"-1-1", ::SiteType"S=1") =
 [0 0 0
  0 0 0
  0 0 1]
ITensors.op(::OpName"11", ::SiteType"S=1") =
 [1 0 0
  0 0 0
  0 0 0]
ITensors.op(::OpName"00", ::SiteType"S=1") = 
    [0 0 0
     0 1 0
     0 0 0]
ITensors.op(::OpName"E", ::SiteType"S=1") = 
    [1 0 0
     0 0 0
     0 0 -1]

ITensors.op(::OpName"Zero", ::SiteType"S=1") = 
    [0 0 0
     0 0 0
     0 0 0]


function hamiltonian(s, L, sites; PBC=false, potential=false, g=false,
    gauge_breaking_sites=[], gauge_breaking_strengths=[], orientation="vertical", 
    constraint_breaking=false)
    symmetry_breaking_sites = [i for i in 1:L+1 if i in gauge_breaking_sites || i-1 in gauge_breaking_sites]
    #println("Symmetry breaking sites : ", symmetry_breaking_sites)
    os = OpSum()

    # Start with the kinetic and potential perturbation_terms
    # Counter always points to the horizontal link

    # The left boundary term is treated differently
    if 1 ∉ symmetry_breaking_sites && 2 ∉ symmetry_breaking_sites
        os += 1, "U01", 1, "P01", 2
        os += 1, "U0-1", 1, "P-10", 2
        if potential != false
            os += 3potential, "E2", 1
        end
        if g != false
            os += -g, "E", 1
        end
        if constraint_breaking != false
            os += constraint_breaking, "U", 1
            os += constraint_breaking, "Ud", 1
        end
        counter = 2
        #println("Added U01 P01 at site 1 and counter 1")
    elseif 1 ∉ symmetry_breaking_sites && 2 ∈ symmetry_breaking_sites
        os += 1, "U", 1, "U", 2
        os += 1, "Ud", 1, "Ud", 2
        if potential != false
            os += 3potential, "E2", 1
        end
        if g != false
            os += -g, "E", 1
        end
        if constraint_breaking != false
            prinln("Not supported yet")
        end
        counter = 3
        #println("Added U U at site 1 and counter 1")
    elseif 1 ∈ symmetry_breaking_sites && 2 ∉ symmetry_breaking_sites
        println("Not supported yet")
    else
        os += 1, "U", 1, "Ud", 2, "Ud", 3
        os += 1, "Ud", 1, "U", 2, "U", 3
        if potential != false
            os += potential, "E2", 1
            os += 2potential, "E2", 2
        end
        if g != false
            os += g, "E", 1
        end
        if constraint_breaking != false
            println("Not supported yet")
        end
        counter = 4
        #println("Added U Ud Ud at site 1 and counter 2")
    end

    # All the middle sites

    for i in 2:L-1
        l = counter-1
        r = counter+1

        if i ∈ symmetry_breaking_sites && i+1 ∈ symmetry_breaking_sites
            os += 1, "U", l, "Ud", counter, "Ud", r
            os += 1, "Ud", l, "U", counter, "U", r
            if potential != false
                os += potential, "E2", l
                os += 2potential, "E2", counter
                os += potential, "E2", r
            end
            if g != false
                os += (-1)^(i-1)*g, "E", l
            end
            if constraint_breaking != false
                println("Not supported yet")
            end
            #println("Added U Ud Ud at site ", i, " and counter ", counter)
            counter += 2
            
        elseif i ∈ symmetry_breaking_sites && i+1 ∉ symmetry_breaking_sites
            os += 1, "Ud", l, "0-1", counter, "P-10", r
            os += 1, "U", l, "-10", counter, "P-10", r
            os += 1, "U", l, "01", counter, "P01", r
            os += 1, "Ud", l, "10", counter, "P01", r
            if potential != false
                os += 2potential, "E2", counter
                os += potential, "E2", l
            end
            if g != false
                os += (-1)^(i-1)*g, "E", l
            end
            if constraint_breaking != false
                println("Not supported yet")
            end
            #println("Added Ud 0-1 P-10 at site ", i, " and counter ", counter)
            counter += 1
        elseif i ∉ symmetry_breaking_sites && i+1 ∈ symmetry_breaking_sites
            os += 1, "P-10", l, "0-1", counter, "U", r
            os += 1, "P-10", l, "-10", counter, "Ud", r
            os += 1, "P01", l, "01", counter, "Ud", r
            os += 1, "P01", l, "10", counter, "U", r
            if potential != false
                os += 2potential, "E2", counter
                os += potential, "00", l, "11", counter
                os += potential, "11", l, "00", counter
                os += potential, "-1-1", l, "00", counter
                os += potential, "00", l, "-1-1", counter
            end
            if g != false
                os += (-1)^(i-1)*g, "E", l, "00", counter
                os += -(-1)^(i-1)*g, "00", l, "E", counter
            end
            if constraint_breaking != false
                println("Not supported yet")
            end
            #println("Added P-10 0-1 U at site ", i, " and counter ", counter)
            counter += 2
        else
            os += 1, "P01", l, "U01", counter, "P01", r
            os += 1, "P-10", l, "U0-1", counter, "P-10", r
            if potential != false
                os += 2potential, "E2", counter
                os += potential, "00", l, "11", counter
                os += potential, "00", l, "-1-1", counter
                os += potential, "11", l, "00", counter
                os += potential, "-1-1", l, "00", counter
            end
            if g != false
                os += (-1)^(i-1)*g, "E", l, "00", counter
                os += -(-1)^(i-1)*g, "00", l, "E", counter
            end
            if constraint_breaking != false
                os += constraint_breaking, "U", counter
                os += constraint_breaking, "Ud", counter
            end
            #println("Added P01 U01 P01 at site ", i, " and counter ", counter)
            counter += 1
        end
    end

    # Right boundary is treated differently
    l = counter - 1
    r = counter + 1
    if L ∈ symmetry_breaking_sites && L+1 ∈ symmetry_breaking_sites
        os += 1, "U", l, "Ud", counter, "Ud", r
        os += 1, "Ud", l, "U", counter, "U", r
        if potential != false
            os += potential, "E2", l
            os += 2potential, "E2", counter
            os += potential, "E2", r
        end
        if g != false
            os += (-1)^(i-1)*g, "E", l
            os += (-1)^i*g, "E", r
        end
        if constraint_breaking != false
            println("Not supported yet")
        end
        
        #println("Added U Ud Ud at site ", L, " and counter ", counter)
    elseif L ∈ symmetry_breaking_sites && L+1 ∉ symmetry_breaking_sites  
        os += 1, "U", l, "Ud", counter
        os += 1, "Ud", l, "U", counter
        if potential != false
            os += 3potential, "E2", counter
            os += potential, "E2", l
        end  
        if g != false
            os += (-1)^(L-1)*g, "E", l
            os += (-1)^L*g, "E", counter
        end
        if constraint_breaking != false
            println("Not supported yet")
        end
        #println("Added U Ud at site ", L, " and counter ", counter)
    elseif L ∉ symmetry_breaking_sites && L+1 ∈ symmetry_breaking_sites
        prinln("Not supported yet")    
    else
        os += 1, "P01", l, "U01", counter
        os += 1, "P-10", l, "U0-1", counter
        if potential != false
            os += 3potential, "E2", counter
            os += potential, "00", l, "11", counter
            os += potential, "11", l, "00", counter
            os += potential, "-1-1", l, "00", counter
            os += potential, "00", l, "-1-1", counter
        end 
        if g != false
            os += (-1)^(L-1)*g, "E", l, "00", counter
            os += -(-1)^(L-1)*g, "00", l, "E", counter
            os += g, "E", counter  # (-1)^0 = 1 but negative sign comes from horizontal being opposite of vertical
        end
        if constraint_breaking != false
            os += constraint_breaking, "U", counter
            os += constraint_breaking, "Ud", counter
        end
        #println("Added P01 U01 at site ", L, " and counter ", counter)
    end

    # Adding the gauge breaking terms
    if 1 ∈ symmetry_breaking_sites
        counter = 2
    else
        counter = 1
    end

    gauge_counter = 1

    for i in 1:L
        if i ∈ gauge_breaking_sites
            if orientation == "vertical"
                l = counter - 1
                r = counter + 1
                os += gauge_breaking_strengths[gauge_counter], "U", l, "Ud", r
                os += gauge_breaking_strengths[gauge_counter], "Ud", l, "U", r
                #println("Added gauge breaking term at site ", i, " and counter ", counter)
                counter += 2
                gauge_counter += 1
            elseif orientation == "horizontal"
                os += gauge_breaking_strengths[gauge_counter], "U", counter
                os += gauge_breaking_strengths[gauge_counter], "Ud", counter
                counter += 2
                gauge_counter += 1
            end

        else
            if i+1 ∈ symmetry_breaking_sites
                counter += 2
            else
                counter += 1
            end
        end
    end
    return MPO_new(os, sites) 
end


function projector_hamiltonian(s, L, sites; PBC=false, potential=false, g=0.0,
    gauge_breaking_sites=[], gauge_breaking_strengths=[])

    symmetry_breaking_sites = [i for i in 1:L if i in gauge_breaking_sites || mod1(i-1, L) in gauge_breaking_sites]
    os = OpSum()
    counter = 1
    gauge_counter = 1
    for i in 1:L
        l = counter-1
        r = counter+1

        if i ∈ symmetry_breaking_sites && i+1 ∈ symmetry_breaking_sites
            # |   |
            # . _ .
            os += 1, "Ud", l, "U", counter, "U", r
            os += 1, "U", l, "Ud", counter, "Ud", r
            counter += 2
        elseif i ∈ symmetry_breaking_sites
            # |
            # . _ .
            os += 1, "Ud", l, "U", counter
            os += 1, "U", l, "Ud", counter
            counter += 1
        elseif i+1 ∈ symmetry_breaking_sites
            #     |
            # . _ .
            os += 1, "U", counter, "U", r
            os += 1, "Ud", counter, "Ud", r
            counter += 2
        else
            # . _ .
            os += 1, "U", counter
            os += 1, "Ud", counter
            if potential != false
                if i == 1
                    os += 3potential, "E2", counter
                elseif i == L
                    os += 3potential, "E2", counter
                    os += potential, "11", l, "00", counter
                    os += potential, "00", l, "11", counter
                    os += potential, "-1-1", l, "00", counter
                    os += potential, "00", l, "-1-1", counter
                else
                    os += 2potential, "E2", counter
                    os += potential, "11", l, "00", counter
                    os += potential, "00", l, "11", counter
                    os += potential, "-1-1", l, "00", counter
                    os += potential, "00", l, "-1-1", counter
                end
            end
            counter += 1
        end

        if i ∈ gauge_breaking_sites
            os += gauge_breaking_strengths[gauge_counter], "U", l, "Ud", r
            os += gauge_breaking_strengths[gauge_counter], "Ud", l, "U", r
            gauge_counter += 1
        end

    end
    
    return MPO_new(os, sites)
end

function charge_conjugation(s, L, sites; gauge_breaking_sites=[])
    symmetry_breaking_sites = [i for i in 1:L+1 if i in gauge_breaking_sites || i-1 in gauge_breaking_sites]
    #println("Symmetry breaking sites : ", symmetry_breaking_sites)
    os = OpSum()

    # Start with the kinetic and potential perturbation_terms
    # Counter always points to the horizontal link

    # The left boundary term is treated differently
    if 1 ∉ symmetry_breaking_sites && 2 ∉ symmetry_breaking_sites
        os += -1, "E", 1
        counter = 2
        #println("Added U01 P01 at site 1 and counter 1")
    elseif 1 ∉ symmetry_breaking_sites && 2 ∈ symmetry_breaking_sites
        os += -1, "E", 1
        counter = 3
        #println("Added U U at site 1 and counter 1")
    elseif 1 ∈ symmetry_breaking_sites && 2 ∉ symmetry_breaking_sites
        println("Not supported yet")
    else
        os += 1, "E", 1
        counter = 4
        #println("Added U Ud Ud at site 1 and counter 2")
    end

    # All the middle sites

    for i in 2:L-1
        l = counter-1
        r = counter+1

        if i ∈ symmetry_breaking_sites && i+1 ∈ symmetry_breaking_sites
            os += (-1)^(i-1), "E", l
            #println("Added U Ud Ud at site ", i, " and counter ", counter)
            counter += 2
            
        elseif i ∈ symmetry_breaking_sites && i+1 ∉ symmetry_breaking_sites
            os += (-1)^(i-1), "E", l
            
            #println("Added Ud 0-1 P-10 at site ", i, " and counter ", counter)
            counter += 1
        elseif i ∉ symmetry_breaking_sites && i+1 ∈ symmetry_breaking_sites
            os += (-1)^(i-1), "E", l, "00", counter
            os += -(-1)^(i-1), "00", l, "E", counter
            #println("Added P-10 0-1 U at site ", i, " and counter ", counter)
            counter += 2
        else
            
            os += (-1)^(i-1), "E", l, "00", counter
            os += -(-1)^(i-1), "00", l, "E", counter
            #println("Added P01 U01 P01 at site ", i, " and counter ", counter)
            counter += 1
        end
    end

    # Right boundary is treated differently
    l = counter - 1
    r = counter + 1
    if L ∈ symmetry_breaking_sites && L+1 ∈ symmetry_breaking_sites
        os += (-1)^(L-1), "E", l
        os += (-1)^L, "E", r
        #println("Added U Ud Ud at site ", L, " and counter ", counter)
    elseif L ∈ symmetry_breaking_sites && L+1 ∉ symmetry_breaking_sites  
            os += (-1)^(L-1), "E", l
            os += (-1)^L, "E", counter
        #println("Added U Ud at site ", L, " and counter ", counter)
    elseif L ∉ symmetry_breaking_sites && L+1 ∈ symmetry_breaking_sites
        prinln("Not supported yet")    
    else
        os += (-1)^(L-1), "E", l, "00", counter
        os += -(-1)^(L-1), "00", l, "E", counter
        os += 1, "E", counter  # (-1)^0 = 1 but negative sign comes from horizontal being opposite of vertical
        #println("Added P01 U01 at site ", L, " and counter ", counter)
    end

    return MPO_new(os, sites) 
end

function potential(s, L, sites; gauge_breaking_sites=[])
    symmetry_breaking_sites = [i for i in 1:L+1 if i in gauge_breaking_sites || i-1 in gauge_breaking_sites]
    #println("Symmetry breaking sites : ", symmetry_breaking_sites)
    os = OpSum()

    # Start with the kinetic and potential perturbation_terms
    # Counter always points to the horizontal link

    # The left boundary term is treated differently
    if 1 ∉ symmetry_breaking_sites && 2 ∉ symmetry_breaking_sites
        os += 3, "E2", 1
        counter = 2
        #println("Added U01 P01 at site 1 and counter 1")
    elseif 1 ∉ symmetry_breaking_sites && 2 ∈ symmetry_breaking_sites
        os += 3, "E2", 1
        counter = 3
        #println("Added U U at site 1 and counter 1")
    elseif 1 ∈ symmetry_breaking_sites && 2 ∉ symmetry_breaking_sites
        println("Not supported yet")
    else 
        os += 1, "E2", 1
        os += 2, "E2", 2
        counter = 4
        #println("Added U Ud Ud at site 1 and counter 2")
    end

    # All the middle sites

    for i in 2:L-1
        l = counter-1
        r = counter+1

        if i ∈ symmetry_breaking_sites && i+1 ∈ symmetry_breaking_sites  
            os += 1, "E2", l
            os += 2, "E2", counter
            os += 1, "E2", r
            
            #println("Added U Ud Ud at site ", i, " and counter ", counter)
            counter += 2
            
        elseif i ∈ symmetry_breaking_sites && i+1 ∉ symmetry_breaking_sites
            
            os += 2, "E2", counter
            os += 1, "E2", l
            #println("Added Ud 0-1 P-10 at site ", i, " and counter ", counter)
            counter += 1
        elseif i ∉ symmetry_breaking_sites && i+1 ∈ symmetry_breaking_sites
            
            os += 2, "E2", counter
            os += 1, "00", l, "11", counter
            os += 1, "11", l, "00", counter
            os += 1, "-1-1", l, "00", counter
            os += 1, "00", l, "-1-1", counter
            #println("Added P-10 0-1 U at site ", i, " and counter ", counter)
            counter += 2
        else
            os += 2, "E2", counter
            os += 1, "00", l, "11", counter
            os += 1, "00", l, "-1-1", counter
            os += 1, "11", l, "00", counter
            os += 1, "-1-1", l, "00", counter
            #println("Added P01 U01 P01 at site ", i, " and counter ", counter)
            counter += 1
        end
    end

    # Right boundary is treated differently
    l = counter - 1
    r = counter + 1
    if L ∈ symmetry_breaking_sites && L+1 ∈ symmetry_breaking_sites
        
        os += 1, "E2", l
        os += 2, "E2", counter
        os += 1, "E2", r
         
        #println("Added U Ud Ud at site ", L, " and counter ", counter)
    elseif L ∈ symmetry_breaking_sites && L+1 ∉ symmetry_breaking_sites  
        os += 3, "E2", counter
        os += 1, "E2", l
        
        #println("Added U Ud at site ", L, " and counter ", counter)
    elseif L ∉ symmetry_breaking_sites && L+1 ∈ symmetry_breaking_sites
        println("Not supported yet")    
    else
        os += 3, "E2", counter
        os += 1, "00", l, "11", counter
        os += 1, "11", l, "00", counter
        os += 1, "-1-1", l, "00", counter
        os += 1, "00", l, "-1-1", counter
        
        #println("Added P01 U01 at site ", L, " and counter ", counter)
    end

    return MPO_new(os, sites) 
end


function gauge_breaking(s, L, sites; gauge_breaking_sites=[])
    symmetry_breaking_sites = [i for i in 1:L+1 if i in gauge_breaking_sites || i-1 in gauge_breaking_sites]
    #println("Symmetry breaking sites : ", symmetry_breaking_sites)
    os = OpSum()
    
    # Adding the gauge breaking terms
    if 1 ∈ symmetry_breaking_sites
        counter = 2
    else
        counter = 1
    end

    gauge_counter = 1

    for i in 1:L
        if i ∈ gauge_breaking_sites
            l = counter - 1
            r = counter + 1
            os += 1, "U", l, "Ud", r
            os += 1, "Ud", l, "U", r
            #println("Added gauge breaking term at site ", i, " and counter ", counter)
            counter += 2
            gauge_counter += 1

        else
            os += 1, "Zero", counter  # Seems to be required for things to work
            if i+1 ∈ symmetry_breaking_sites
                counter += 2
            else
                counter += 1
            end
        end
    end

    return MPO_new(os, sites) 
end

function hamiltonian_0(s, L, sites; gauge_breaking_sites=[])
    symmetry_breaking_sites = [i for i in 1:L+1 if i in gauge_breaking_sites || i-1 in gauge_breaking_sites]
    #println("Symmetry breaking sites : ", symmetry_breaking_sites)
    os = OpSum()

    # Start with the kinetic and potential perturbation_terms
    # Counter always points to the horizontal link

    # The left boundary term is treated differently
    if 1 ∉ symmetry_breaking_sites && 2 ∉ symmetry_breaking_sites
        os += 1, "U01", 1, "P01", 2
        os += 1, "U0-1", 1, "P-10", 2
        
        counter = 2
        #println("Added U01 P01 at site 1 and counter 1")
    elseif 1 ∉ symmetry_breaking_sites && 2 ∈ symmetry_breaking_sites
        os += 1, "U", 1, "U", 2
        os += 1, "Ud", 1, "Ud", 2
        
        counter = 3
        #println("Added U U at site 1 and counter 1")
    elseif 1 ∈ symmetry_breaking_sites && 2 ∉ symmetry_breaking_sites
        println("Not supported yet")
    else
        os += 1, "U", 1, "Ud", 2, "Ud", 3
        os += 1, "Ud", 1, "U", 2, "U", 3
        
        counter = 4
        #println("Added U Ud Ud at site 1 and counter 2")
    end

    # All the middle sites

    for i in 2:L-1
        l = counter-1
        r = counter+1

        if i ∈ symmetry_breaking_sites && i+1 ∈ symmetry_breaking_sites
            os += 1, "U", l, "Ud", counter, "Ud", r
            os += 1, "Ud", l, "U", counter, "U", r
            
            #println("Added U Ud Ud at site ", i, " and counter ", counter)
            counter += 2
            
        elseif i ∈ symmetry_breaking_sites && i+1 ∉ symmetry_breaking_sites
            os += 1, "Ud", l, "0-1", counter, "P-10", r
            os += 1, "U", l, "-10", counter, "P-10", r
            os += 1, "U", l, "01", counter, "P01", r
            os += 1, "Ud", l, "10", counter, "P01", r
            
            #println("Added Ud 0-1 P-10 at site ", i, " and counter ", counter)
            counter += 1
        elseif i ∉ symmetry_breaking_sites && i+1 ∈ symmetry_breaking_sites
            os += 1, "P-10", l, "0-1", counter, "U", r
            os += 1, "P-10", l, "-10", counter, "Ud", r
            os += 1, "P01", l, "01", counter, "Ud", r
            os += 1, "P01", l, "10", counter, "U", r
            
            #println("Added P-10 0-1 U at site ", i, " and counter ", counter)
            counter += 2
        else
            os += 1, "P01", l, "U01", counter, "P01", r
            os += 1, "P-10", l, "U0-1", counter, "P-10", r
            
            #println("Added P01 U01 P01 at site ", i, " and counter ", counter)
            counter += 1
        end
    end

    # Right boundary is treated differently
    l = counter - 1
    r = counter + 1
    if L ∈ symmetry_breaking_sites && L+1 ∈ symmetry_breaking_sites
        os += 1, "U", l, "Ud", counter, "Ud", r
        os += 1, "Ud", l, "U", counter, "U", r
        
        #println("Added U Ud Ud at site ", L, " and counter ", counter)
    elseif L ∈ symmetry_breaking_sites && L+1 ∉ symmetry_breaking_sites  
        os += 1, "U", l, "Ud", counter
        os += 1, "Ud", l, "U", counter
        
        #println("Added U Ud at site ", L, " and counter ", counter)
    elseif L ∉ symmetry_breaking_sites && L+1 ∈ symmetry_breaking_sites
        prinln("Not supported yet")    
    else
        os += 1, "P01", l, "U01", counter
        os += 1, "P-10", l, "U0-1", counter
        
        #println("Added P01 U01 at site ", L, " and counter ", counter)
    end

    # Adding the gauge breaking terms
    
    return MPO_new(os, sites) 
end

function hamiltonian_middle_perturbations(s, L, sites; type="potential", strength=0)
    os = OpSum()

    # Start with the kinetic and potential perturbation_terms
    # Counter always points to the horizontal link

    os += 1, "U01", 1, "P01", 2
    os += 1, "U0-1", 1, "P-10", 2
        
    for i in 2:L-1
        l = i-1
        r = i+1

        os += 1, "P01", l, "U01", i, "P01", r
        os += 1, "P-10", l, "U0-1", i, "P-10", r
        
        #println("Added P01 U01 P01 at site ", i, " and counter ", counter)

    end

    # Right boundary is treated differently
    l = L - 1
       
    os += 1, "P01", l, "U01", L
    os += 1, "P-10", l, "U0-1", L

    m = L ÷ 2

    if type == "potential"
        os += 2strength, "E2", m
        os += strength, "00", m-1, "11", m
        os += strength, "00", m-1, "-1-1", m
        os += strength, "11", m-1, "00", m
        os += strength, "-1-1", m-1, "00", m
        m += 1
        os += 2strength, "E2", m
        os += strength, "00", m-1, "11", m
        os += strength, "00", m-1, "-1-1", m
        os += strength, "11", m-1, "00", m
        os += strength, "-1-1", m-1, "00", m
        m += 1
        os += strength, "00", m-1, "11", m
        os += strength, "00", m-1, "-1-1", m
        os += strength, "11", m-1, "00", m
        os += strength, "-1-1", m-1, "00", m
    elseif type == "cc_breaking"
        m = L ÷ 2
        os += (-1)^(m-1)*strength, "E", m-1, "00", m
        os += -(-1)^(m-1)*strength, "00", m-1, "E", m
        m += 1
        os += (-1)^(m-1)*strength, "E", m-1, "00", m
        os += -(-1)^(m-1)*strength, "00", m-1, "E", m
        m += 1
        os += (-1)^(m-1)*strength, "E", m-1, "00", m
        os += -(-1)^(m-1)*strength, "00", m-1, "E", m
    end
    
    return MPO_new(os, sites) 
end