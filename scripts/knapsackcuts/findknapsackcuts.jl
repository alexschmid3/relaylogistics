
using Combinatorics

function findjourneypaircuts(z)

    cutindex = 0
    cutvars, cutrhs = Dict(), Dict()
    for d in drivers
        l, s = driverHomeLocs[d], drivershift[d]
        usedjourneys = [f for f in workingfragments[l,s] if value(z[d,f]) > 0]
        journeypairs = [(i,j) for i in usedjourneys, j in usedjourneys if !(i==j) & (value(z[d,i]) + value(z[d,j]) >= 1.0001) & (fragworkinghours[l,s,i]+fragworkinghours[l,s,j] > maxweeklydriverhours)]
        for (i,j) in journeypairs
            cutindex += 1
            cutvars[cutindex] = [(d,i), (d,j)]
            cutrhs[cutindex] = 1
        end
    end

    return cutindex, cutvars, cutrhs 

end

#---------------------------------------------------------------------------------------#

function findjourneytripletcuts(z)

    cutindex = 0
    cutvars, cutrhs = Dict(), Dict()
    for d in drivers
        l, s = driverHomeLocs[d], drivershift[d]
        usedjourneys = [f for f in workingfragments[l,s] if (value(z[d,f]) > 0) ] #& (fragworkinghours[l,s,f] < 24)]
        journeypairs = [(i,j,k) for i in usedjourneys, j in usedjourneys, k in usedjourneys if (i<j) & (j<k)  & (value(z[d,i]) + value(z[d,j]) + value(z[d,k]) >= 2.0001) & (fragworkinghours[l,s,i]+fragworkinghours[l,s,j]+fragworkinghours[l,s,k] > maxweeklydriverhours)]
        for (i,j,k) in journeypairs
            cutindex += 1
            cutvars[cutindex] = [(d,i), (d,j), (d,k)]
            cutrhs[cutindex] = 2
        end
    end

    return cutindex, cutvars, cutrhs 

end

#---------------------------------------------------------------------------------------#

function checkminimalcover(l,s,journeysubset)

    if sum(fragworkinghours[l,s,j] for j in journeysubset) <= maxweeklydriverhours
        return false
    end

    covercount = 0
    for f in journeysubset
        if sum(fragworkinghours[l,s,j] for j in setdiff(journeysubset,f)) <= maxweeklydriverhours
            covercount += 1
        end 
    end

    if covercount == length(journeysubset)
        return true
    else
        return false
    end

end

#---------------------------------------------------------------------------------------#

function findliftcoeffs_opt(ss, j, rhs, l, s)

    model = Model(Gurobi.Optimizer) #Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_optimizer_attribute(model, "OutputFlag", 0)
    @variable(model, z[1:numfragments[l,s]], Bin)
    @objective(model, Max, sum(z[f] for f in ss))
    @constraint(model, sum(fragworkinghours[l,s,f] * z[f] for f in 1:numfragments[l,s]) <= maxweeklydriverhours)
    @constraint(model, z[j] == 1)
    optimize!(model)

    return rhs - objective_value(model)

end

#---------------------------------------------------------------------------------------#

#Iterate over all possible lifting coefficients for journey j in the min cover cut described by variables in the set ss and a right-hand side of rhs
#This is equivalent to the optimization above, but is faster and doesn't call Gurobi as much
function findliftcoeffs(ss, j, rhs, l, s)

    sortedss = sort(ss, by=x->fragdrivinghours[l,s,x])
    for n in rhs:-1:0 #Max coeff is rhs, min coeff is 0

        #Check whether we could choose j plus the n journeys with the fewest working hours and still satisfy the cut   
        if (n == 0) & (fragdrivinghours[l,s,j] <= maxweeklydriverhours)
            return rhs
        elseif (n>0) & (sum(fragdrivinghours[l,s,f] for f in sortedss[1:n]) + fragdrivinghours[l,s,j] <= maxweeklydriverhours)
            return rhs - n
        end

    end

    #If journey j does not satify the knapsack, it's lifting coeff should be ∞ 
    #Note: shouldn't happen bc of pre-processing
    return 1e10

end

#---------------------------------------------------------------------------------------#

function sequentiallifting(ss, coeff, usedjourneys, journeysubset, bestrhs, l, s)
    
    #Iterate over all journeys and try to lift each coefficient
    for j2 in union(setdiff(usedjourneys, journeysubset), setdiff(workingfragments[l,s], usedjourneys))
        liftcoeff = findliftcoeffs(ss, j2, bestrhs, l, s)
        if liftcoeff > 1e-4
            push!(ss, j2)
            coeff[j2] = liftcoeff
        end
    end

    return ss, coeff

end

#---------------------------------------------------------------------------------------#

#Parameters described in section 4 of Gu, Nemhauser, and Savelsbergh (2000)
function calcparametersforseqindeplift(cover, covermap, l, s)

    r = length(cover)
    λ = sum(fragworkinghours[l,s,j] for j in cover) - maxweeklydriverhours

    μ, ρ = Dict(), Dict()
    μ[0] = 0
    for i in 1:r
        μ[i] = sum(fragworkinghours[l,s,covermap[j]] for j in 1:i)
    end
    for i in 0:r-1
        ρ[i] = max(0, fragworkinghours[l,s,covermap[i+1]] - (fragworkinghours[l,s,covermap[1]] - λ))
    end

    return μ, λ, ρ, r

end

#---------------------------------------------------------------------------------------#

#Function described in section 4 of Gu, Nemhauser, and Savelsbergh (2000)
function g(z, r, μ, λ, ρ)

    if z == 0
        return 0
    else
        for h in 0:r-1
            if μ[h] - λ + ρ[h] < z <= μ[h+1] - λ
                return h
            end
        end
        for h in 1:r-1
            if μ[h] - λ < z <= μ[h] - λ + ρ[h]
                return h - (μ[h] - λ + ρ[h] - z) / ρ[1]
            end
        end
    end

    #Catch all for any journeys that are longer than the max hours
    if z >= (maxnightsaway + 1) * shiftlength
        return 0
    end

end

#---------------------------------------------------------------------------------------#

function sequenceindependentlifting(ss, coeff, journeysubset, l, s)

    #Sort and format the cover
    cover = sort(journeysubset, by=x->fragworkinghours[l,s,x], rev=true)
    jindex, covermap = 1, Dict()
    for j in cover
        covermap[jindex] = j
        jindex += 1
    end
    
    #Calculate lifting coefficients
    μ, λ, ρ, r = calcparametersforseqindeplift(cover, covermap, l, s)
    liftcoeffforhours = Dict()
    for z in tstep:tstep:shiftlength*(maxnightsaway+1)
        liftcoeffforhours[z] = g(z, r, μ, λ, ρ)
    end

    #Lift all coefficients at the same time
    for j2 in setdiff(workingfragments[l,s], journeysubset)
        if roundeddrivinghours_flag == 1
            if liftcoeffforhours[fragworkinghours[l,s,j2]] > 1e-4
                push!(ss, j2)
                coeff[j2] = liftcoeffforhours[fragworkinghours[l,s,j2]]
            end
        elseif roundeddrivinghours_flag == 0
            mycoeff = g(fragworkinghours[l,s,j2], r, μ, λ, ρ)
            if mycoeff > 1e-4
                push!(ss, j2)
                coeff[j2] = mycoeff
            end
        end
    end

    return ss, coeff

end

#---------------------------------------------------------------------------------------#

function findminimalcovercuts(z, lifting_flag)

    cutindex = 0
    cutvars, cutrhs = Dict(), Dict()
    cutcoeff = Dict()

    #Search for cuts over each driver
    for d in drivers
        l, s = driverHomeLocs[d], drivershift[d]

        #List the journeys used in current solutions
        usedjourneys = [f for f in workingfragments[l,s] if value(z[d,f]) > 0]

        #Enumerate over all subsets of used journeys to identify minimal covers
        for journeysubset in powerset(usedjourneys)
            if length(journeysubset) >= 2 
                
                #If minimal cover found, check whether the min cover cut is violated by current solution
                if checkminimalcover(l,s,journeysubset)

                    bestrhs = length(journeysubset) - 1

                    #If min cover cut is violated, format and lift the cut
                    if sum(value(z[d,f]) for f in journeysubset) > bestrhs + 1e-4

                        #Save min cover cut info (all coefficients are 1)
                        ss = copy(journeysubset)
                        coeff = Dict()
                        for j in ss
                            coeff[j] = 1
                        end

                        #Lifting
                        if lifting_flag == 1
                            ss, coeff = sequentiallifting(ss, coeff, usedjourneys, journeysubset, bestrhs, l, s)
                        elseif lifting_flag == 2
                            ss, coeff = sequenceindependentlifting(ss, coeff, journeysubset, l, s)
                        end

                        #Add the cut for all drivers identical to driver d to avoid adding the same cut for another driver in a future iteration
                        for d_sym in driversets[l,s]
                           
                            #Format and save the cut information (variables, coeffs, rhs)
                            cutindex += 1
                            cutvars[cutindex] = [(d_sym,i) for i in ss]
                            cutrhs[cutindex] = bestrhs
                            cutcoeff[cutindex] = Dict()
                            for i in ss
                                cutcoeff[cutindex][d_sym,i] = coeff[i]
                            end  
                        end
                    end
                end
            end
        end
    end

    return cutindex, cutvars, cutrhs, cutcoeff

end

#---------------------------------------------------------------------------------------#

function findstrongestcutrhs(l, s, journeysubset)

    maxhrs = maximum([fragworkinghours[l,s,j] for j in journeysubset])
    target_coeff = floor(maxhrs/2)+1
    adj_rhs = maxweeklydriverhours + sum(max(0,target_coeff - fragworkinghours[l,s,j]) for j in journeysubset)
    rdd_rhs = floor(adj_rhs/target_coeff)

    return rdd_rhs

end

#---------------------------------------------------------------------------------------#

function findallcuts(z)

    cutindex = 0
    cutvars, cutrhs = Dict(), Dict()
    for d in drivers
        #println("Searching driver $d...")
        l, s = driverHomeLocs[d], drivershift[d]
        usedjourneys = [f for f in workingfragments[l,s] if value(z[d,f]) > 0]
        for journeysubset in powerset(usedjourneys)
            if length(journeysubset) >= 2
                bestrhs = findstrongestcutrhs(l, s, journeysubset)
                
                if sum(value(z[d,f]) for f in journeysubset) > bestrhs
                    for d_sym in driversets[l,s]
                        cutindex += 1
                        cutvars[cutindex] = [(d_sym,i) for i in journeysubset]
                        cutrhs[cutindex] = bestrhs
                    end
                end
            end
        end
    end

    return cutindex, cutvars, cutrhs 

end

#---------------------------------------------------------------------------------------#

function findknapsackcuts(z, cuttype)

    if cuttype == 0
        numcuts, cutvars, cutrhs, cutcoeff = 0, (), (), ()
    elseif cuttype == 1        #Find pairs of journeys that can be rounded for a cut
        numcuts, cutvars, cutrhs = findjourneypaircuts(z)
    elseif cuttype == 2    #Find journey triplets that can be rounded for a cut
        numcuts, cutvars, cutrhs = findjourneytripletcuts(z)
    elseif cuttype == 3    #Find all subsets of journeys that can be rounded for a cut
        numcuts, cutvars, cutrhs = findallcuts(z)
    elseif cuttype == 4    #Minimal cover cuts, no lifting
        numcuts, cutvars, cutrhs, cutcoeff = findminimalcovercuts(z, 0)
    elseif cuttype == 5    #Minimal cover cuts + sequential lifting
        numcuts, cutvars, cutrhs, cutcoeff = findminimalcovercuts(z, 1)
    elseif cuttype == 6    #Minimal cover cuts + sequence independent lifting
        numcuts, cutvars, cutrhs, cutcoeff = findminimalcovercuts(z, 2)
    end

    cuts = (numcuts=numcuts, vars=cutvars, rhs=cutrhs, coeff=cutcoeff)

    return cuts

end
