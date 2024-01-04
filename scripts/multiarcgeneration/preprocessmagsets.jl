

function preprocesssubproblemsets(orderArcSet)

	#Create copies of all important sets (so they can be modified for the shortest path problem)
	nodesLookupSP = deepcopy(nodesLookup)
	arcsSP, arcLookupSP, nodesLookupSP = Dict(), Dict(), Dict()
	for i in orders
		arcsSP[i], arcLookupSP[i] = deepcopy(arcs), deepcopy(arcLookup)
	end

	#Create arclists
	arclistSP = Dict()
	for i in orders
		arclistSP[i] = deepcopy(setdiff(orderArcSet[i], dummyarc))
	end

	#Create dummy origin and destination nodes + arcs 
	dummyorig = extendednumnodes + 1
	dummydest = extendednumnodes + 2
	numnodesSP = extendednumnodes + 2
	nodesLookupSP[dummyorig], nodesLookupSP[dummydest] = (-1,-1), (-1,99999)

	#Add dummy arcs
	for i in orders
		index = extendednumarcs + 2

		#Add arcs from dummy origin to all origin nodes in pickup window
		for j in 1:length(Origin[i])
			n = Origin[i][j]
			arcsSP[i][dummyorig, n] = index
			arcLookupSP[i][index] = (dummyorig, n)
			pushfirst!(arclistSP[i], index)
			index += 1
		end

		#Add arcs from all destination nodes in delivery window to the dummy destination
		for n in Destination[i]
			arcsSP[i][n, dummydest] = index
			arcLookupSP[i][index] = (n, dummydest) 
			push!(arclistSP[i], index)
			index += 1
		end
	end

	#Remove arcs that are not feasible because they are longer than the driver shift
	for i in orders, a in 1:numarcs
		arctime = nodesLookup[arcLookup[a][2]][2] - nodesLookup[arcLookup[a][1]][2]
		if arctime > shiftlength
			remove!(arclistSP[i], a)
		end
	end

	return numnodesSP, arclistSP, nodesLookupSP, arcLookupSP, dummyorig, dummydest

end

#----------------------------------------------------------------------------------------#

function alphatransformation()

    M_alpha = Dict()
    for i in orders
        M_alpha[i] = spzeros(extendednumarcs, extendednumnodes)
        for n in setdiff([n2 for n2 in 1:numnodes], union(Origin[i], Destination[i])), a in A_minus[n]
            M_alpha[i][a,n] -= 1
        end
        for n in setdiff([n2 for n2 in 1:numnodes], union(Origin[i], Destination[i])), a in A_plus[n]
            M_alpha[i][a,n] += 1
        end
    end

    M_beta = Dict()
    for i in orders
        M_beta[i] = spzeros(extendednumarcs)
        for n in Origin[i], a in intersect(union(A_space, dummyarc), A_plus[n])
            M_beta[i][a] -= 1
        end
        if !(i in ordersinprogress)
            a = extendedarcs[last(Origin[i]), last(Destination[i])]
            M_beta[i][a] -= 1
        else
            for n in Origin[i], a in setdiff(A_plus[n], union(A_space, dummyarc))
                M_beta[i][a] -= 1
            end
        end
    end

    M_gamma = Dict()
    for i in orders
        M_gamma[i] = spzeros(extendednumarcs)
        for n in Destination[i], a in A_minus[n]
            if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))
                M_gamma[i][a] -= 1
            end
        end
	end

    M_theta = spzeros(extendednumarcs, length(N_0))
    for n in N_0, a in setdiff(A_plus[n], dummyarc)  
        M_theta[a,n] -= 1
    end

    M_nu = spzeros(extendednumarcs, length(N_end))
    adjuster = N_end[1] - 1
    for n in N_end, a in A_minus[n]
		M_nu[a,n-adjuster] -= 1
	end

    M_mu = spzeros(extendednumarcs, numnodes - length(N_0) - length(N_end))
    adjuster = last(N_0)
    for n in N_flow_t, a in A_minus[n]
		M_mu[a,n-adjuster] -= 1
	end
	for n in N_flow_t, a in A_plus[n]
		M_mu[a,n-adjuster] += 1
	end

    M_psi = Dict() 
    for i in orders
        M_psi[i] = spzeros(extendednumarcs)
        for n in Destination[i], a in A_minus[n]
		    M_psi[i][a] += arcfinishtime[a]
        end
	end

    M_xi = spzeros(extendednumarcs, length(A_space))
    for a in A_space
		M_xi[a,a] += 1
	end

    M = (alpha=M_alpha, beta=M_beta, gamma=M_gamma, theta=M_theta, nu=M_nu, mu=M_mu, psi=M_psi, xi=M_xi);

    return M

end

#----------------------------------------------------------------------------------------#

function preprocessmagsets(orderArcSet)

    M = alphatransformation()
    numnodesSP, arclistSP, nodesLookupSP, arcLookupSP, dummyorig, dummydest = preprocesssubproblemsets(orderArcSet)
    subproblemsets = (numnodes=numnodesSP, arclist=arclistSP, nodelookup=nodesLookupSP, arclookup=arcLookupSP, dummyorig=dummyorig, dummydest=dummydest)

    return M, subproblemsets

end
