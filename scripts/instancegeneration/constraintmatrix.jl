
#----------------------------------------------------------------------------------------#

function generatecoeffs_orderflowbalance(A_minus_i, A_plus_i)

	xcoeffs_alpha = Dict()

	for i in orders
		alpha_row, alpha_col, alpha_val = [], [], []
		for n in setdiff([n2 for n2 in 1:numnodes], union(Origin[i], Destination[i]))
			for a in A_minus_i[i,n]
				push!(alpha_row, a)
				push!(alpha_col, n)
				push!(alpha_val, 1)
			end
			for a in A_plus_i[i,n]
				push!(alpha_row, a)
				push!(alpha_col, n)
				push!(alpha_val, -1)
			end
		end
		xcoeffs_alpha[i] = sparse(alpha_row, alpha_col, alpha_val)
	end

	return xcoeffs_alpha

end

#----------------------------------------------------------------------------------------#

function generatecoeffs_departorigin(A_plus_i)

	beta_row, beta_col, beta_val = [], [], []
	for i in orders
		for n in Origin[i], a in intersect(union(A_space, dummyarc), A_plus_i[i,n])
			push!(beta_row, a)
			push!(beta_col, i)
			push!(beta_val, 1)
		end
		if !(i in ordersinprogress)
			extendedorderarc = extendedarcs[last(Origin[i]), last(Destination[i])]
			if extendedorderarc in orderArcSet[i]
				push!(beta_row, extendedorderarc)
				push!(beta_col, i)
				push!(beta_val, 1)
			end
		else  #online implementation - allow order to stay at "origin" if we've already departed from pickup location
			for n in Origin[i], a in setdiff(A_plus[n], union(A_space, dummyarc))
				if a in orderArcSet[i]
					push!(beta_row, a)
					push!(beta_col, i)
					push!(beta_val, 1)
				end
			end
		end
	end

	xcoeffs_beta = sparse(beta_row, beta_col, beta_val)

	return xcoeffs_beta

end

#----------------------------------------------------------------------------------------#

function generatecoeffs_arrivedestin(A_minus_i)

	xcoeffs_gamma = Dict()

	gamma_row, gamma_col, gamma_val = [], [], []
	for i in orders
		for n in Destination[i], a in A_minus_i[i,n] 
			if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))
				push!(gamma_row, a)
				push!(gamma_col, i)
				push!(gamma_val, 1)
			end
		end
	end

	xcoeffs_gamma = sparse(gamma_row, gamma_col, gamma_val)

	return xcoeffs_gamma

end

#----------------------------------------------------------------------------------------#

function generatecoeffs_initialtrucks(A_plus_i)

	xcoeffs_theta = Dict()

	for i in orders
		theta_row, theta_col, theta_val = [], [], []
		for n in N_0
			for a in setdiff(A_plus_i[i,n], dummyarc)
				push!(theta_row, a)
				push!(theta_col, n)
				push!(theta_val, 1)
			end
		end
		xcoeffs_theta[i] = sparse(theta_row, theta_col, theta_val)
	end

	return xcoeffs_theta

end

#----------------------------------------------------------------------------------------#

function generatecoeffs_finaltrucks(A_minus_i)

	xcoeffs_nu = Dict()

	for i in orders
		nu_row, nu_col, nu_val = [], [], []
		for n in N_end
			for a in setdiff(A_minus_i[i,n], dummyarc)
				push!(nu_row, a)
				push!(nu_col, n)
				push!(nu_val, 1)
			end
		end
		xcoeffs_nu[i] = sparse(nu_row, nu_col, nu_val)
	end

	return xcoeffs_nu

end

#----------------------------------------------------------------------------------------#

function generatecoeffs_truckflowbalance(A_minus_i, A_plus_i)

	xcoeffs_mu = Dict()

	for i in orders
		mu_row, mu_col, mu_val = [], [], []
		for n in N_flow_t
			for a in setdiff(A_minus_i[i,n], dummyarc)
				push!(mu_row, a)
				push!(mu_col, n)
				push!(mu_val, 1)
			end
			for a in setdiff(A_plus_i[i,n], dummyarc)
				push!(mu_row, a)
				push!(mu_col, n)
				push!(mu_val, -1)
			end
		end
		xcoeffs_mu[i] = sparse(mu_row, mu_col, mu_val)
	end

	return xcoeffs_mu

end

#----------------------------------------------------------------------------------------#

function generatecoeffs_driveravailability(A_minus_i)

	xcoeffs_xi = Dict()

	xi_row, xi_col, xi_val = [], [], []
	for i in orders
		for a in orderArcSet_space[i]
			push!(xi_row, a)
			push!(xi_col, i)
			push!(xi_val, -1)
		end
	end

	xcoeffs_xi = sparse(xi_row, xi_col, xi_val)

	return xcoeffs_xi

end

#----------------------------------------------------------------------------------------#

function generatecoeffs_deliverytime(A_minus_i)

	xcoeffs_psi = Dict()

	psi_row, psi_col, psi_val = [], [], []
	for i in orders
		for n in Destination[i], a in A_minus_i[i,n]
			push!(psi_row, a)
			push!(psi_col, i)
			push!(psi_val, -1 * arcfinishtime[a])
		end
	end

	xcoeffs_psi = sparse(psi_row, psi_col, psi_val)

	return xcoeffs_psi

end

#----------------------------------------------------------------------------------------#

function findxcoefficientmatrix(A_minus_i, A_plus_i)

	xcoeffs_alpha = generatecoeffs_orderflowbalance(A_minus_i, A_plus_i)
	xcoeffs_beta = generatecoeffs_departorigin(A_plus_i)
	xcoeffs_gamma = generatecoeffs_arrivedestin(A_minus_i)
	xcoeffs_theta = generatecoeffs_initialtrucks(A_plus_i)
	xcoeffs_nu = generatecoeffs_finaltrucks(A_minus_i)
	xcoeffs_mu = generatecoeffs_truckflowbalance(A_minus_i, A_plus_i)
	xcoeffs_xi = generatecoeffs_driveravailability(A_minus_i)
	xcoeffs_psi = generatecoeffs_deliverytime(A_minus_i)

	return xcoeffs_alpha, xcoeffs_beta, xcoeffs_gamma, xcoeffs_theta, xcoeffs_nu, xcoeffs_mu, xcoeffs_xi, xcoeffs_psi

end

