

function createdriversets()
		
	#Find the effective shifts 
	#i.e. For a three day horizon starting on a Monday, a MTWTF driver looks the same as a MTWR+Su driver
	effectiveshifts = unique!([T_off[s] for s in 1:numshifts])
	numeffshifts = length(effectiveshifts)
	effshift, shiftsincluded = [], [[] for es in 1:numeffshifts]
	for s in 1:numshifts
		es = findfirst(x->x==T_off[s], effectiveshifts)
		push!(effshift, es)
		push!(shiftsincluded[es], s)
	end

	#Create driver sets
	driversets = Dict()
	for l in 1:numlocs, s in 1:numeffshifts
		driversets[l,s] = []
	end
	for d in drivers
		push!(driversets[driverHomeLocs[d], effshift[drivershift[d]]], d)
	end

	driverSetStartNodes = Dict()
	for l in 1:numlocs, s in 1:numeffshifts
		driverSetStartNodes[l,s] = nodes[l,0]
	end

	return driversets, driverSetStartNodes, numeffshifts, effshift, shiftsincluded

end

#----------------------------------------------------------------------------------------#

function addnextsegment(currentfrags, completedfrags, adjacentlocs, deadline, homeloc)
	
	currfragnum = length(currentfrags)

	for f in 1:currfragnum
		frag = popfirst!(currentfrags)
		currloc, currtime = last(frag)
		if mod(currtime, 24) < shiftlength #Working hours
			for a in union(adjacentlocs[currloc], currloc)
				if (a == currloc) & (currtime + tstep <= deadline - traveltimebetweenlocs_rdd[a, homeloc]) 
					newfrag = deepcopy(frag)
					push!(newfrag, (a, currtime + tstep))
					push!(currentfrags, newfrag)
					if currtime + tstep == deadline
						push!(completedfrags, newfrag)
					end
				elseif (a != currloc) && (mod(currtime + arcLength[currloc, a],24) <= shiftlength) && (currtime + arcLength[currloc, a] <= deadline - traveltimebetweenlocs_rdd[a, homeloc]) 
					newfrag = deepcopy(frag)
					push!(newfrag, (a, currtime + arcLength[currloc, a]))
					push!(currentfrags, newfrag)
					if currtime + arcLength[currloc, a] == deadline
						push!(completedfrags, newfrag)
					end
				end
			end
		elseif currloc != homeloc
			newfrag = deepcopy(frag)
			push!(newfrag, (currloc, currtime + tstep))
			push!(currentfrags, newfrag)
			if currtime + tstep == deadline
				push!(completedfrags, newfrag)
			end
		end
	end
	
	currentfrags = setdiff(currentfrags, completedfrags)

	return currentfrags, completedfrags

end

#----------------------------------------------------------------------------------------#

function findbasefragments(homeloc, maxnightsaway, adjacentlocs)
		
	basefragments = []

	for deadline in shiftlength:24:shiftlength+24*maxnightsaway
		currentfrags = [[(homeloc, 0)]] 
		completedfrags = []
		
		#Loop to enumerate paths
		while currentfrags != []
			currentfrags, completedfrags = addnextsegment(currentfrags, completedfrags, adjacentlocs, deadline, homeloc)
	 	end

	 	#Add the paths of length "deadline" to the list of fragments 
	 	basefragments = union(basefragments, completedfrags)
	end

	return basefragments

end

#----------------------------------------------------------------------------------------#

function initializejourneymodel(maxnightsaway)

	driversets, driverSetStartNodes, numeffshifts, effshift, shiftsincluded = createdriversets()

	adjacentlocs = [[] for l in 1:numlocs]
	for pa in prearcs
		if pa[3] <= shiftlength
			push!(adjacentlocs[pa[1]], pa[2])
		end
	end

	numfragments = Dict()
	fragmentscontaining = Dict()
	F_plus_ls, F_minus_ls = Dict(), Dict()
	fragdrivinghours, fragworkinghours = Dict(), Dict()
	for l in 1:numlocs, s in 1:numeffshifts
		numfragments[l,s] = 0
	end
	for l in 1:numlocs, s in 1:numeffshifts, a in 1:numarcs
		fragmentscontaining[l,s,a] = []
	end
	for l in 1:numlocs, s in 1:numeffshifts, n in 1:numnodes
		F_plus_ls[l,s,n] = []
		F_minus_ls[l,s,n] = []
	end

	for l in 1:numlocs
		basefragments = findbasefragments(l, maxnightsaway, adjacentlocs)
		for s in 1:numeffshifts
			#println("------ $l, $s ------")
			fragindex = 1
			if driversets[l,s] != []
				d_ex = driversets[l,s][1]
				#Dummy fragment to stay home during off hours at time 0
				if !(0 in T_on_0[d_ex])
					frag = [(l, t) for t in 0:tstep:T_on_0[d_ex][1]]

					#Create official timed fragment 
					currnode = nodes[frag[1][1], frag[1][2]]
					push!(F_plus_ls[l,s,currnode], fragindex)

					for n in 2:length(frag)
						newnode = nodes[frag[n][1], frag[n][2]]
						newarc = arcs[currnode, newnode]
						push!(fragmentscontaining[l,s,newarc], fragindex)
						currnode = newnode
					end

					push!(F_minus_ls[l,s,currnode], fragindex)

					#println(fragindex, " ==> ", frag, " @ t = 0 (pre-shift)")

					fragdrivinghours[l,s,fragindex] = 0
					fragworkinghours[l,s,fragindex] = 0

					fragindex += 1
					numfragments[l,s] += 1

				end	

				for t in [t2 for t2 in T_on_0[d_ex] if t2 <= horizon] #All possible start times of plans
					for frag in basefragments
						fragmentendtime = t + last(frag)[2]
						
						if (fragmentendtime - shiftlength in T_on_0[d_ex]) || (fragmentendtime - shiftlength in T_on_0[d_ex])
							#Create official timed fragment 
							currnode = nodes[frag[1][1], frag[1][2] + t]
							push!(F_plus_ls[l,s,currnode], fragindex)

							fragdrivinghours[l,s,fragindex] = 0
							fragworkinghours[l,s,fragindex] = 0

							fragarcs = []
							for n in 2:length(frag)
								if frag[n][2] + t <= horizon
									newnode = nodes[frag[n][1], frag[n][2] + t]
									newarc = arcs[currnode, newnode]
									push!(fragmentscontaining[l,s,newarc], fragindex)
									push!(fragarcs, newarc)
									currnode = newnode
								end
							end

							#Add the off hours rest period to the end of each fragment to connect to the start of the next fragment
							for trest in fragmentendtime+tstep:tstep:min(horizon, fragmentendtime+24-shiftlength)
								newnode = nodes[l, trest]
								newarc = arcs[currnode, newnode]
								push!(fragmentscontaining[l,s,newarc], fragindex)
								push!(fragarcs, newarc)
								currnode = newnode
							end

							push!(F_minus_ls[l,s,currnode], fragindex)

							for a in fragarcs
								l1,t1 = nodesLookup[arcLookup[a][1]]
								l2,t2 = nodesLookup[arcLookup[a][2]]
								if !(t1 in T_off[s]) & (l1 != l2)
									fragdrivinghours[l,s,fragindex] += t2-t1
									fragworkinghours[l,s,fragindex] += t2-t1
								elseif !(t1 in T_off[s]) & (l1 == l2) & !(l1 == l)
									fragworkinghours[l,s,fragindex] += t2-t1
								end
							end	

							#println(fragindex, " ==> ", frag, " @ t = $t")

							fragindex += 1
							numfragments[l,s] += 1
						end
					end
				end
			end
		end
	end

	N_flow_ls = Dict()
	for l in 1:numlocs, s in 1:numeffshifts
		N_flow_ls[l,s] = []
		for n in setdiff(1:numnodes, union(driverSetStartNodes[l,s], N_end))
			#if union(F_plus_ls[l,s,n], F_minus_ls[l,s,n]) != []
			push!(N_flow_ls[l,s], n)
			#end
		end
	end

	workingfragments = Dict()
	for l in 1:numlocs, s in 1:numeffshifts
		workingfragments[l,s] = []
		for f in 1:numfragments[l,s] 
			if fragworkinghours[l,s,f] > 1e-4
				push!(workingfragments[l,s], f)
			end
		end
	end

	return driversets, driverSetStartNodes, numfragments, fragmentscontaining, F_plus_ls, F_minus_ls, N_flow_ls, numeffshifts, effshift, shiftsincluded, fragdrivinghours, fragworkinghours, workingfragments

end 
