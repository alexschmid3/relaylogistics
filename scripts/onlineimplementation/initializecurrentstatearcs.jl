
function finddriversets_online(T_off_online, driverStartNodes, awaylastnight)
		
	#Find the effective shifts 
	#i.e. For a three day horizon starting on a Monday, a MTWTF driver looks the same as a MTWR+Su driver
	effectiveshifts = unique!(T_off_online)
	numeffshifts = length(effectiveshifts)
	effshift, shiftsincluded = [], [[] for es in 1:numeffshifts]
	for s in 1:numshifts
		es = findfirst(x->x==T_off_online[s], effectiveshifts)
		push!(effshift, es)
		push!(shiftsincluded[es], s)
	end

	#Create driver sets
	driversingroup = Dict()
	driversets = []
	for homeloc in 1:numlocs, shiftsched in 1:numshifts, driverawaylastnight in 0:1, startnode in 1:numnodes
		driversingroup[homeloc, shiftsched, startnode, driverawaylastnight] = []
	end
	for d in drivers
		push!(driversingroup[driverHomeLocs[d], drivershift[d], driverStartNodes[d], awaylastnight[d]], d)
	end
	for item in driversingroup
		if item[2] != []
			push!(driversets, item[1])
		end
	end

	#Group lookups
	drivergroupnum, drivergroupdesc = Dict(), Dict()
	index = 1
	for (hl,ss,sn,aln) in driversets
		drivergroupnum[hl,ss,sn,aln] = index
		drivergroupdesc[index] = (hl,ss,sn,aln)
		index += 1
	end
	numdrivergroups = length(drivergroupdesc)

	return driversets, driversingroup, numdrivergroups, drivergroupnum, drivergroupdesc, numeffshifts, effshift, shiftsincluded

end

#----------------------------------------------------------------------------------------#

function addnextsegment_online(g, ss, currentfrags, completedfrags, deadline, offtimes, homeloc, driversingroup, drivergroupdesc, driverarcs)
	
	currfragnum = length(currentfrags)

	println(keys(driverarcs))

	#Note this does not work for irregular work schedules!!
	offtimes_ext = union(offtimes, [t+24 for t in intersect([t for t in horizon-3*tstep:tstep:horizon],offtimes)])
	T_on_0_ext = [t for t in setdiff(0:tstep:horizon+24, offtimes_ext) if !(t-tstep in setdiff(0:tstep:horizon+24, offtimes_ext))]
	offtimes_0 = [t for t in offtimes_ext if !(t-tstep in offtimes_ext)]	

	#Add next segment to each fragment
	for f in 1:currfragnum
		frag = popfirst!(currentfrags)
		currloc, currtime = last(frag)
		d_ex = driversingroup[drivergroupdesc[g]][1]
		endofcurrshift = minimum(union([t for t in offtimes_0 if t > currtime], horizon+24))
		startofnextshift = minimum(union([t for t in T_on_0_ext if t > currtime], horizon+24))
		println("   $currloc, $currtime")
		for a in driverarcs.A_plus[d_ex, nodes[currloc, currtime]]
			arcendloc, arcendtime = nodesLookup[arcLookup[a][2]]
			if traveltimebetweenlocs_rdd[arcendloc, homeloc] > 0
				effectivedeadline = maximum([t2 for t2 in offtimes_0 if t2 <= deadline])
			else
				effectivedeadline = deadline
			end
			if (arcendtime <= min(effectivedeadline - traveltimebetweenlocs_rdd[arcendloc, homeloc], endofcurrshift)) & (!((arcendloc == homeloc) & (currtime in offtimes_ext) & (setdiff([t for t in startofnextshift:tstep:deadline-tstep], offtimes) != [])) || startofnextshift > horizon)
				newfrag = deepcopy(frag)
				push!(newfrag, (arcendloc, arcendtime))
				if (arcendtime == deadline) || (arcendtime == horizon)
					push!(completedfrags, newfrag)
				else 
					push!(currentfrags, newfrag)
				end
			end
		end
	end

	return currentfrags, completedfrags

end

#----------------------------------------------------------------------------------------#

function createfragmentsets_online(currstate, hl, ss, sn, aln, drivergroupnum, driversingroup, drivergroupdesc, driverarcs)

	journeys = []

	#Get driver group information
	startloc, starttime = nodesLookup[sn]
	g = drivergroupnum[hl,ss,sn,aln]
	offtimes = currstate.T_off[ss]
	upcomingfragmentstarttime = [t for t in currstate.T_on_0[ss] if t > starttime]

	#Add fragments starting from the drivers' start node
	timeuntilnextshift = upcomingfragmentstarttime[1]
	if 0 in offtimes
		effective_aln = 0
	else
		effective_aln = aln
	end
	for deadline in timeuntilnextshift:24:timeuntilnextshift+24*(maxnightsaway - effective_aln)
		currentfrags = [[(startloc, starttime)]] 
		completedfrags = []
		while currentfrags != []
	 		currentfrags, completedfrags = addnextsegment_online(g, ss, currentfrags, completedfrags, deadline, offtimes, hl, driversingroup, drivergroupdesc, driverarcs)
	 	end
	 	journeys = union(journeys, completedfrags)
	end

	#Add fragments starting from future times
	for t in upcomingfragmentstarttime, deadline in 24:24:24+24*maxnightsaway
		currentfrags = [[(hl, t)]] 
		completedfrags = []
		while currentfrags != []
	 		currentfrags, completedfrags = addnextsegment_online(g, ss, currentfrags, completedfrags, t+deadline, offtimes, hl, driversingroup, drivergroupdesc, driverarcs)
	 	end
	 	journeys = union(journeys, completedfrags)
	end

	return journeys

end 


#-----------------------------------------------------------------------------------------#

function initializedriversetjourneys(currstate, driversets, drivergroupnum, driversingroup, drivergroupdesc, driverarcs)
		
	numfragments = Dict()
	fragmentscontaining = Dict()
	fragmentarcs = Dict()
	F_plus_g, F_minus_g = Dict(), Dict()
	N_flow_g = Dict()

    for (hl,ss,sn,aln) in driversets
        numfragments[hl,ss,sn,aln] = 0
        for a in 1:numarcs
            fragmentscontaining[hl,ss,sn,aln,a] = []
        end
        for n in 1:numnodes
            F_plus_g[hl,ss,sn,aln,n] = []
            F_minus_g[hl,ss,sn,aln,n] = []
        end
    
        N_flow_g[hl,ss,sn,aln] = []
        for n in setdiff(1:numnodes, union(sn, N_end))
            push!(N_flow_g[hl,ss,sn,aln], n)
        end

        journeys = createfragmentsets_online(currstate, hl,ss,sn,aln, drivergroupnum, driversingroup, drivergroupdesc, driverarcs)

        #Process the list of journeys
        for journey in journeys
            
            #Get journey number
            numfragments[hl,ss,sn,aln] += 1
            j = numfragments[hl,ss,sn,aln]
            fragmentarcs[hl,ss,sn,aln,j] = []

            #Find arcs contained on journey
            currnode = journey[1]
            for n in journey[2:length(journey)]
                a = arcs[nodes[currnode],nodes[n]]
                push!(fragmentscontaining[hl,ss,sn,aln,a], j)
                push!(fragmentarcs[hl,ss,sn,aln,j], a)
                currnode = n
            end
            #Get journey flow balance
            push!(F_plus_g[hl,ss,sn,aln, nodes[first(journey)]], j)
            push!(F_minus_g[hl,ss,sn,aln, nodes[last(journey)]], j)
        end

	end

    return numfragments, fragmentscontaining, fragmentarcs, F_plus_g, F_minus_g, N_flow_g

end

#-----------------------------------------------------------------------------------------#

function getfragmentstats(currstate, driversets, numfragments, fragmentarcs)

    fragdrivinghours, fragworkinghours, workingfragments, fragmentnightsaway = Dict(), Dict(), Dict(), Dict()

    for (hl,ss,sn,aln) in driversets, f in numfragments[hl,ss,sn,aln] 
        
        drivinghours, workinghours = 0, 0
        for a in fragmentarcs[hl,ss,sn,aln,f]
            l1,t1 = nodesLookup[arcLookup[a][1]]
            l2,t2 = nodesLookup[arcLookup[a][2]]
            if !(t1 in currstate.T_off[s]) & (l1 != l2)
                drivinghours += t2-t1
                workinghours += t2-t1
            elseif !(t1 in currstate.T_off[s]) & (l1 == l2) & !(l1 == hl)
                workinghours += t2-t1
            end
			#Check me!!
			if (t1 in currstate.T_off_0[s]) & (l1 != hl)
				fragmentnightsaway[hl,ss,sn,aln,f] += 1
			end
        end

        fragdrivinghours[hl,ss,sn,aln,f] = drivinghours
        fragworkinghours[hl,ss,sn,aln,f] = workinghours

        if workinghours > 1e-4
            push!(workingfragments[hl,ss,sn,aln], f)
        end
	end

    return fragdrivinghours, fragworkinghours, workingfragments, fragmentnightsaway

end

#-----------------------------------------------------------------------------------------#

function initializecurrentstatearcs(currstate)

    primaryarcs, extendedtimearcs, orderarcs, driverarcs, hasdriverarcs = initializearcsets(A_space, A_plus, A_minus, currstate.orders, currstate.Origin, currstate.Destination, currstate.driverStartNodes, currstate.T_off)
	R_off = findreturnhomearcsets(driverarcs, currstate.T_off_constr)
    magarcs = initializeorderarcsets(k, currstate.orders, originloc, destloc, currstate.Origin, currstate.Destination, currstate.shortesttriptimes)
    #driversets, driverSetStartNodes, numfragments, fragmentscontaining, F_plus_ls, F_minus_ls, N_flow_ls, numeffshifts, effshift, shiftsincluded, fragdrivinghours, fragworkinghours, workingfragments, fragmentnightsaway = initializejourneymodel(maxnightsaway, currstate.T_off, currstate.T_on_0)
    
    #Create driver sets and journeys
    driversets, driversingroup, numdrivergroups, drivergroupnum, drivergroupdesc, numeffshifts, effshift, shiftsincluded = finddriversets_online(currstate.T_off, currstate.driverStartNodes, currstate.awaylastnight)
    numfragments, fragmentscontaining, fragmentarcs, F_plus_g, F_minus_g, N_flow_g = initializedriversetjourneys(currstate, driversets, drivergroupnum, driversingroup, drivergroupdesc, driverarcs)
	fragdrivinghours, fragworkinghours, workingfragments, fragmentnightsaway = getfragmentstats(currstate, driversets, numfragments, fragmentarcs)

    currarcs = (orderarcs=orderarcs, driverarcs=driverarcs, hasdriverarcs=hasdriverarcs, magarcs=magarcs, R_off=R_off)
    #currfragments = (driversets=driversets, driverSetStartNodes=driverSetStartNodes, numfragments=numfragments, 
    #            fragmentscontaining=fragmentscontaining, F_plus_ls=F_plus_ls, F_minus_ls=F_minus_ls, N_flow_ls=N_flow_ls, 
    #            effshift=effshift, shiftsincluded=shiftsincluded, fragdrivinghours=fragdrivinghours, 
    #            fragworkinghours=fragworkinghours, workingfragments=workingfragments, fragmentnightsaway=fragmentnightsaway)
    currfragments = (driversets=driversets, driversingroup=driversingroup, numdrivergroups=numdrivergroups, drivergroupnum=drivergroupnum, 
                drivergroupdesc=drivergroupdesc, numeffshifts=numeffshifts, effshift=effshift, shiftsincluded=shiftsincluded,
                numfragments=numfragments, fragmentscontaining=fragmentscontaining, fragmentarcs=fragmentarcs, 
                F_plus_g=F_plus_g, F_minus_g=F_minus_g, N_flow_g=N_flow_g, 
                #effshift=effshift, shiftsincluded=shiftsincluded, 
                fragdrivinghours=fragdrivinghours, fragworkinghours=fragworkinghours, 
                workingfragments=workingfragments, fragmentnightsaway=fragmentnightsaway)

    return currarcs, currfragments, primaryarcs, extendedtimearcs, numeffshifts

end
