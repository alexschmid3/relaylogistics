
function updatedriverarcsets()

	if mod(timedelta,24) == 0
		#No need to update because the shift schedule didn't change
		1+1
	else
		1+1
		#throw(DomainError(timedelta, "timedelta not valid: must be multiple of 24"))
	end

end

#---------------------------------------------------------------------------------------#

function updatedriverjourneys(enumeratestandardjourneys_flag)

	#Create driver sets and journeys
    @time driversets, driversingroup, numdrivergroups, drivergroupnum, drivergroupdesc, numeffshifts, effshift, shiftsincluded = finddriversets_online(currstate.T_off, currstate.driverStartNodes, currstate.lasttimehome) 
	@time numfragments, fragmentscontaining, fragmentarcs, F_plus_g, F_minus_g, N_flow_g = initializedriversetjourneys(currstate, driversets, drivergroupnum, driversingroup, drivergroupdesc, currarcs.driverarcs, currarcs.ghostdriverarcs, enumeratestandardjourneys_flag)
	@time fragdrivinghours, fragworkinghours, workingfragments, fragmentnightsaway = getfragmentstats(currstate, driversets, numfragments, fragmentarcs)

	#=
	for item in currfragments.driversets
		remove!(currfragments.driversets, item)
	end
	for item in driversets
		push!(currfragments.driversets, item)
	end

	for ky in keys(currfragments.driversingroup)
		currfragments.driversingroup[ky] = driversingroup[ky]
	end
	for ky in keys(drivergroupnum)
		currfragments.drivergroupnum[ky] = drivergroupnum[ky]
		drivergroupdesc[drivergroupnum[ky]] = ky
	end
	for ky in driversets
		currfragments.fragments[ky] = [j for j in 1:numfragments[ky]]
	end


	newfragments = (#driversets=driversets, driversingroup=driversingroup, numdrivergroups=numdrivergroups, drivergroupnum=drivergroupnum, 
                #drivergroupdesc=drivergroupdesc, numeffshifts=numeffshifts, 
				#effshift=effshift, shiftsincluded=shiftsincluded,
                #numfragments=numfragments, 
				fragmentscontaining=fragmentscontaining, 
				fragmentarcs=fragmentarcs, 
                F_plus_g=F_plus_g, F_minus_g=F_minus_g, N_flow_g=N_flow_g, 
                #effshift=effshift, shiftsincluded=shiftsincluded, 
                #fragdrivinghours=fragdrivinghours, fragworkinghours=fragworkinghours, 
                #workingfragments=workingfragments, fragmentnightsaway=fragmentnightsaway)
				=#

	fragments = Dict()
	for item in driversets
		fragments[item] = [j for j in 1:numfragments[item]]
	end

	newfragments = (driversets=driversets, driversingroup=driversingroup, numdrivergroups=numdrivergroups, drivergroupnum=drivergroupnum, 
                drivergroupdesc=drivergroupdesc, numeffshifts=numeffshifts, effshift=effshift, shiftsincluded=shiftsincluded,
                fragments=fragments, numfragments=numfragments, fragmentscontaining=fragmentscontaining, fragmentarcs=fragmentarcs, 
                F_plus_g=F_plus_g, F_minus_g=F_minus_g, N_flow_g=N_flow_g, 
                #effshift=effshift, shiftsincluded=shiftsincluded, 
                fragdrivinghours=fragdrivinghours, fragworkinghours=fragworkinghours, 
                workingfragments=workingfragments, fragmentnightsaway=fragmentnightsaway)

    return newfragments

end

#---------------------------------------------------------------------------------------#

function updateorderarcs()

	if operations == "relay"
		@time orderarcset_full, orderarcset_space_full, A_plus_i_full, A_minus_i_full = orderarcreduction(currstate.orders, currstate.Origin, currstate.Destination)
	elseif operations == "ptp"
		@time orderarcset_full, orderarcset_space_full, A_plus_i_full, A_minus_i_full = pointotpointarcs(currstate.orders, currstate.Origin, currstate.Destination)
	end

	orderarcs = (A=orderarcset_full, A_space=orderarcset_space_full, A_minus=A_minus_i_full, A_plus=A_plus_i_full, available=[], closelocs=[]);

	if operations == "relay"
		@time magarcs = initializeorderarcsets(k, currstate.orders, originloc, destloc, currstate.Origin, currstate.Destination, currstate.shortesttriptimes)
	elseif operations == "ptp"
		@time magarcs = orderarcs
	end

	return orderarcs, magarcs

end

#---------------------------------------------------------------------------------------#
