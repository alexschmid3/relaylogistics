
function initializecurrentstatearcs(currstate)

    primaryarcs, extendedtimearcs, orderarcs, driverarcs, hasdriverarcs = initializearcsets(A_space, A_plus, A_minus, currstate.orders, currstate.Origin, currstate.Destination, currstate.driverStartNodes, currstate.T_off)
    R_off = findreturnhomearcsets(driverarcs, currstate.T_off_constr)
    magarcs = initializeorderarcsets(k, currstate.orders, originloc, destloc, currstate.Origin, currstate.Destination, currstate.shortesttriptimes)
    driversets, driverSetStartNodes, numfragments, fragmentscontaining, F_plus_ls, F_minus_ls, N_flow_ls, numeffshifts, effshift, shiftsincluded, fragdrivinghours, fragworkinghours, workingfragments, fragmentnightsaway = initializejourneymodel(maxnightsaway, currstate.T_off, currstate.T_on_0)

    currarcs = (orderarcs=orderarcs, driverarcs=driverarcs, hasdriverarcs=hasdriverarcs, magarcs=magarcs, R_off=R_off)
    currfragments = (driversets=driversets, driverSetStartNodes=driverSetStartNodes, numfragments=numfragments, 
                fragmentscontaining=fragmentscontaining, F_plus_ls=F_plus_ls, F_minus_ls=F_minus_ls, N_flow_ls=N_flow_ls, 
                effshift=effshift, shiftsincluded=shiftsincluded, fragdrivinghours=fragdrivinghours, 
                fragworkinghours=fragworkinghours, workingfragments=workingfragments, fragmentnightsaway=fragmentnightsaway)

    return currarcs, currfragments, primaryarcs, extendedtimearcs, numeffshifts

end
