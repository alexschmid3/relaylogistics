
function writeresults_onlineiteration(filename, currtime)

    df = DataFrame(
        experiment_id = [experiment_id],
        instance = [ex],
        weekstart = [weekstart],
        horizon = [horizon],
        tstep = [tstep],
        lambda = [lambda],
        numdrivers = [length(drivers)],
        operations = [operations],
        onlineweeks = [onlineweeks],
        maxnightsaway = [maxnightsaway],
        maxrepositioningdistance = [maxrepositioningdistance],
        method = [solutionmethod],
        currtime = [currtime],
        currorders = [length(currstate.orders)],
        totalorders = [highestorderindex],
        orderscompleted = [highestorderindex - length(currstate.orders)],
        totalcost = [totalpastcost[currtime]],
        ordertrips = [totalordertrips[currtime]],
        emptytrips = [totalemptytrips[currtime]],
        taxitrips = [totaltaxitrips[currtime]],
        ordermiles = [totalordermiles[currtime]],
        penaltyordermiles = [0],
        emptymiles = [totalemptymiles[currtime]],
        taximiles = [totaltaximiles[currtime]],
        shortesttripmiles = [shortestpossible_ordermiles[currtime]],
        totaldeliverytime = [total_delivtime[currtime]],
        penaltydeliverytime = [0],
        maxdeliverytime = [max_delivtime[currtime]],
        shortestdeliverytime_abs = [0],
        shortestdeliverytime_rest = [shortestpossible_delivtime[currtime]],
        drivetime_raw = [0],
        drivetime_rdd = [0],
        drivertrips = [0],
        drivermiles = [0],
        drivernightsaway = [0],
        drivernightshome = [0],
        driverjourneyssatisfyingcondition = [0],
        driverdayssatisfyingcondition = [0],
        driverutilization = [0]
    )

    if currtime == 0
        CSV.write(filename, df)
    else
        CSV.write(filename, df, append=true)
    end

end
  
#--------------------------------------------------------------------------#

function calculatedriverstatistics()

    totaldrivertrips, totaldriverhours_raw, totaldriverhours_rdd, totaldrivermiles, totalnightsaway, totalnightshome = 0, 0, 0, 0, 0, 0 
    for sgmt in pastdriversegments_space 
        totaldrivermiles += distbetweenlocs[sgmt[4],sgmt[5]]
        totaldriverhours_raw += traveltimebetweenlocs_raw[sgmt[4],sgmt[5]]
        totaldriverhours_rdd += traveltimebetweenlocs_rdd[sgmt[4],sgmt[5]]
        totaldrivertrips += 1
    end
    driverjourneyssatisfyingcondition, driverdayssatisfyingcondition = 0, 0

    driversegments = Dict()
    for d in drivers
        driversegments[d] = []
    end
    for sgmt in pastdriversegments
        push!(driversegments[sgmt[1]], sgmt)
    end
    T_off_Monday8am_0 = []
    for item in T_off_Monday8am
        newlist = []
        for hr in intersect(item, 0:tstep:onlinetimehorizon)
            if !(hr-tstep in item)
                push!(newlist, hr)
            end
        end
        push!(T_off_Monday8am_0, newlist)
    end
    for d in drivers
        hl = driverHomeLocs[d]
        homeeachnight = []
        for t in T_off_Monday8am_0[drivershift[d]]
            if (d,t,t+tstep,hl,hl) in driversegments[d]
                totalnightshome += 1
                push!(homeeachnight, 1)
            else
                totalnightsaway += 1
                push!(homeeachnight, 0)
            end
        end
        currnightsaway, journeynightsaway = 0, []
        while homeeachnight != []
            n = popfirst!(homeeachnight)
            if n == 1
                push!(journeynightsaway, currnightsaway)
                currnightsaway = 0
            else
                currnightsaway += 1
            end
        end
        push!(journeynightsaway, currnightsaway)
        standardjourneynightsaway = [j for j in journeynightsaway if j <= maxnightsaway]
        if standardjourneynightsaway != []
            driverjourneyssatisfyingcondition += length(standardjourneynightsaway) / length(journeynightsaway)
            driverdayssatisfyingcondition += sum(j+1 for j in standardjourneynightsaway) / sum(j+1 for j in journeynightsaway)
        end
    end

    driverjourneyssatisfyingcondition = driverjourneyssatisfyingcondition / length(drivers)
    driverdayssatisfyingcondition = driverdayssatisfyingcondition / length(drivers)

    return totaldrivertrips, totaldriverhours_raw, totaldriverhours_rdd, totaldrivermiles, totalnightsaway, totalnightshome, driverjourneyssatisfyingcondition, driverdayssatisfyingcondition

end

#--------------------------------------------------------------------------#

function writeresults_onlinefinal(filename, currtime)

    totaldrivertrips, totaldriverhours_raw, totaldriverhours_rdd, totaldrivermiles, totalnightsaway, totalnightshome, driverjourneyssatisfyingcondition, driverdayssatisfyingcondition = calculatedriverstatistics()

    df = DataFrame(
        experiment_id = [experiment_id],
        instance = [ex],
        weekstart = [weekstart],
        horizon = [horizon],
        tstep = [tstep],
        lambda = [lambda],
        numdrivers = [length(drivers)],
        operations = [operations],
        onlineweeks = [onlineweeks],
        maxnightsaway = [maxnightsaway],
        maxrepositioningdistance = [maxrepositioningdistance],
        method = [solutionmethod],
        currtime = [currtime],
        currorders = [length(currstate.orders)],
        totalorders = [highestorderindex],
        orderscompleted = [highestorderindex - length(currstate.orders)],
        totalcost = [sum(values(totalpastcost))],
        ordertrips = [sum(values(totalordertrips))],
        emptytrips = [sum(values(totalemptytrips))],
        taxitrips = [sum(values(totaltaxitrips))],
        ordermiles = [sum(values(totalordermiles))],
        penaltyordermiles = [sum(values(ordermilespenalty))],
        emptymiles = [sum(values(totalemptymiles))],
        taximiles = [sum(values(totaltaximiles))],
        shortesttripmiles = [sum(values(shortestpossible_ordermiles))],
        totaldeliverytime = [sum(values(total_delivtime))], 
        penaltydeliverytime = [sum(values(orderdelaypenalty))],
        maxdeliverytime = [maximum(values(max_delivtime))],
        shortestdeliverytime_abs = [0],
        shortestdeliverytime_rest = [sum(values(shortestpossible_delivtime))],
        drivetime_raw = [totaldriverhours_raw],
        drivetime_rdd = [totaldriverhours_rdd],
        drivertrips = [totaldrivertrips],
        drivermiles = [totaldrivermiles],
        drivernightsaway = [totalnightsaway],
        drivernightshome = [totalnightshome],
        driverjourneyssatisfyingcondition = [driverjourneyssatisfyingcondition],
        driverdayssatisfyingcondition = [driverdayssatisfyingcondition],
        driverutilization = [totaldriverhours_rdd / (length(drivers) * 0.5 * onlinetimehorizon)]
    )
    
    CSV.write(filename, df, append=true)

end

#--------------------------------------------------------------------------#

function writedeliverytimetofile(filename, currtime, i, actualdelivtime, bestdelivtime)

    df = DataFrame(
        experiment_id = [experiment_id],
        currtime = [currtime],
        ordernumber = [i],
        totaldeliverytime = [actualdelivtime],
        shortestdeliverytime = [bestdelivtime]
    )
    
    CSV.write(filename, df, append=true)

end