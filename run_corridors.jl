
using CSV, Luxor, Colors, Random, DataFrames, Dates, StatsBase, JuMP, Gurobi

include("scripts/instancegeneration/readrivigodata.jl")

#--------------------------------------------------------------------------------------------------#

ORIENTATION = "NS"
target_ei, target_ni, target_gi = 0.9, 0.3, 0
weeks = 4
N = 800 * weeks
n = N
outputfilename = "orders_$(target_ei)_$(target_ni)_$(target_gi).csv"

if ORIENTATION == "EW"
    #East/West
    #group1 = [63,65,66,58,52,41,42,46,51,59,49,48,57,34,38, 58] #East
    #group2 = [50, 56, 55, 53, 28, 27, 31, 36, 40, 39, 35, 32, 29, 26, 45] #West
    group1 = [63, 65, 66, 58, 52, 41, 42, 46, 51]
    group2 = [56, 55, 53, 28, 27, 31, 36, 40, 39, 35, 32, 29, 26, 45]
    loclist = union(group1, group2)
    region1 = [56, 55, 53, 63, 65, 66, 51, 46, 45, 58] #North
    region2 = [58, 52, 41, 42, 31, 36, 40, 39, 28, 27, 35, 32, 29, 26, 58] #South
    DIVIDEINTOREGIONS = false

elseif ORIENTATION == "NS"
    #North/South
    group1 = [64, 61, 60, 62, 56, 55, 53, 54, 45]
    group2 = [25, 19, 20, 17, 16, 13, 14, 8, 7, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12]
    loclist = [64, 61, 60, 62, 56, 55, 53, 54, 45, 37, 30, 25, 19, 20, 17, 16, 13, 14, 8, 7, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12]
    DIVIDEINTOREGIONS = false

elseif ORIENTATION == "NS2"
    #North/South
    group1 = [60, 45, 54, 36]
    group2 = [26, 27, 15, 1, 10, 7, 9]
    loclist = union(group1, group2)
    DIVIDEINTOREGIONS = false

elseif ORIENTATION == "all"
    #All
    group1 = 1:numlocs
    group2 = 1:numlocs
    loclist = 1:numlocs
    DIVIDEINTOREGIONS = false

end

#--------------------------------------------------------------------------------------------------#

thickest, thinnest = 20, 3
pixelshift = 22
maxlocs = 66
tstep = 6
roundup_flag = 1
excludeoutliers_flag = 1
googlemapstraveltimes_flag = 1
includesymmetricarcs_flag = 1
ensureconnectivity_flag = 1
shiftlength = 12
lhdataisbfilename = "data/lh_data_isb_connect_clean.csv"
hubdistancesfilename = "data/hubdistances.csv"
hubdataisbfilename = "data/hub_data_isb_connect.csv"
traveltimesfilename = "data/traveltimes_outliers.csv"
operations = "relay"

hubCoords, hubsLookup, hubsReverseLookup, hubsTravelTimeIndex, numlocs = readlocations(hubdataisbfilename, maxlocs)
prearcs, arcLength, arcLength_raw = readandprocessarcs(operations, traveltimesfilename, hubdistancesfilename, tstep, numlocs, hubsTravelTimeIndex, roundup_flag, excludeoutliers_flag, hubsReverseLookup, googlemapstraveltimes_flag)
distbetweenlocs, shortestpatharclists = cacheShortestDistance(numlocs, prearcs)

#--------------------------------------------------------------------------------------------------#

#Pulls all historical Rivigo data
function readrivigoorders(lhdataisbfilename, loclist)

    data_agg = CSV.read(lhdataisbfilename, DataFrame)
    rivigodata = DataFrame(id=[], dayofweek=[], timeofday=[], origin=[], destination=[], stopsequence=[])

    hubsList = collect(values(hubsLookup))
    for i in 1:size(data_agg)[1]
        orig, dest = data_agg[!, 26][i], data_agg[!, 27][i]
        psseq_raw = data_agg[i, 8]
        orderdatetime = DateTime(1970, 1, 1, 0, 0, 0) + Millisecond(data_agg[i, "departure_timestamp"]) #dates in csv are in milliseconds since 1/1/1970, lol
        orderdayofweek = dayofweek(orderdatetime)
        ordertimeofday = Time(orderdatetime)
        psseq = split(psseq_raw, "-")

        #Check whether all intermediate nodes from the Rivigo pitstop sequence are included in the subset of locs
        intermedlocs_flag = 0
        stopsequence = []
        for ps in psseq
            if ps in hubsList
                loc = hubsReverseLookup[ps]
                if loc > numlocs
                    intermedlocs_flag = 1
                    break
                else
                    push!(stopsequence, loc)
                end
            else
                intermedlocs_flag = 1
                break
            end
        end

        #if (orig != dest) & (orig in loclist) & (dest in loclist) & (intermedlocs_flag == 0) #& (orderwindowstart <= start_avail_ts <= orderwindowend) 
        #    push!(rivigodata,[data_agg[!,1][i], orderdayofweek, ordertimeofday, orig, dest, stopsequence])
        #end
        if !DIVIDEINTOREGIONS
            if (orig != dest) & (orig in group1) & (dest in group2) & (intermedlocs_flag == 0) #& (orderwindowstart <= start_avail_ts <= orderwindowend) 
                push!(rivigodata, [data_agg[!, 1][i], orderdayofweek, ordertimeofday, orig, dest, stopsequence])
            elseif (orig != dest) & (orig in group2) & (dest in group1) & (intermedlocs_flag == 0) #& (orderwindowstart <= start_avail_ts <= orderwindowend) 
                push!(rivigodata, [data_agg[!, 1][i], orderdayofweek, ordertimeofday, orig, dest, stopsequence])
            end
        else
            if (orig != dest) & (orig in group1) & (dest in group2) & (intermedlocs_flag == 0) & (orig in region1) & (dest in region1) #& (orderwindowstart <= start_avail_ts <= orderwindowend) 
                push!(rivigodata, [data_agg[!, 1][i], orderdayofweek, ordertimeofday, orig, dest, stopsequence])
            elseif (orig != dest) & (orig in group2) & (dest in group1) & (intermedlocs_flag == 0) & (orig in region1) & (dest in region1) #& (orderwindowstart <= start_avail_ts <= orderwindowend) 
                push!(rivigodata, [data_agg[!, 1][i], orderdayofweek, ordertimeofday, orig, dest, stopsequence])
            elseif (orig != dest) & (orig in group1) & (dest in group2) & (intermedlocs_flag == 0) & (orig in region2) & (dest in region2) #& (orderwindowstart <= start_avail_ts <= orderwindowend) 
                push!(rivigodata, [data_agg[!, 1][i], orderdayofweek, ordertimeofday, orig, dest, stopsequence])
            elseif (orig != dest) & (orig in group2) & (dest in group1) & (intermedlocs_flag == 0) & (orig in region2) & (dest in region2) #& (orderwindowstart <= start_avail_ts <= orderwindowend) 
                push!(rivigodata, [data_agg[!, 1][i], orderdayofweek, ordertimeofday, orig, dest, stopsequence])
            end
        end

    end

    return rivigodata

end

#--------------------------------------------------------------------------------------------------#

function calculatebalance(ordersbetween)

    edgeimbalance, nodeimbalance, totalflow = 0, 0, 0
    for i in group1, j in group2
        totalflow += ordersbetween[i, j] + ordersbetween[j, i]
        edgeimbalance += abs(ordersbetween[i, j] - ordersbetween[j, i])
    end
    for i in group1
        nodeimbalance += abs(sum(ordersbetween[i, j] for j in group2) - sum(ordersbetween[j, i] for j in group2))
    end
    for i in group2
        nodeimbalance += abs(sum(ordersbetween[i, j] for j in group1) - sum(ordersbetween[j, i] for j in group1))
    end
    globalimbalance = abs(sum(sum(ordersbetween[i, j] for j in group2) for i in group1) - sum(sum(ordersbetween[j, i] for j in group2) for i in group1))
    println("Edge imbalance = ", edgeimbalance / totalflow)
    println("Node imbalance = ", nodeimbalance / totalflow / 2)
    println("Global imbalance = ", globalimbalance / totalflow)
    println("Total flow = ", totalflow)

    return edgeimbalance / totalflow, nodeimbalance / totalflow / 2, globalimbalance / totalflow

end

#--------------------------------------------------------------------------------------------------#

function restorebalance(N, target_ei, target_ni, target_gi, ordersbetween_temp)

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "OutputFlag", 1)
    set_optimizer_attribute(model, "TimeLimit", 120)
    @variable(model, x[loclist, loclist] >= 0, Int)
    @variable(model, diff_ei[group1, group2] >= 0)
    @variable(model, diff_ni[loclist] >= 0)
    @variable(model, diff_gi >= 0)
    @variable(model, z_ei[group1, group2], Bin)
    @variable(model, z_ni[loclist], Bin)
    @variable(model, z_gi, Bin)
    @variable(model, ei_diff >= 0)
    @variable(model, ni_diff >= 0)
    @variable(model, gi_diff >= 0)
    @objective(model, Min, 100 * ei_diff + ni_diff) #diff_gi + sum(diff_ni[i] for i in loclist) + sum(sum(diff_ei[i,j] for j in group2) for i in group1) )
    #@objective(model, Min, sum(diff_ni[i] for i in loclist) ) 
    @constraint(model, demand[i in loclist, j in loclist], x[i, j] >= ordersbetween_temp[i, j])
    @constraint(model, totalflow[i in loclist, j in loclist], sum(sum(x[i, j] + x[j, i] for j in group2) for i in group1) == N)
    @constraint(model, edgediff1[i in group1, j in group2], diff_ei[i, j] >= x[i, j] - x[j, i])
    @constraint(model, edgediff2[i in group1, j in group2], diff_ei[i, j] >= x[j, i] - x[i, j])
    @constraint(model, edgediff3[i in group1, j in group2], diff_ei[i, j] <= x[i, j] - x[j, i] + N * z_ei[i, j])
    @constraint(model, edgediff4[i in group1, j in group2], diff_ei[i, j] <= x[j, i] - x[i, j] + N * (1 - z_ei[i, j]))
    @constraint(model, nodediff1[i in group1], diff_ni[i] >= sum(x[i, j] - x[j, i] for j in group2))
    @constraint(model, ndoediff2[i in group1], diff_ni[i] >= sum(x[j, i] - x[i, j] for j in group2))
    @constraint(model, nodediff3[i in group2], diff_ni[i] >= sum(x[i, j] - x[j, i] for j in group1))
    @constraint(model, ndoediff4[i in group2], diff_ni[i] >= sum(x[j, i] - x[i, j] for j in group1))
    @constraint(model, nodediff5[i in group1], diff_ni[i] <= sum(x[i, j] - x[j, i] for j in group2) + N * z_ni[i])
    @constraint(model, ndoediff6[i in group1], diff_ni[i] <= sum(x[j, i] - x[i, j] for j in group2) + N * (1 - z_ni[i]))
    @constraint(model, nodediff7[i in group2], diff_ni[i] <= sum(x[i, j] - x[j, i] for j in group1) + N * z_ni[i])
    @constraint(model, ndoediff8[i in group2], diff_ni[i] <= sum(x[j, i] - x[i, j] for j in group1) + N * (1 - z_ni[i]))
    @constraint(model, globaldiff1, diff_gi >= sum(sum(x[i, j] - x[j, i] for j in group2) for i in group1))
    @constraint(model, globaldiff2, diff_gi >= sum(sum(x[j, i] - x[i, j] for j in group2) for i in group1))
    @constraint(model, globaldiff3, diff_gi <= sum(sum(x[i, j] - x[j, i] for j in group2) for i in group1) + N * z_gi)
    @constraint(model, globaldiff4, diff_gi <= sum(sum(x[j, i] - x[i, j] for j in group2) for i in group1) + N * (1 - z_gi))
    @constraint(model, ei_diff >= sum(sum(diff_ei[i, j] for j in group2) for i in group1) - N * target_ei)
    @constraint(model, ei_diff >= N * target_ei - sum(sum(diff_ei[i, j] for j in group2) for i in group1))
    @constraint(model, ni_diff >= 2 * N * target_ni - sum(diff_ni[i] for i in loclist))
    @constraint(model, ni_diff >= sum(diff_ni[i] for i in loclist) - 2 * N * target_ni)
    #@constraint(model, gi_diff >= diff_gi - N * target_gi)
    #@constraint(model, gi_diff >= N * target_gi - diff_gi)
    @constraint(model, N * target_gi == diff_gi)
    #@constraint(model, sum(sum(diff_ei[i,j] for j in group2) for i in group1) == N * target_ei)

    optimize!(model)

    additionaltrips, checkbalance = Dict(), Dict()
    for i in group1, j in group2
        additionaltrips[i, j] = value(x[i, j]) - ordersbetween_temp[i, j]
        additionaltrips[j, i] = value(x[j, i]) - ordersbetween_temp[j, i]
        #checkbalance[i,j] = value(x[i,j])
        #checkbalance[j,i] = value(x[j,i])
        #@assert(value(x[i,j]) >= ordersbetween_temp[i,j])
        #@assert(value(x[j,i]) >= ordersbetween_temp[j,i])
    end
    println("Targets = $target_ei, $target_ni, $target_gi")
    #calculatebalance(checkbalance)
    #calculatebalance(ordersbetween_temp)
    #@assert(N == sum(sum(additionaltrips[i,j] for j in group2) for i in group1) + sum(sum(additionaltrips[j,i] for j in group2) for i in group1))

    return additionaltrips

end

#--------------------------------------------------------------------------------------------------#

function restorebalance_changedemand(N, target_ei, target_ni, target_gi, ordersbetween_temp)

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "OutputFlag", 1)
    set_optimizer_attribute(model, "TimeLimit", 60)
    @variable(model, x[loclist, loclist] >= 0, Int)
    @variable(model, diff_ei[group1, group2] >= 0)
    @variable(model, diff_ni[loclist] >= 0)
    @variable(model, diff_gi >= 0)
    @variable(model, z_ei[group1, group2], Bin)
    @variable(model, z_ni[loclist], Bin)
    @variable(model, z_gi, Bin)
    @variable(model, demand_diff[loclist, loclist] >= 0, Int)
    @objective(model, Min, sum(sum(demand_diff[i, j] + demand_diff[j, i] for j in group2) for i in group1))
    #@constraint(model, demand[i in loclist, j in loclist], x[i,j] >= ordersbetween_temp[i,j])
    @constraint(model, totalflow[i in loclist, j in loclist], sum(sum(x[i, j] + x[j, i] for j in group2) for i in group1) == N)
    @constraint(model, demanddiff1[i in group1, j in group2], demand_diff[i, j] >= x[i, j] - ordersbetween_temp[i, j])
    @constraint(model, demanddiff2[i in group1, j in group2], demand_diff[i, j] >= ordersbetween_temp[i, j] - x[i, j])
    @constraint(model, demanddiff3[i in group1, j in group2], demand_diff[i, j] >= x[j, i] - ordersbetween_temp[j, i])
    @constraint(model, demanddiff4[i in group1, j in group2], demand_diff[i, j] >= ordersbetween_temp[j, i] - x[j, i])
    @constraint(model, edgediff1[i in group1, j in group2], diff_ei[i, j] >= x[i, j] - x[j, i])
    @constraint(model, edgediff2[i in group1, j in group2], diff_ei[i, j] >= x[j, i] - x[i, j])
    @constraint(model, edgediff3[i in group1, j in group2], diff_ei[i, j] <= x[i, j] - x[j, i] + N * z_ei[i, j])
    @constraint(model, edgediff4[i in group1, j in group2], diff_ei[i, j] <= x[j, i] - x[i, j] + N * (1 - z_ei[i, j]))
    @constraint(model, nodediff1[i in group1], diff_ni[i] >= sum(x[i, j] - x[j, i] for j in group2))
    @constraint(model, ndoediff2[i in group1], diff_ni[i] >= sum(x[j, i] - x[i, j] for j in group2))
    @constraint(model, nodediff3[i in group2], diff_ni[i] >= sum(x[i, j] - x[j, i] for j in group1))
    @constraint(model, ndoediff4[i in group2], diff_ni[i] >= sum(x[j, i] - x[i, j] for j in group1))
    @constraint(model, nodediff5[i in group1], diff_ni[i] <= sum(x[i, j] - x[j, i] for j in group2) + N * z_ni[i])
    @constraint(model, ndoediff6[i in group1], diff_ni[i] <= sum(x[j, i] - x[i, j] for j in group2) + N * (1 - z_ni[i]))
    @constraint(model, nodediff7[i in group2], diff_ni[i] <= sum(x[i, j] - x[j, i] for j in group1) + N * z_ni[i])
    @constraint(model, ndoediff8[i in group2], diff_ni[i] <= sum(x[j, i] - x[i, j] for j in group1) + N * (1 - z_ni[i]))
    @constraint(model, globaldiff1, diff_gi >= sum(sum(x[i, j] - x[j, i] for j in group2) for i in group1))
    @constraint(model, globaldiff2, diff_gi >= sum(sum(x[j, i] - x[i, j] for j in group2) for i in group1))
    @constraint(model, globaldiff3, diff_gi <= sum(sum(x[i, j] - x[j, i] for j in group2) for i in group1) + N * z_gi)
    @constraint(model, globaldiff4, diff_gi <= sum(sum(x[j, i] - x[i, j] for j in group2) for i in group1) + N * (1 - z_gi))
    @constraint(model, N * target_gi == diff_gi)
    @constraint(model, sum(sum(diff_ei[i, j] for j in group2) for i in group1) == N * target_ei)
    @constraint(model, 2 * N * target_ni == sum(diff_ni[i] for i in loclist))

    optimize!(model)

    println("Achieved edge imbalance = ", sum(sum(value(diff_ei[i, j]) for j in group2) for i in group1) / N)
    println("Achieved node imbalance = ", sum(value(diff_ni[i]) for i in loclist) / (2 * N))
    println("Achieved global imbalance = ", value(diff_gi) / N)

    additionaltrips, checkbalance = Dict(), Dict()
    for i in group1, j in group2
        additionaltrips[i, j] = convert(Int64, value(x[i, j]) - ordersbetween_temp[i, j])
        additionaltrips[j, i] = convert(Int64, value(x[j, i]) - ordersbetween_temp[j, i])
    end
    println("Targets = $target_ei, $target_ni, $target_gi")

    return additionaltrips

end

#--------------------------------------------------------------------------------------------------#

#Pulls one week worth of Rivigo data, randomly sampled according to parameters: target_ei, target_ni, target_gi
function pullrivigosample(lhdataisbfilename, loclist, target_ei, target_ni, target_gi)

    data_agg = CSV.read(lhdataisbfilename, DataFrame)
    empirical_times_hours = Hour.(Millisecond.(data_agg[:, "departure_timestamp"]) .+ DateTime(1970, 1, 1, 0, 0, 0))
    empirical_times_dayofweek = dayofweek.(Millisecond.(data_agg[:, "departure_timestamp"]) .+ DateTime(1970, 1, 1, 0, 0, 0))

    rivigodata = readrivigoorders(lhdataisbfilename, loclist)

    #Sample data
    sampleddata = DataFrame(id=[], dayofweek=[], timeofday=[], origin=[], destination=[], stopsequence=[])
    for i in 1:n
        push!(sampleddata, rivigodata[rand(1:size(rivigodata)[1]), :])
    end
    ordersbetween_temp = Dict()
    for i in 1:numlocs, j in 1:numlocs
        ordersbetween_temp[i, j] = 0
    end
    for i in 1:size(sampleddata)[1]
        orig, dest, stopsequence = sampleddata[i, "origin"], sampleddata[i, "destination"], sampleddata[i, "stopsequence"]
        ordersbetween_temp[orig, dest] += 1
    end

    additionaltrips = restorebalance_changedemand(N, target_ei, target_ni, target_gi, ordersbetween_temp)

    for i in group1, j in group2
        if additionaltrips[i, j] > 1e-4
            filtereddata = filter(row -> (row.origin == i) & (row.destination == j), rivigodata)
            for newrow in 1:additionaltrips[i, j]
                try
                    push!(sampleddata, filtereddata[rand(1:size(filtereddata)[1]), :])
                catch
                    push!(sampleddata, [-1, rand(empirical_times_dayofweek), Time(Dates.value(rand(empirical_times_hours)), rand(0:59), 0), i, j, []])
                end
            end
        elseif additionaltrips[i, j] < -1e-4
            filtereddata = filter(row -> (row.origin == i) & (row.destination == j), sampleddata)
            subsetdata = filtereddata[randperm(size(filtereddata)[1])[1:-1*additionaltrips[i,j]], :][:,"id"]
            for item in subsetdata
                deleteat!(sampleddata, findfirst(in([item]),sampleddata.id))
            end
        end
    end
    for i in group2, j in group1
        if additionaltrips[i, j] > 1e-4
            filtereddata = filter(row -> (row.origin == i) & (row.destination == j), rivigodata)
            for newrow in 1:additionaltrips[i, j]
                try
                    push!(sampleddata, filtereddata[rand(1:size(filtereddata)[1]), :])
                catch
                    push!(sampleddata, [-1, rand(empirical_times_dayofweek), Time(Dates.value(rand(empirical_times_hours)), rand(0:59), 0), i, j, []])
                end
            end
        elseif additionaltrips[i, j] < -1e-4
            filtereddata = filter(row -> (row.origin == i) & (row.destination == j), sampleddata)
            subsetdata = filtereddata[randperm(size(filtereddata)[1])[1:-1*additionaltrips[i,j]], :][:,"id"]
            for item in subsetdata
                deleteat!(sampleddata, findfirst(in([item]),sampleddata.id))
            end
        end
    end

    #Get sets for visualization 
    tripson, group1tripson, group2tripson, origincount, destinationcount, ordersbetween = Dict(), Dict(), Dict(), Dict(), Dict(), Dict()
    for i in 1:numlocs, j in 1:numlocs
        tripson[i, j] = 0
        group1tripson[i, j] = 0
        group2tripson[i, j] = 0
        ordersbetween[i, j] = 0
    end
    for i in 1:numlocs
        origincount[i] = 0
        destinationcount[i] = 0
    end

    numpitstops, tripdists = [], []

    for i in 1:size(sampleddata)[1]
        orig, dest, stopsequence = sampleddata[i, "origin"], sampleddata[i, "destination"], sampleddata[i, "stopsequence"]
        ordersbetween[orig, dest] += 1
        tripdist = 0
        for i in 1:length(stopsequence)-1
            tripson[stopsequence[i], stopsequence[i+1]] += 1
            tripdist += distbetweenlocs[stopsequence[i], stopsequence[i+1]]
            if orig in group1
                group1tripson[stopsequence[i], stopsequence[i+1]] += 1
            else
                group2tripson[stopsequence[i], stopsequence[i+1]] += 1
            end
        end
        origincount[orig] += 1
        destinationcount[dest] += 1
        push!(numpitstops, length(stopsequence))
        push!(tripdists, tripdist)
    end

    eb, nb, gb = calculatebalance(ordersbetween)

    return sampleddata, tripson, origincount, destinationcount, ordersbetween, group1tripson, group2tripson

end

#--------------------------------------------------------------------------------------------------#

function spatialnetwork_downsampled(drawingname, drawingname2, lhdataisbfilename, xdim, ydim, loclist, tripson, ordersbetween, group1tripson, group2tripson)

    #Get correct scale
    maxlat, minlat = 0, 100
    for l in 1:numlocs
        maxlat = max(hubCoords[l, 1], maxlat)
        minlat = min(hubCoords[l, 1], minlat)
    end

    latmult = -(xdim - 200) / (maxlat - minlat)
    latshift = -(xdim - 200) / 2 + (xdim - 200) * maxlat / (maxlat - minlat)
    longmult = -1 * latmult * 24 / 29
    maxlongcoord, minlongcoord = -100000, 100000
    for l in 1:numlocs
        maxlongcoord = max(longmult * hubCoords[l, 2], maxlongcoord)
        minlongcoord = min(longmult * hubCoords[l, 2], minlongcoord)
    end
    longshift = -(maxlongcoord + minlongcoord) / 2

    #Format and transform latitude and longitude coordinates of each pit stop
    pointDict = Dict()
    listofpoints = []
    listofpoints_labels = []
    for l in 1:numlocs
        longitude, latitude = hubCoords[l, 2], hubCoords[l, 1]
        transformedcoords = (longmult * longitude + longshift, latmult * latitude + latshift)
        pointDict[l] = Point(transformedcoords)
        push!(listofpoints, transformedcoords)
        push!(listofpoints_labels, [transformedcoords, string(l)])
    end
    locationPoints = Point.(listofpoints)

    #--------------------------------------------------------#

    #Calculate thickness of each arc
    mintrips, maxtrips = 1, maximum(values(tripson))
    arcList = []
    for i in 1:numlocs, j in setdiff(1:numlocs, i)
        if group1tripson[i, j] >= 1
            startPoint = locationPoints[i]
            endPoint = locationPoints[j]
            thickness = round(thinnest + (tripson[i, j] - mintrips) / (maxtrips - mintrips) * (thickest - thinnest))
            if ORIENTATION[1:2] in ["NS", "al"]
                push!(arcList, (startPoint + Point(10, 0), endPoint + Point(10, 0), (200, 0, 0), thickness, "solid"))
            elseif ORIENTATION[1:2] == "EW"
                push!(arcList, (startPoint + Point(0, 10), endPoint + Point(0, 10), (200, 0, 0), thickness, "solid"))
            end
        end
        if group2tripson[i, j] >= 1
            startPoint = locationPoints[i]
            endPoint = locationPoints[j]
            thickness = round(thinnest + (tripson[i, j] - mintrips) / (maxtrips - mintrips) * (thickest - thinnest))
            if ORIENTATION[1:2] in ["NS", "al"]
                push!(arcList, (startPoint - Point(10, 0), endPoint - Point(10, 0), (0, 0, 200), thickness, "solid"))
            elseif ORIENTATION[1:2] == "EW"
                push!(arcList, (startPoint - Point(0, 10), endPoint - Point(0, 10), (0, 0, 200), thickness, "solid"))
            end
        end
    end

    #=mintrips, maxtrips = 1, maximum(values(ordersbetween))
    arcList = []	
    for i in 1:numlocs, j in setdiff(1:numlocs, i)
        if (i in group1) & (ordersbetween[i,j] >= 1)
            startPoint = locationPoints[i]
            endPoint = locationPoints[j]
            thickness = round(thinnest + (ordersbetween[i,j] - mintrips)/(maxtrips - mintrips) * (thickest - thinnest) )
            push!(arcList, (startPoint + Point(10,0), endPoint + Point(10,0), (200,0,0), thickness, "solid"))
        end
        if (i in group2) & (ordersbetween[i,j] >= 1)
            startPoint = locationPoints[i]
            endPoint = locationPoints[j]
            thickness = round(thinnest + (ordersbetween[i,j] - mintrips)/(maxtrips - mintrips) * (thickest - thinnest) )
            push!(arcList, (startPoint - Point(10,0), endPoint - Point(10,0), (0,0,200), thickness, "solid"))
        end
    end=#

    #--------------------------------------------------------#

    #Create new drawing
    Drawing(xdim, ydim, drawingname)
    origin()
    background("white")

    #Draw the arcs
    for i in arcList
        #Set arc attributes
        setline(i[4])
        setcolor(i[3])
        setdash(i[5])

        #Draw the arc line
        line(i[1], i[2], :stroke)

        #Calculate the rotation and placement of the arrowhead
        theta = atan((i[2][2] - i[1][2]) / (i[2][1] - i[1][1]))
        dist = distance(i[1], i[2])
        arrowhead = (1 - pixelshift / dist) * i[2] + (pixelshift / dist) * i[1] #center of arrowhead positioned 8 pixels from the end node

        #Rotate the arrowhead appropriately
        if i[1][1] >= i[2][1]
            local p = ngon(arrowhead, min(pixelshift, i[4] * 2), 3, theta - pi, vertices=true)
        else
            local p = ngon(arrowhead, min(pixelshift, i[4] * 2), 3, theta, vertices=true)
        end

        #Draw the arrowhead
        poly(p, :fill, close=true)
    end

    #Draw the pit stop nodes
    setcolor("black")
    circle.(locationPoints, 16, :fill)
    setcolor("black")
    setline(3)
    circle.(locationPoints, 16, :stroke)

    #Add pit stop labels
    fontsize(22)
    setcolor("white")
    for item in listofpoints_labels
        #label(item[2], :0, Point(item[1]))
        Luxor.text(item[2], Point(item[1]), halign=:center, valign=:middle)
    end
    setcolor("black")

    #Legend box
    #=setline(4)
    legendstartx = 0.5*xdim - 0.43*xdim
    legendstarty = 0.5*ydim - 0.3*ydim
    rect(legendstartx, legendstarty, 0.4*xdim, 0.25*ydim, :stroke)

    #Arcs for the legend
    fontsize(70)
    numlegendarcs = 4
    meantrips = convert(Int,round(mean([k for k in values(tripson) if k > 0]), digits=0))
    trips90 = convert(Int,round(percentile([k for k in values(tripson) if k > 0], 90), digits=0))
    legendthicknesses = [mintrips, meantrips, trips90, maxtrips]
    legendlabels = ["$mintrips trip (min)", "$meantrips trips (mean)", "$trips90 trips (p90)", "$maxtrips trips (max)"]
    for legendarc in 1:numlegendarcs
        startPoint = Point(legendstartx + 0.03*xdim, legendstarty + (legendarc-0.5)/numlegendarcs * 0.25*ydim)
        endPoint = startPoint + Point(xdim/20, 0)
        thickness = round(thinnest + (legendthicknesses[legendarc] - mintrips)/(maxtrips - mintrips) * (thickest - thinnest) )

        #Draw the arc line
        setline(thickness)
    line(startPoint, endPoint , :stroke)

    #Calculate the rotation and placement of the arrowhead
    theta = atan((endPoint[2] - startPoint[2])/(endPoint[1] - startPoint[1]))
    dist = distance(startPoint, endPoint)
    arrowhead = (1-0/dist)*endPoint + (0/dist)*startPoint #center of arrowhead positioned 8 pixels from the end node

    #Rotate the arrowhead appropriately
    if startPoint[1] >= endPoint[1]
    local p = ngon(arrowhead, min(pixelshift, thickness*2), 3, theta - pi , vertices=true)
    else
    local p = ngon(arrowhead, min(pixelshift, thickness*2), 3, theta , vertices=true)
    end

    #Draw the arrowhead
    poly(p, :fill,  close=true)

        #Add the label
        label(legendlabels[legendarc], :E , endPoint + Point(xdim/40, 0))
    end=#

    #--------------------------------------------------------#

    finish()
    preview()

    #--------------------------------------------------------#

    #Calculate thickness of each arc
    #=mintrips, maxtrips = 1, maximum(values(tripson))
    arcList = []	
    for i in 1:numlocs, j in setdiff(1:numlocs, i)
        if group1tripson[i,j] >= 1
            startPoint = locationPoints[i]
            endPoint = locationPoints[j]
            thickness = round(thinnest + (tripson[i,j] - mintrips)/(maxtrips - mintrips) * (thickest - thinnest) )
            push!(arcList, (startPoint + Point(10,0), endPoint + Point(10,0), (200,0,0), thickness, "solid"))
        end
        if group2tripson[i,j] >= 1
            startPoint = locationPoints[i]
            endPoint = locationPoints[j]
            thickness = round(thinnest + (tripson[i,j] - mintrips)/(maxtrips - mintrips) * (thickest - thinnest) )
            push!(arcList, (startPoint - Point(10,0), endPoint - Point(10,0), (0,0,200), thickness, "solid"))
        end
    end=#

    mintrips, maxtrips = 1, maximum(values(ordersbetween))
    arcList = []
    for i in 1:numlocs, j in setdiff(1:numlocs, i)
        if (i in group1) & (ordersbetween[i, j] >= 1)
            startPoint = locationPoints[i]
            endPoint = locationPoints[j]
            thickness = round(thinnest + (ordersbetween[i, j] - mintrips) / (maxtrips - mintrips) * (thickest - thinnest))
            if ORIENTATION[1:2] in ["NS", "al"]
                push!(arcList, (startPoint + Point(10, 0), endPoint + Point(10, 0), (200, 0, 0), thickness, "solid"))
            elseif ORIENTATION[1:2] == "EW"
                push!(arcList, (startPoint + Point(0, 10), endPoint + Point(0, 10), (200, 0, 0), thickness, "solid"))
            end
        end
        if (i in group2) & (ordersbetween[i, j] >= 1)
            startPoint = locationPoints[i]
            endPoint = locationPoints[j]
            thickness = round(thinnest + (ordersbetween[i, j] - mintrips) / (maxtrips - mintrips) * (thickest - thinnest))
            if ORIENTATION[1:2] in ["NS", "al"]
                push!(arcList, (startPoint - Point(10, 0), endPoint - Point(10, 0), (0, 0, 200), thickness, "solid"))
            elseif ORIENTATION[1:2] == "EW"
                push!(arcList, (startPoint - Point(0, 10), endPoint - Point(0, 10), (0, 0, 200), thickness, "solid"))
            end
        end
    end

    #--------------------------------------------------------#

    #Create new drawing
    Drawing(xdim, ydim, drawingname2)
    origin()
    background("white")

    #Draw the arcs
    for i in arcList
        #Set arc attributes
        setline(i[4])
        setcolor(i[3])
        setdash(i[5])

        #Draw the arc line
        line(i[1], i[2], :stroke)

        #Calculate the rotation and placement of the arrowhead
        theta = atan((i[2][2] - i[1][2]) / (i[2][1] - i[1][1]))
        dist = distance(i[1], i[2])
        arrowhead = (1 - pixelshift / dist) * i[2] + (pixelshift / dist) * i[1] #center of arrowhead positioned 8 pixels from the end node

        #Rotate the arrowhead appropriately
        if i[1][1] >= i[2][1]
            local p = ngon(arrowhead, min(pixelshift, i[4] * 2), 3, theta - pi, vertices=true)
        else
            local p = ngon(arrowhead, min(pixelshift, i[4] * 2), 3, theta, vertices=true)
        end

        #Draw the arrowhead
        poly(p, :fill, close=true)
    end

    #Draw the pit stop nodes
    setcolor("black")
    circle.(locationPoints, 16, :fill)
    setcolor("black")
    setline(3)
    circle.(locationPoints, 16, :stroke)

    #Add pit stop labels
    fontsize(22)
    setcolor("white")
    for item in listofpoints_labels
        #label(item[2], :0, Point(item[1]))
        Luxor.text(item[2], Point(item[1]), halign=:center, valign=:middle)
    end
    setcolor("black")

    #Legend box
    #=setline(4)
    legendstartx = 0.5*xdim - 0.43*xdim
    legendstarty = 0.5*ydim - 0.3*ydim
    rect(legendstartx, legendstarty, 0.4*xdim, 0.25*ydim, :stroke)

    #Arcs for the legend
    fontsize(70)
    numlegendarcs = 4
    meantrips = convert(Int,round(mean([k for k in values(tripson) if k > 0]), digits=0))
    trips90 = convert(Int,round(percentile([k for k in values(tripson) if k > 0], 90), digits=0))
    legendthicknesses = [mintrips, meantrips, trips90, maxtrips]
    legendlabels = ["$mintrips trip (min)", "$meantrips trips (mean)", "$trips90 trips (p90)", "$maxtrips trips (max)"]
    for legendarc in 1:numlegendarcs
        startPoint = Point(legendstartx + 0.03*xdim, legendstarty + (legendarc-0.5)/numlegendarcs * 0.25*ydim)
        endPoint = startPoint + Point(xdim/20, 0)
        thickness = round(thinnest + (legendthicknesses[legendarc] - mintrips)/(maxtrips - mintrips) * (thickest - thinnest) )

        #Draw the arc line
        setline(thickness)
    line(startPoint, endPoint , :stroke)

    #Calculate the rotation and placement of the arrowhead
    theta = atan((endPoint[2] - startPoint[2])/(endPoint[1] - startPoint[1]))
    dist = distance(startPoint, endPoint)
    arrowhead = (1-0/dist)*endPoint + (0/dist)*startPoint #center of arrowhead positioned 8 pixels from the end node

    #Rotate the arrowhead appropriately
    if startPoint[1] >= endPoint[1]
    local p = ngon(arrowhead, min(pixelshift, thickness*2), 3, theta - pi , vertices=true)
    else
    local p = ngon(arrowhead, min(pixelshift, thickness*2), 3, theta , vertices=true)
    end

    #Draw the arrowhead
    poly(p, :fill,  close=true)

        #Add the label
        label(legendlabels[legendarc], :E , endPoint + Point(xdim/40, 0))
    end=#

    #--------------------------------------------------------#

    finish()
    preview()

end

function writeorderfile(filename, orderdata)

    df = DataFrame(id=[],
        departure_timestamp=[],
        arrival_timestamp=[],
        ps_seq=[],
        Origin_PS=[],
        Destination_PS=[])

    currindex = 1
    for row in 1:size(orderdata)[1]
        push!(df, [currindex,
            Date(2026, 1, 4) + Day(orderdata[row, "dayofweek"] + 7*rand(0:weeks-1)) + orderdata[row, "timeofday"],
            Date(2026, 1, 4) + Day(orderdata[row, "dayofweek"] + 7*rand(0:weeks-1)) + orderdata[row, "timeofday"] + Millisecond(round(1000 * 3600 * 24 * 1.5 * distbetweenlocs[48, 61] / 500, digits=0)),
            orderdata[row, "stopsequence"],
            orderdata[row, "origin"],
            orderdata[row, "destination"]])
        currindex += 1
    end

    sort!(df, :departure_timestamp)

    CSV.write(filename, df)

end

sampleddata, tripson, origincount, destinationcount, ordersbetween, group1tripson, group2tripson = pullrivigosample(lhdataisbfilename, loclist, target_ei, target_ni, target_gi)
#spatialnetwork_downsampled("figures/alex/ordermap_$(target_ei)_$(target_ni)_$(target_gi)_flow.png", "figures/alex/ordermap_$(target_ei)_$(target_ni)_$(target_gi)_OD.png", lhdataisbfilename, 2000, 1900, loclist, tripson, ordersbetween, group1tripson, group2tripson)
writeorderfile(outputfilename, sampleddata)

