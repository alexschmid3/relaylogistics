
using CSV, DataFrames

include("../relaylogistics/scripts/instancegeneration/readrivigodata.jl")

numdrivers = 3316
tstep, horizon, tstep_ol = 6, 162, 6 
maxlocs = 66
shiftlength = 12

hubdistancesfilename = "../relaylogistics/data/hubdistances.csv"
traveltimesfilename = "../relaylogistics/data/traveltimes_outliers.csv"
hubdataisbfilename = "../relaylogistics/data/hub_data_isb_connect.csv"

roundup_flag = 1	
excludeoutliers_flag = 1		 				
googlemapstraveltimes_flag = 1
includesymmetricarcs_flag = 1						
traveltimefordelay_flag = 2 						
ensureconnectivity_flag = 1

hubCoords, hubsLookup, hubsReverseLookup, hubsTravelTimeIndex, numlocs = readlocations(hubdataisbfilename, maxlocs)
prearcs, arcLength, arcLength_raw = readarcs(traveltimesfilename, hubdistancesfilename, tstep, numlocs, hubsTravelTimeIndex, roundup_flag, excludeoutliers_flag, hubsReverseLookup, googlemapstraveltimes_flag)
distbetweenlocs, shortestpatharclists = cacheShortestDistance(numlocs, prearcs)

#Read x_sol
pastordersegments = []
for t in 0:tstep_ol:horizon
    data = CSV.read(string("outputs/online/ex4_rr_rundate2022-03-23/solution/x_soln",t,".csv"), DataFrame)
    for row in 1:size(data)[1]
        if data[row,9] > 1e-4
            l1, l2 = data[row,5], data[row,6]
            t1, t2 = data[row,7], data[row,8]
            if (l1!=l2) & (t1 < t+tstep_ol-1e-4)
                push!(pastordersegments, (l1,l2,t1,t2))
            end
        end
    end
end

#Read z_sol
pastdriversegments = Dict()
for d in 1:3318
    pastdriversegments[d] = []
end
for t in 0:tstep_ol:horizon
    data = CSV.read(string("outputs/online/ex4_rr_rundate2022-03-23/solution/z_soln",t,".csv"), DataFrame)
    for row in 1:size(data)[1]
        if data[row,9] > 1e-4
            l1, l2 = data[row,5], data[row,6]
            t1, t2 = data[row,7], data[row,8]
            d = data[row,3]
            if (l1!=l2) & (t1 < t+tstep_ol-1e-4)
                push!(pastdriversegments[d], (l1,l2,t1,t2))
            end
        end
    end
end

#Find the empty versus full
em, om = Dict(), Dict()
for d in 1:numdrivers
    empty, full = 0, 0
    for (l1,l2,t1,t2) in pastdriversegments[d]
        if (l1,l2,t1,t2) in pastordersegments
            full += distbetweenlocs[l1,l2]
            deleteat!(pastordersegments, findfirst(x->x==(l1,l2,t1,t2), pastordersegments))
        else
            empty += distbetweenlocs[l1,l2]
        end
    end
    em[d] = empty
    om[d] = full
end

df = DataFrame(
    method = ["rr" for d in 1:numdrivers],
    example = [4 for d in 1:numdrivers],
    driverid = [d for d in 1:numdrivers],
    emptymiles = [em[d] for d in 1:numdrivers],
    ordermiles = [om[d] for d in 1:numdrivers],
    percentempty = [em[d]+om[d] == 0 ? "" : em[d]/(em[d]+om[d]) for d in 1:numdrivers]
)

CSV.write("../relaylogistics/figures/onlinehistograms/driverempty_rr.csv", df)