
using CSV, Luxor, Colors, Random, DataFrames, Dates, StatsBase, JuMP, Gurobi

include("scripts/instancegeneration/readrivigodata.jl")

#--------------------------------------------------------------------------------------------------#

runtype = "all" # "subntwk" "all"
numlocs = 66
basefilename = "trafficfigs"

#--------------------------------------------------------------------------------------------------#

#Figure design parameters
thickest, thinnest = 40, 2
pixelshift = 40
pitstopsize = 35
NSarcshift = 21
pitstoplabels = false
subnetworklegend = true

#Paper colors
fullflowcolor = (50,50,50)
subNcolor = (100,143,255)
subScolor = (254,97,0)

#Slides colors
#fullflowcolor = (28,41,173)
#subNcolor = (51,43,136) 
#subScolor = (204,102,119) 

fullflowcolor = fullflowcolor ./ 255
subNcolor = subNcolor ./ 255
subScolor = subScolor ./ 255

#Full traffic figure
if runtype == "all"
    fulldata, tripson, origincount, destinationcount, ordersbetween, group1tripson, group2tripson = pullrivigosample(lhdataisbfilename, loclist)
    drawingname, xdim, ydim = "fulltraffic.png", 3000, 3000
    spatialnetwork_fullrivigo(drawingname, lhdataisbfilename, xdim, ydim, tripson)
end

#Subnetwork figure
if runtype == "subntwk"
    subdata, tripson, origincount, destinationcount, ordersbetween, group1tripson, group2tripson = pullrivigosample(lhdataisbfilename, loclist)
    drawingname, xdim, ydim = "corridor.png", 3000, 3000
    spatialnetwork_subntwk(drawingname, xdim, ydim, tripson, group1tripson, group2tripson)
end

#--------------------------------------------------------------------------------------------------#

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

#--------------------------------------------------------------------------------------------------#

if runtype == "subntwk"
    ORIENTATION = "NS"
    outputfilename = "data/orderbalance/"*basefilename*".csv"
    flowpngname = "figures/alex/"*basefilename*"_flow.png"
    ODpngname = "figures/alex/"*basefilename*"_OD.png"
elseif runtype == "all"
    ORIENTATION = "all"
    outputfilename = "none.csv"
    flowpngname = "figures/alex/"*basefilename*"_flow.png"
    ODpngname = "figures/alex/"*basefilename*"_OD.png"
end

if ORIENTATION == "EW"
    #East/West
    #group1 = [63,65,66,58,52,41,42,46,51,59,49,48,57,34,38, 58] #East
    #group2 = [50, 56, 55, 53, 28, 27, 31, 36, 40, 39, 35, 32, 29, 26, 45] #West
    group1 = [63, 65, 66, 58, 52, 41, 42, 46, 51]
    group2 = [56, 55, 53, 28, 27, 31, 36, 40, 39, 35, 32, 29, 26, 45]
    loclist = union(group1, group2)
    region1 = [56, 55, 53, 63, 65, 66, 51, 46, 45, 58] #North
    region2 = [52, 41, 42, 31, 36, 40, 39, 28, 27, 35, 32, 29, 26] #South
    DIVIDEINTOREGIONS = true

elseif ORIENTATION == "NS"
    #North/South
    group1 = [64, 61, 60, 62, 56, 55, 53, 54, 45] #n=9
    group2 = [25, 19, 20, 17, 16, 13, 14, 8, 7, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12] #n=18
    loclist = [64, 61, 60, 62, 56, 55, 53, 54, 45, 37, 30, 25, 19, 20, 17, 16, 13, 14, 8, 7, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12]
    DIVIDEINTOREGIONS = false

elseif ORIENTATION == "NS2"
    #North/South
    group1 = [60, 45, 54, 36, 47, 44, 53, 55, 56, 43] #North
    group2 = [26, 22, 18, 15, 12, 21, 1, 10, 7, 9, 5] #South
    region1 = [26, 22, 18, 15, 12, 21, 54, 36, 47, 44, 43] #West
    region2 = [1, 7, 9, 10, 45, 60, 5, 53, 55, 56] #East
    loclist = union(group1, group2)
    g1r1, g2r1, g1r2, g2r2 = intersect(group1, region1), intersect(group2, region1), intersect(group1, region2), intersect(group2, region2)

    DIVIDEINTOREGIONS = true

elseif ORIENTATION == "all"
    #All
    group1 = 1:numlocs
    group2 = 1:numlocs
    loclist = 1:numlocs
    DIVIDEINTOREGIONS = false

end

#--------------------------------------------------------------------------------------------------#

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
        if ORIENTATION == "all" 
            if (orig != dest) & (intermedlocs_flag == 0) 
                push!(rivigodata, [data_agg[!, 1][i], orderdayofweek, ordertimeofday, orig, dest, stopsequence])
            end
        elseif !DIVIDEINTOREGIONS
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
    #println("Edge imbalance = ", edgeimbalance / totalflow)
    #println("Node imbalance = ", nodeimbalance / totalflow / 2)
    #println("Global imbalance = ", globalimbalance / totalflow)
    #println("Total flow = ", totalflow)

    if DIVIDEINTOREGIONS
        region1imbalance = abs(sum(sum(ordersbetween[i, j] for j in intersect(region1,group2)) for i in intersect(region1,group1)) - sum(sum(ordersbetween[j, i] for j in intersect(region1,group2)) for i in intersect(region1,group1)))
        region1flow = sum(sum(ordersbetween[i, j] for j in intersect(region1,group2)) for i in intersect(region1,group1)) + sum(sum(ordersbetween[j, i] for j in intersect(region1,group2)) for i in intersect(region1,group1))
        region2flow = sum(sum(ordersbetween[i, j] for j in intersect(region2,group2)) for i in intersect(region2,group1)) + sum(sum(ordersbetween[j, i] for j in intersect(region2,group2)) for i in intersect(region2,group1))
        region2imbalance = abs(sum(sum(ordersbetween[i, j] for j in intersect(region2,group2)) for i in intersect(region2,group1)) - sum(sum(ordersbetween[j, i] for j in intersect(region2,group2)) for i in intersect(region2,group1)))

        println("Flow 1 = $region1flow")
        println("Flow 2 = $region2flow")
        
        return edgeimbalance / totalflow, nodeimbalance / totalflow / 2, globalimbalance / totalflow, region1imbalance / region1flow, region2imbalance / region2flow

    else

        return edgeimbalance / totalflow, nodeimbalance / totalflow / 2, globalimbalance / totalflow, 0, 0

    end
    
end

#--------------------------------------------------------------------------------------------------#

#Pulls one week worth of Rivigo data, randomly sampled according to parameters: target_ei, target_ni, target_gi
function pullrivigosample(lhdataisbfilename, loclist)

    rivigodata = readrivigoorders(lhdataisbfilename, loclist)
    sampleddata = rivigodata

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
            elseif orig in group2
                group2tripson[stopsequence[i], stopsequence[i+1]] += 1
            else
                println("huh")
                println( orig, dest, stopsequence)
            end
        end
        origincount[orig] += 1
        destinationcount[dest] += 1
        push!(numpitstops, length(stopsequence))
        push!(tripdists, tripdist)
    end

    return sampleddata, tripson, origincount, destinationcount, ordersbetween, group1tripson, group2tripson

end

#--------------------------------------------------------------------------------------------------#

function spatialnetwork_subntwk(drawingname, xdim, ydim, tripson, group1tripson, group2tripson)

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

    group1tripson[12, 16] = 5
    group1tripson[16, 12] = 5
    group1tripson[45, 54] = 5
    group1tripson[54, 45] = 5
    group2tripson[12, 16] = 5
    group2tripson[16, 12] = 5
    group2tripson[45, 54] = 5
    group2tripson[54, 45] = 5

    #Calculate thickness of each arc
    mintrips, maxtrips = 1, max(maximum(values(group1tripson)), maximum(values(group2tripson)))
    arcList = []
    for i in loclist, j in setdiff(loclist, i)
        if group1tripson[i, j] >= 1
            startPoint = locationPoints[i]
            endPoint = locationPoints[j]
            thickness = round(thinnest + (group1tripson[i, j] - mintrips) / (maxtrips - mintrips) * (thickest - thinnest))
            if ORIENTATION[1:2] in ["NS", "al"]
                push!(arcList, (startPoint + Point(NSarcshift, 0), endPoint + Point(NSarcshift, 0), subScolor, thickness, "solid"))
            elseif ORIENTATION[1:2] == "EW"
                push!(arcList, (startPoint + Point(0, NSarcshift), endPoint + Point(0, NSarcshift), subScolor, thickness, "solid"))
            end
        end
        if group2tripson[i, j] >= 1
            startPoint = locationPoints[i]
            endPoint = locationPoints[j]
            thickness = round(thinnest + (group2tripson[i, j] - mintrips) / (maxtrips - mintrips) * (thickest - thinnest))
            if ORIENTATION[1:2] in ["NS", "al"]
                push!(arcList, (startPoint - Point(NSarcshift, 0), endPoint - Point(NSarcshift, 0), subNcolor, thickness, "solid"))
            elseif ORIENTATION[1:2] == "EW"
                push!(arcList, (startPoint - Point(0, NSarcshift), endPoint - Point(0, NSarcshift), subNcolor, thickness, "solid"))
            end
        end
    end

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
        if i[1][1] > i[2][1]
            local p = ngon(arrowhead, min(pixelshift, i[4] * 2), 3, theta - pi, vertices=true)
        else
            local p = ngon(arrowhead, min(pixelshift, i[4] * 2), 3, theta, vertices=true)
        end

        #Draw the arrowhead
        poly(p, :fill, close=true)
    end

    #Draw the pit stop nodes
    for pt in 1:length(locationPoints)
        if pt in union(group1, group2, [37,30])
            setcolor("black")
            circle(locationPoints[pt], pitstopsize, :fill)
            setcolor("black")
            setline(3)
            circle(locationPoints[pt], pitstopsize, :stroke)
        else
            setcolor((0.8,0.8,0.8))
            circle(locationPoints[pt], pitstopsize, :fill)
            setcolor((0.5,0.5,0.5))
            setline(3)
            circle(locationPoints[pt], pitstopsize, :stroke)
        end 
    end

    #Add pit stop labels
    fontsize(22)
    setcolor("white")
    if pitstoplabels
        for item in listofpoints_labels
            #label(item[2], :0, Point(item[1]))
            Luxor.text(item[2], Point(item[1]), halign=:center, valign=:middle)
        end
    end
    setcolor("black")

    #Legend
    if subnetworklegend
        #Box
        setline(4)
        legendstartx = 0.5*xdim - 0.43*xdim
        legendstarty = 0.5*ydim - 0.3*ydim
        rect(legendstartx, legendstarty, 0.4*xdim, 0.25*ydim, :stroke)

        #Labels
        fontsize(140)
        Luxor.text("Northbound", Point(legendstartx + 0.1*xdim, legendstarty + 0.25/3*ydim), halign=:left, valign=:middle)
        Luxor.text("Southbound", Point(legendstartx + 0.1*xdim, legendstarty + 0.5/3*ydim), halign=:left, valign=:middle)
    
        #Northbound
        startPoint_x = legendstartx + 0.05*xdim
        startPoint_y = legendstarty + 0.25/3*ydim + 0.025*ydim
        endPoint_x = legendstartx + 0.05*xdim
        endPoint_y = legendstarty + 0.25/3*ydim - 0.025*ydim
        setline(thickest)
        setdash("solid")
        setcolor(subNcolor)
        line(Point(startPoint_x, startPoint_y), Point(endPoint_x, endPoint_y), :stroke)
        theta = atan((endPoint_y - startPoint_y) / (endPoint_x - startPoint_x))
        arrowhead = Point(endPoint_x, endPoint_y)
        local p = ngon(arrowhead, thickest, 3, theta, vertices=true)
        poly(p, :fill, close=true)

        #Southbound
        startPoint_x = legendstartx + 0.05*xdim
        startPoint_y = legendstarty + 2*0.25/3*ydim - 0.025*ydim
        endPoint_x = legendstartx + 0.05*xdim
        endPoint_y = legendstarty + 2*0.25/3*ydim + 0.025*ydim
        setline(thickest)
        setdash("solid")
        setcolor(subScolor)
        line(Point(startPoint_x, startPoint_y), Point(endPoint_x, endPoint_y), :stroke)
        theta = atan((endPoint_y - startPoint_y) / (endPoint_x - startPoint_x))
        arrowhead = Point(endPoint_x, endPoint_y)
        local p = ngon(arrowhead, thickest, 3, theta, vertices=true)
        poly(p, :fill, close=true)
    end

    #--------------------------------------------------------#

    finish()
    preview()

end

#--------------------------------------------------------------------------------------------------#

function spatialnetwork_fullrivigo(drawingname, lhdataisbfilename, xdim, ydim, tripson)

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
    arcList = []
    totalorders = Dict()
    for i in 1:numlocs, j in 1:numlocs
        totalorders[i,j] = 0
    end
    data = CSV.read(lhdataisbfilename, DataFrame)
    
    tripson[12, 16] = 5
    tripson[16, 12] = 5
    tripson[45, 54] = 5
    tripson[54, 45] = 5

    mintrips, maxtrips = 1, maximum(values(tripson))
    arcList = []
    for (i,j) in keys(tripson)
        if tripson[i,j] >= 1
            startPoint = locationPoints[i]
            endPoint = locationPoints[j]
            thickness = round(thinnest + (tripson[i, j] - mintrips) / (maxtrips - mintrips) * (thickest - thinnest))
            push!(arcList, (startPoint, endPoint, fullflowcolor, thickness, "solid"))
        end
    end

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
        pixelshift_scaled = pixelshift
        arrowhead = (1 - pixelshift_scaled / dist) * i[2] + (pixelshift_scaled / dist) * i[1] #center of arrowhead positioned "pixelshift_scaled" pixels from the end node

        #Rotate the arrowhead appropriately
        if i[1][1] > i[2][1]
            local p = ngon(arrowhead, min(pixelshift_scaled, i[4] * 2), 3, theta - pi, vertices=true)
        else
            local p = ngon(arrowhead, min(pixelshift_scaled, i[4] * 2), 3, theta, vertices=true)
        end

        #Draw the arrowhead
        poly(p, :fill, close=true)
    end

    #Draw the pit stop nodes
    for pt in 1:length(locationPoints)
        setcolor("black")
        circle(locationPoints[pt], pitstopsize, :fill)
    end
    setcolor("black")
    setline(3)
    circle.(locationPoints, pitstopsize, :stroke)

    #Add pit stop labels
    fontsize(40)
    setcolor("white")
    if pitstoplabels
        for item in listofpoints_labels
            #label(item[2], :0, Point(item[1]))
            Luxor.text(item[2], Point(item[1]), halign=:center, valign=:middle)
        end
    end
    setcolor("black")

    #--------------------------------------------------------#

    finish()
    preview()

end

#Full traffic figure
if runtype == "all"
    fulldata, tripson, origincount, destinationcount, ordersbetween, group1tripson, group2tripson = pullrivigosample(lhdataisbfilename, loclist)
    drawingname, xdim, ydim = "fulltraffic.png", 3000, 3000
    spatialnetwork_fullrivigo(drawingname, lhdataisbfilename, xdim, ydim, tripson)
end

#Subnetwork figure
if runtype == "subntwk"
    subdata, tripson, origincount, destinationcount, ordersbetween, group1tripson, group2tripson = pullrivigosample(lhdataisbfilename, loclist)
    drawingname, xdim, ydim = "corridor.png", 3000, 3000
    spatialnetwork_subntwk(drawingname, xdim, ydim, tripson, group1tripson, group2tripson)
end


