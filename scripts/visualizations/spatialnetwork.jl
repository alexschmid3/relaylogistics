
using CSV, Luxor, Colors, Random, DataFrames, Dates, StatsBase

#include("scripts/instancegeneration/readrivigodata.jl")

#--------------------------------------------------------------------------------------------------#

thickest, thinnest = 20, 1
pixelshift = 22
maxlocs = 66
hubCoords, hubsLookup, hubsReverseLookup, hubsTravelTimeIndex, numlocs = readlocations(hubdataisbfilename, maxlocs)

#--------------------------------------------------------------------------------------------------#

function getrivigotriphistory(lhdataisbfilename)

	data_agg = CSV.read(lhdataisbfilename, DataFrame)

	tripson, origincount, destinationcount = Dict(), Dict(), Dict()
    for i in 1:numlocs, j in 1:numlocs
        tripson[i,j] = 0
    end
	for i in 1:numlocs
        origincount[i] = 0
		destinationcount[i] = 0
    end

    hubsList = collect(values(hubsLookup))
	for i in 1:size(data_agg)[1]
		#if data_agg[i,1] in [7882,8502,8388,7221,2723,5581,2606,8,737,1219]
			orig, dest = data_agg[!,26][i], data_agg[!,27][i]
			psseq_raw = data_agg[i,8]
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

			if (orig != dest) & (1 <= orig <= numlocs) & (1 <= dest <= numlocs) & (intermedlocs_flag == 0) #& (orderwindowstart <= start_avail_ts <= orderwindowend) 
				for i in 1:length(stopsequence)-1
					tripson[stopsequence[i], stopsequence[i+1]] += 1
				end
				origincount[orig] += 1
				destinationcount[dest] += 1
			end
		#end
	end

	return tripson, origincount, destinationcount

end

#--------------------------------------------------------------------------------------------------#

function spatialnetwork(drawingname, lhdataisbfilename, xdim, ydim)

    tripson, origincount, destinationcount = getrivigotriphistory(lhdataisbfilename)

    #--------------------------------------------------------#

    #Get correct scale
	maxlat, minlat = 0, 100
	for l in 1:numlocs
        maxlat = max(hubCoords[l,1], maxlat)
        minlat = min(hubCoords[l,1], minlat)
    end
    
    latmult = -(xdim-200) / (maxlat - minlat)	 
    latshift = -(xdim-200)/2 + (xdim-200) * maxlat / (maxlat - minlat)
	longmult = -1 * latmult * 24/29
	maxlongcoord, minlongcoord = -100000, 100000
	for l in 1:numlocs
        maxlongcoord = max(longmult*hubCoords[l,2], maxlongcoord)
        minlongcoord = min(longmult*hubCoords[l,2], minlongcoord)
    end
	longshift = -(maxlongcoord + minlongcoord)/2
	 
	#Format and transform latitude and longitude coordinates of each pit stop
	pointDict = Dict()
	listofpoints = []
	listofpoints_labels = []
	for l in 1:numlocs
		longitude, latitude = hubCoords[l,2], hubCoords[l,1]
		transformedcoords = (longmult*longitude+longshift, latmult*latitude+latshift)
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
        if tripson[i,j] >= 1
            startPoint = locationPoints[i]
            endPoint = locationPoints[j]
            thickness = round(thinnest + (tripson[i,j] - mintrips)/(maxtrips - mintrips) * (thickest - thinnest) )
            push!(arcList, (startPoint, endPoint, (0,0,0), thickness, "solid"))
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
		line(i[1], i[2] , :stroke)
		
		#Calculate the rotation and placement of the arrowhead
		theta = atan((i[2][2] - i[1][2])/(i[2][1] - i[1][1]))
		dist = distance(i[1], i[2])
		arrowhead = (1-pixelshift/dist)*i[2] + (pixelshift/dist)*i[1] #center of arrowhead positioned 8 pixels from the end node

		#Rotate the arrowhead appropriately
		if i[1][1] >= i[2][1]
			local p = ngon(arrowhead, min(pixelshift, i[4]*2), 3, theta - pi , vertices=true)
		else
			local p = ngon(arrowhead, min(pixelshift, i[4]*2), 3, theta , vertices=true)
		end

		#Draw the arrowhead
		poly(p, :fill,  close=true)
	end

	#Draw the pit stop nodes
	setcolor("red")
	circle.(locationPoints, 16, :fill)
    setcolor("black")
    setline(3)
    circle.(locationPoints, 16, :stroke)

	#Add pit stop labels
	fontsize(22)
	setcolor("white")
	for item in listofpoints_labels
 		#label(item[2], :0, Point(item[1]))
		Luxor.text(item[2],  Point(item[1]), halign=:center,   valign = :middle)
	end
	setcolor("black")

    #Legend box
    setline(4)
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
    end

	#--------------------------------------------------------#

	finish()
	preview()

end


spatialnetwork("figures/orderhistorymap_numbers.png", lhdataisbfilename, 2000, 1900)

#spatialnetwork("figures/truck18path.png", lhdataisbfilename, 2000, 1900)

#--------------------------------------------------------#

