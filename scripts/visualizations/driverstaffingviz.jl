
using Luxor, CSV, DataFrames

include("scripts/instancegeneration/readrivigodata.jl")
include("scripts/visualizations/spatialnetwork.jl")

hubdataisbfilename = "data/hub_data_isb_connect.csv"
lhdataisbfilename = "data/lh_data_isb_connect_clean.csv"
maxlocs = 66
pixelshift = 25

#=
#Generate input file
tripson, origincount, destinationcount = getrivigotriphistory(lhdataisbfilename)

allcorridors = [c for c in collect(keys(tripson)) if tripson[c] > 1e-4]
sortedcorridors = reverse(sort(allcorridors, by=x->tripson[x]))

allpitstoptraffic = zeros(numlocs)
for (l1,l2) in allcorridors
	allpitstoptraffic[l1] += tripson[l1,l2]
	allpitstoptraffic[l2] += tripson[l1,l2]
end
sortedpitstops = reverse(sort(1:numlocs, by=x->allpitstoptraffic[x]))

allpitstopdemand = zeros(numlocs)

for i in 1:numlocs
	println("$i,",origincount[i],",",destinationcount[i],",",allpitstoptraffic[i])
end

CSV.write("data/driversensitivity/flowlocations.csv", DataFrame(
    loc=[i for i in 1:numlocs],
    orig=[origincount[i] for i in 1:numlocs],
    dest=[destinationcount[i] for i in 1:numlocs],
    total=[allpitstoptraffic[i] for i in 1:numlocs])
=#

data = CSV.read("data/driversensitivity/flowlocations.csv", DataFrame)

hubCoords, hubsLookup, hubsReverseLookup, hubsTravelTimeIndex, numlocs = readlocations(hubdataisbfilename, maxlocs)

#--------------------------------------------------------------------------------------------------#

function driverstaffingnetwork(drawingname, data, lhdataisbfilename, xdim, ydim)

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

    #Sort locs into categories
    loccatgeory = Dict()
    for row in 1:size(data,1)
        loccatgeory[convert(Int, data[row,1])] = convert(String, data[row,7])
    end

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
    for i in 1:maxlocs
	    if loccatgeory[i] == "corridor"
            r_val, g_val, b_val = 255, 0, 0
        elseif loccatgeory[i] == "hybrid"
            r_val, g_val, b_val = 153, 29, 29
        elseif loccatgeory[i] == "demand"
            r_val, g_val, b_val = 0, 0, 0
        end
        setcolor(convert(Colors.HSV, Colors.RGB(r_val/255, g_val/255, b_val/255)))
        circle(locationPoints[i], 20, :fill)
        setcolor("black")
        setline(3)
        circle(locationPoints[i], 20, :stroke)
    end

	#Add pit stop labels
	fontsize(22)
	setcolor("white")
	for item in listofpoints_labels
 		#label(item[2], :0, Point(item[1]))
		Luxor.text(item[2],  Point(item[1]), halign=:center,   valign = :middle)
	end
	setcolor("black")

	#--------------------------------------------------------#

    #Legend box
    setline(4)
    legendstartx = 0.5*xdim - 0.43*xdim
    legendstarty = 0.5*ydim - 0.3*ydim
    rect(legendstartx, legendstarty, 0.4*xdim, 0.25*ydim, :stroke)

    #Arcs for the legend
    fontsize(70)
    numlegendarcs = 3
    legendlabels = ["High flow location", "Mid flow location", "Low flow location"]
    for legendarc in 1:numlegendarcs
        centerPoint = Point(legendstartx + 0.03*xdim, legendstarty + (legendarc-0.5)/numlegendarcs * 0.25*ydim)
                
        if legendarc == 1
            r_val, g_val, b_val = 255, 0, 0
        elseif legendarc == 2
            r_val, g_val, b_val = 153, 29, 29
        elseif legendarc == 3
            r_val, g_val, b_val = 0, 0, 0
        end
        setcolor(convert(Colors.HSV, Colors.RGB(r_val/255, g_val/255, b_val/255)))
        circle(centerPoint, 20, :fill)
        setcolor("black")
        setline(3)
        circle(centerPoint, 20, :stroke)
        
        #Add the label
        label(legendlabels[legendarc], :E , centerPoint + Point(xdim/40, 0))
    end

    #--------------------------------------------------------#    

	finish()
	preview()

end

driverstaffingnetwork("figures/driverstaffingviz.png", data, lhdataisbfilename, 2000, 1900)
