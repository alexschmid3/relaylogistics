
using Luxor, Colors

tempHueList = collect(Colors.color_names)

#Exclude colors that are too light
hueList = []
for item in tempHueList
	if item[2][1] + item[2][2] + item[2][3] <= 615
		push!(hueList, item) 
	end
end

#drawingname, arclistlist, colorlist, thicknesslist, fractlist, x_size, y_size = string("outputs/viz/order", i,".png"), arclistlist, colorlist, thicknesslist, fractlist, 2000, 1200

function timespacenetwork(drawingname, arclistlist, colorlist, thicknesslist, dashlist, fractlist, x_size, y_size)

	#Find coordinates for each time-space node
	nodelist = []
	x_size_trimmed, y_size_trimmed = x_size*0.9, y_size*0.9
	k1 = x_size_trimmed/(horizon/tstep + 2) 
	k2 = y_size_trimmed/(numlocs + 2)
	for i in 1:numnodes
		ycoord = nodesLookup[i][1]
		xcoord = (nodesLookup[i][2]/tstep)+1

		#Scaling to image size
		tup = (-x_size_trimmed/2 + xcoord*k1, -y_size_trimmed/2 + ycoord*k2)   
		
		push!(nodelist,tup)
	end

	#Create actual points as a Luxor object
	nodePoints = Point.(nodelist)

	#---------------------------------------------------------------------------------------#

	#Arcs for visualization
	#Duplicate for multiple input arc lists with different colors/thickness/dash if you're trying to show m
	arcinfo = []
    for j in 1:length(arclistlist), a in intersect(1:numarcs,arclistlist[j])
        startPoint = nodePoints[arcLookup[a][1]]
        endPoint = nodePoints[arcLookup[a][2]]
        
        #Set arc attributes
        arcColor = colorlist[j] # (0,0,255) #RGB tuple 
        arcDash = dashlist[j] #"solid" #"solid", "dashed"			
        arcThickness = thicknesslist[j]
		#if fractlist[j] != []
		#	arcLabel = fractlist[j][a]
		#else
		arcLabel = ""
		#end
        
        #Add to arcinfo list to be used in the drawing 
        push!(arcinfo, (startPoint, endPoint, arcColor, arcDash, arcThickness, arcLabel))
    end

	#Legend arcs
	#legendarcs = []
	#push!(legendarcs, (Point(200,800), Point(200,800)+Point(120,0), colorlist[1], dashlist[1], thicknesslist[1]*2, ""))
	#push!(legendarcs, (Point(-700,800), Point(-700,800)+Point(120,0), colorlist[2], dashlist[2], thicknesslist[2]*2, ""))

	#-------------------------------------------------------------------------#

	#Initiailize drawing
	Drawing(x_size, y_size, drawingname)
	origin()
	background("white")

	#Add in gray boxes for night shifts
	#sethue(convert(Colors.HSV, Colors.RGB(220/255, 220/255, 220/255)))
	#rect(nodePoints[nodes[1,12]]+Point(-3,-6), 246, 1490, :fill)
	#rect(nodePoints[nodes[1,36]]+Point(-3,-6), 246, 1490, :fill)
	#rect(nodePoints[nodes[1,60]]+Point(-3,-6), 246, 1490, :fill)
	#rect(nodePoints[nodes[1,84]]+Point(-3,-6), 246, 1490, :fill)
	#fontsize(40)
	#setcolor("black")
	#for i in 1:4
	#	Luxor.text("Night $i", nodePoints[nodes[1,12+24*(i-1)]]+Point(120,1453), halign=:center, valign = :middle)
	#end

	#Draw arcs
	for i in arcinfo #union(arcinfo,legendarcs)
		
		#Set arc attributes from the arcinfo
		r_val, g_val, b_val = i[3][1]/255, i[3][2]/255, i[3][3]/255
		setcolor(convert(Colors.HSV, Colors.RGB(r_val, g_val, b_val)))  #You can also use setcolor("colorname")
		setdash(i[4])
		setline(i[5])

		#Draw the line from the start node to end node
		line(i[1], i[2] , :stroke)
		
		#Figure out the angle of the arrow head
		theta = atan((i[2][2] - i[1][2])/(i[2][1] - i[1][1]))
		dist = distance(i[1], i[2])
		if i in legendarcs
			arrowhead = (1-1/dist)*i[2] + (1/dist)*i[1] #1 pixels from the end node
			local p = ngon(arrowhead, max(10, i[5]*2), 3, theta, vertices=true)
		else
			arrowhead = (1-12/dist)*i[2] + (12/dist)*i[1] #12 pixels from the end node
			local p = ngon(arrowhead, max(10, i[5]), 3, theta, vertices=true)
		end
		
		#Draw the arrow head
		poly(p, :fill,  close=true)

	end

	#Draw node points
	setcolor("black")
	circle.(nodePoints, 9, :fill)

	#Set font size for labels
	fontsize(40)

	#Add location labels
	for l in 1:numlocs
		coord = nodePoints[nodes[(l,0.0)]]
		label("Loc $l       ", :W , coord)
	end

	#Add time labels
	for t in 0:tstep*2:horizon
		coord = nodePoints[nodes[(1,t)]] + Point(0,-30)
		label("t = $t", :N , coord)
	end

	#Add arc labels
	for i in arcinfo
		if i[6] != ""
			fontsize(20)
			coord = Point((i[1][1] + i[2][1])/2, (i[1][2] + i[2][2])/2)
			sethue("black")
			rect(coord+Point(-30,-25), 60, 25, :fill)
			sethue("white")
			label(i[6], :N , coord)
		end
	end

	#Legend labels
	#fontsize(60)
	#Luxor.text("Standard journey", Point(-500,800), halign=:left, valign = :middle)
	#Luxor.text("Extended journey", Point(400,800), halign=:left, valign = :middle)

	finish()
	preview()

end

#---------------------------------------------------------------------------------------------------------#

function timespacenetwork_providetsn(tsn, drawingname, arclistlist, colorlist, thicknesslist, dashlist, fractlist, x_size, y_size)

	#Find coordinates for each time-space node
	nodelist = []
	x_size_trimmed, y_size_trimmed = x_size*0.9, y_size*0.9
	k1 = x_size_trimmed/(tsn.horizon/tsn.tstep + 2) 
	k2 = y_size_trimmed/(tsn.numlocs + 2)
	for i in 1:tsn.numnodes
		ycoord = tsn.nodesLookup[i][1]
		xcoord = (tsn.nodesLookup[i][2]/tsn.tstep)+1

		#Scaling to image size
		tup = (-x_size_trimmed/2 + xcoord*k1, -y_size_trimmed/2 + ycoord*k2)   
		
		push!(nodelist,tup)
	end

	#Create actual points as a Luxor object
	nodePoints = Point.(nodelist)

	#---------------------------------------------------------------------------------------#

	#Arcs for visualization
	#Duplicate for multiple input arc lists with different colors/thickness/dash if you're trying to show m
	arcinfo = []
    for j in 1:length(arclistlist), a in intersect(1:tsn.numarcs,arclistlist[j])
        startPoint = nodePoints[tsn.arcLookup[a][1]]
        endPoint = nodePoints[tsn.arcLookup[a][2]]
        
        #Set arc attributes
        arcColor = colorlist[j] # (0,0,255) #RGB tuple 
        arcDash = dashlist[j] #"solid" #"solid", "dashed"			
        arcThickness = thicknesslist[j]
		#if fractlist[j] != []
		#	arcLabel = fractlist[j][a]
		#else
		arcLabel = ""
		#end
        
        #Add to arcinfo list to be used in the drawing 
        push!(arcinfo, (startPoint, endPoint, arcColor, arcDash, arcThickness, arcLabel))
    end

	#-------------------------------------------------------------------------#

	#Initiailize drawing
	Drawing(x_size, y_size, drawingname)
	origin()
	background("white")

	#Draw arcs
	for i in arcinfo
		
		#Set arc attributes from the arcinfo
		r_val, g_val, b_val = i[3][1]/255, i[3][2]/255, i[3][3]/255
		setcolor(convert(Colors.HSV, Colors.RGB(r_val, g_val, b_val)))  #You can also use setcolor("colorname")
		setdash(i[4])
		setline(i[5])

		#Draw the line from the start node to end node
		line(i[1], i[2] , :stroke)
		
		#Figure out the angle of the arrow head
		theta = atan((i[2][2] - i[1][2])/(i[2][1] - i[1][1]))
		dist = distance(i[1], i[2])
		arrowhead = (1-12/dist)*i[2] + (12/dist)*i[1] #8 pixels from the end node
		
		#Draw the arrow head
		local p = ngon(arrowhead, max(10, i[5]), 3, theta, vertices=true)
		poly(p, :fill,  close=true)

	end

	#Draw node points
	setcolor("black")
	circle.(nodePoints, 9, :fill)

	#Set font size for labels
	fontsize(30)

	#Add location labels
	for l in 1:tsn.numlocs
		coord = nodePoints[nodes[(l,0.0)]]
		label("Loc $l       ", :W , coord)
	end

	#Add time labels
	for t in 0:tsn.tstep*2:tsn.horizon
		coord = nodePoints[tsn.nodes[(1,t)]] + Point(0,-30)
		label("t = $t", :N , coord)
	end

	#Add arc labels
	for i in arcinfo
		if i[6] != ""
			fontsize(20)
			coord = Point((i[1][1] + i[2][1])/2, (i[1][2] + i[2][2])/2)
			sethue("black")
			rect(coord+Point(-30,-25), 60, 25, :fill)
			sethue("white")
			label(i[6], :N , coord)
		end
	end

	#Add horizon line
	line_x = -x_size_trimmed/2 + ((horizon/tsn.tstep)+1)*k1 
	line_y1 = -y_size_trimmed/2 + (tsn.numlocs+1)*k2
	line_y2 = -y_size_trimmed/2
	setcolor("black")
	setdash("dashed")
	setline(3)
	line(Point((line_x, line_y1)), Point((line_x, line_y2)), :stroke)

	finish()
	preview()

end


