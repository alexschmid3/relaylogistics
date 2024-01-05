
using Luxor, Colors

tempHueList = collect(Colors.color_names)

#Exclude colors that are too light
hueList = []
for item in tempHueList
	if item[2][1] + item[2][2] + item[2][3] <= 615
		push!(hueList, item) 
	end
end

function timespacenetwork(drawingname, arclistlist, colorlist, thicknesslist, x_size, y_size)

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
        arcDash = "solid" #"solid", "dashed"			
        arcThickness = thicknesslist[j]
        
        #Add to arcinfo list to be used in the drawing 
        push!(arcinfo, (startPoint, endPoint, arcColor, arcDash, arcThickness))
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
		local p = ngon(arrowhead, 10, 3, theta, vertices=true)
		poly(p, :fill,  close=true)
	end

	#Draw node points
	setcolor("black")
	circle.(nodePoints, 7, :fill)

	#Set font size for labels
	fontsize(20)

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

	finish()
	preview()

end