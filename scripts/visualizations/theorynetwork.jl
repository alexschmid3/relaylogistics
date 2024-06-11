
using Luxor, Colors

tempHueList = collect(Colors.color_names)

#Exclude colors that are too light
hueList = []
for item in tempHueList
	if item[2][1] + item[2][2] + item[2][3] <= 615
		push!(hueList, item) 
	end
end

function theorynetwork(drawingname, y, x_size, y_size)

	#Find coordinates for each time-space node
	nodelist, nodelookup = [], Dict()
	x_size_trimmed, y_size_trimmed = x_size*0.95, y_size*0.95
	x_size_orig, y_size_orig = m*w, 2*h
    k_x, shift_x = x_size_trimmed / x_size_orig, -1/2 * x_size_trimmed
	k_y, shift_y = y_size_trimmed / y_size_orig, -1/2 * y_size_trimmed

	for i in pitstops
		ycoord = -1 * (coordinates[i,2] * k_y + shift_y)
		xcoord = coordinates[i,1] * k_x + shift_x
		push!(nodelist, (xcoord, ycoord))
		nodelookup[i] = Point((xcoord, ycoord))
	end

	#Create actual points as a Luxor object
	nodePoints = Point.(nodelist)

	#---------------------------------------------------------------------------------------#

    arccounter = Dict()
    allarcs = []
	for (i,j) in corridors, t in 1:T
        arccounter[i,j,t,Tmod(t+C)] = 0
        arccounter[j,i,t,Tmod(t+C)] = 0
        push!(allarcs, (i,j,t,Tmod(t+C)))
        push!(allarcs, (j,i,t,Tmod(t+C)))
    end
	for (j,i) in corridors, t in 1:T
        arccounter[i,j,t,Tmod(t+C)] = 0
        arccounter[j,i,t,Tmod(t+C)] = 0
        push!(allarcs, (i,j,t,Tmod(t+C)))
        push!(allarcs, (j,i,t,Tmod(t+C)))
    end
    #=for i in W, j in E, t in 1:T
        arccounter[i,j,t,Tmod(t+C)] = 0
        arccounter[j,i,t,Tmod(t+C)] = 0
        push!(allarcs, (i,j,t,Tmod(t+C)))
        push!(allarcs, (j,i,t,Tmod(t+C)))
    end
    for i in W, j in setdiff(W,i), t in 1:T
        arccounter[i,j,t,t] = 0
        push!(allarcs, (i,j,t,t))
    end
    for i in E, j in setdiff(E,i), t in 1:T
        arccounter[i,j,t,t] = 0
        push!(allarcs, (i,j,t,t))
    end=#

    for jindex in journeys
        #if value(y[jindex]) > 1e-4
            for (i,j,t1,t2) in journeyarclookup[jindex]
                arccounter[i,j,t1,t2] += value(y[jindex])
				#println("$i, $j --> ", value(y[jindex]))
            end
        #end
    end

    arcinfo = []
    maxthickness, maxflow = 20, maximum(values(arccounter))
    for (i,j,t1,t2) in allarcs
        if arccounter[i,j,t1,t2] > 1e-4
            startPoint = nodelookup[i]
            endPoint = nodelookup[j]
            arcColor =  (0,0,0) #((1 - demand[i,j,t1] / arccounter[i,j,t1,t2]) * 255, 0, 0)
            arcDash = "solid"
            arcThickness = ceil(maxthickness * arccounter[i,j,t1,t2]/maxflow)
            arcLabel = ""
            push!(arcinfo, (startPoint, endPoint, arcColor, arcDash, arcThickness, arcLabel))
        end
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

        #Rotate the arrowhead appropriately
		if i[1][1] >= i[2][1]
			local p = ngon(arrowhead, min(10, i[5]*2), 3, theta - pi , vertices=true)
		else
			local p = ngon(arrowhead, min(10, i[5]*2), 3, theta , vertices=true)
		end
		poly(p, :fill,  close=true)

	end

	#Draw node points
	setcolor("black")
	circle.(nodePoints, 9, :fill)
    
	#Set font size for labels
	fontsize(30)

	#Add location labels
	for l in pitstops
		coord = nodePoints[l]
		Luxor.text("$l", coord + Point(0,-30), halign=:center, valign = :middle)
	end
	#=
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
	end=#

	finish()
	preview()

end