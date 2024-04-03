
using Luxor, Colors, CSV, DataFrames

tempHueList = collect(Colors.color_names)

#Exclude colors that are too light
hueList = []
for item in tempHueList
	if (item[2][1] + item[2][2] + item[2][3] <= 515) & (max(abs(item[2][1]-item[2][2]), abs(item[2][1]-item[2][3]), abs(item[2][2]-item[2][3])) >= 100)
		push!(hueList, item) 
	end
end
	
#---------------------------------------------------------------------------------------#

function generateheatmap(inputfilename, outputfilename, coastbalance, stepsize, size_x, size_y)

	buffer_xy = 300

	#Find coordinates for each workstation and pod storage location
	boxwidth = (size_x - buffer_xy) / (length(0:stepsize:1)+1)
	boxheight = (size_y - buffer_xy) / (length(0:stepsize:1)+1)

	boxPoints = Dict()
	for ab in 0:stepsize:1, db in 0:stepsize:1
		newx = ab * length(0:stepsize:1)/(length(0:stepsize:1)+1) * (size_x - buffer_xy) - 0.5 * (size_x - buffer_xy)
		newy = (1-db) * length(0:stepsize:1)/(length(0:stepsize:1)+1) * (size_y - buffer_xy) - 0.5 * (size_y - buffer_xy)
		boxPoints[ab,db] = Point((newx, newy))
	end

	#-------------------------------------------------------------------------#

	#Read the input file
	data = CSV.read(inputfilename, DataFrame)
	mileslookup, pctlookup = Dict(), Dict()

	for ab in 0:stepsize:1, db in 0:stepsize:1
		mileslookup[ab,db] = "X"
	end

	allmiles = []
	for line in 1:size(data)[1]
		if data[line,3] == coastbalance
			ab, db = round(stepsize * round(data[line,1] / stepsize), digits=2), round(stepsize * round(data[line,2] / stepsize),digits=2)
			mileslookup[ab,db] = (data[line,5] - data[line,4])  
			pctlookup[ab,db] = 100 * (data[line,4] - data[line,5]) / data[line,5] 
			push!(allmiles, data[line,5] - data[line,4])
		end
	end	

	maxmiles, minmiles = maximum(allmiles), minimum(allmiles)

	#-------------------------------------------------------------------------#
	
	#Create a list of boxes and color types
	thickness = 0
	textfontsize = 10

	counter = 0

	locationsquares, textsquares = [], []
	for ab in 0:stepsize:1, db in 0:stepsize:1
		corner = boxPoints[ab,db]
		center = boxPoints[ab,db] + Point(boxwidth/2, boxheight/2) 
		if mileslookup[ab,db] == "X"
			boxcolor = (200,200,200)
			textcolor = (150,150,150)
			actualtext = "X"
			push!(locationsquares, (corner, boxcolor, boxwidth, boxheight, thickness))
			push!(textsquares, (center, textcolor, actualtext, 36))
		elseif mileslookup[ab,db] >= 0 
			println("$ab, $db")
			lamb = mileslookup[ab,db] / maxmiles
			boxcolor = (98 * lamb + 255 * (1-lamb), 151 * lamb + 255 * (1-lamb), 236 * lamb + 255 * (1-lamb))
			textcolor = ((98 * lamb + 255 * (1-lamb)) / 2, (151 * lamb + 255 * (1-lamb)) / 2, (236 * lamb + 255 * (1-lamb)) / 2)
			actualtext = string(convert(Int, round(pctlookup[ab,db],digits=0)), "%")
			push!(locationsquares, (corner, boxcolor, boxwidth, boxheight, thickness))
			push!(textsquares, (center, textcolor, actualtext, 20))
			counter +=1
		else
			println("$ab, $db")
			lamb = mileslookup[ab,db] / minmiles
			boxcolor = (247 * lamb + 255 * (1-lamb), 91 * lamb + 255 * (1-lamb), 95 * lamb + 255 * (1-lamb))
			textcolor = ((247 * lamb + 255 * (1-lamb)) / 2, (91 * lamb + 255 * (1-lamb)) / 2, (95 * lamb + 255 * (1-lamb)) / 2)
			actualtext = string("+",convert(Int,round(pctlookup[ab,db],digits=0)), "%")
			push!(locationsquares, (corner, boxcolor, boxwidth, boxheight, thickness))
			push!(textsquares, (center, textcolor, actualtext, 20))
			counter +=1
		end		
	end

	#-------------------------------------------------------------------------#

	Drawing(size_x, size_y, outputfilename)
	origin()
	background("white")

	#Heatmap boxes
	for box in locationsquares
		#setline(box[4])
		r_val, g_val, b_val = box[2]
		setcolor(convert(Colors.HSV, Colors.RGB(r_val/255, g_val/255, b_val/255)))
		Luxor.rect(box[1], box[3], box[4], :fill)	
	end

	#Mile labels
	for lbl in textsquares
		fontsize(lbl[4])
		r_val, g_val, b_val = lbl[2]
		setcolor(convert(Colors.HSV, Colors.RGB(r_val/255, g_val/255, b_val/255)))
		text(lbl[3], lbl[1], halign=:center,   valign = :middle)
		#label(lbl[3], :0, lbl[1])
	end

	#Chart Labels
	setcolor("black")
	setdash("solid")
	setline(2)

	#Box
	Luxor.rect(Point(-(size_x - buffer_xy) / 2, - (size_y - buffer_xy) / 2), size_x - buffer_xy, size_y - buffer_xy, :stroke)	

	#Ticks
	for ab in 0:stepsize:1
		center = Point(ab * length(0:stepsize:1)/(length(0:stepsize:1)+1) * (size_x - buffer_xy) - 0.5 * (size_x - buffer_xy),  0.5 * (size_y - buffer_xy))
		Luxor.line(center, center + Point(0, 10), :stroke)
	end
	fontsize(36)
	for ab in 0:stepsize*2:1
		center = Point(ab * length(0:stepsize:1)/(length(0:stepsize:1)+1) * (size_x - buffer_xy) - 0.5 * (size_x - buffer_xy),  0.5 * (size_y - buffer_xy))
		label(string(ab), :S, center + Point(boxwidth/2, 10))
	end

	for db in 0:stepsize:1
		center = Point(- 0.5 * (size_x - buffer_xy), db * length(0:stepsize:1)/(length(0:stepsize:1)+1) * (size_y - buffer_xy) - 0.5 * (size_y - buffer_xy))
		Luxor.line(center - Point(10, 0), center , :stroke)
	end
	fontsize(36)
	for db in 0:stepsize*2:1
		center = Point(- 0.5 * (size_x - buffer_xy), (1-db) * length(0:stepsize:1)/(length(0:stepsize:1)+1) * (size_y - buffer_xy) - 0.5 * (size_y - buffer_xy))
		label(string(db), :W, center + Point(-10, boxheight/2))
	end

	fontsize(60)
	text("Aggregate balance", Point(0.5 * length(0:stepsize:1)/(length(0:stepsize:1)+1) * (size_x - buffer_xy) - 0.5 * (size_x - buffer_xy),  0.5 * (size_y - buffer_xy)) + buffer_xy/2 - 20, halign=:center, valign = :bottom)
	text("Disaggregate balance", Point(- 0.5 * (size_x - buffer_xy) - buffer_xy/2 + 35, 0.5 * length(0:stepsize:1)/(length(0:stepsize:1)+1) * (size_y - buffer_xy) - 0.5 * (size_y - buffer_xy)), halign=:center,   valign = :middle, angle = -pi/2)

	finish()
	preview()

end

