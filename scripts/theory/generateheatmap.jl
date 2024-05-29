
using Luxor, Colors, CSV, DataFrames

tempHueList = collect(Colors.color_names)

#Exclude colors that are too light
hueList = []
for item in tempHueList
	if (item[2][1] + item[2][2] + item[2][3] <= 515) & (max(abs(item[2][1]-item[2][2]), abs(item[2][1]-item[2][3]), abs(item[2][2]-item[2][3])) >= 100)
		push!(hueList, item) 
	end
end

#Parameters
stdev_stepsize = 0.1
n_list = [1,2,3,4,5,6,7,8,9,10]
size_x, size_y = 2000,2000

#File names
inputfilename = "outputs/heatmapdata/heatmap1_repos_outputs.csv"
outputfilename = string("figures/heatmap_repos.png")

#---------------------------------------------------------------------------------------#

function generateheatmap(inputfilename, outputfilename, size_x, size_y)

    buffer_xy = 300

	#Find coordinates for each workstation and pod storage location
	boxwidth = (size_x - buffer_xy) / (length(n_list)+1)
	boxheight = (size_y - buffer_xy) / (length(0:stdev_stepsize:1)+1)
	betweensquares_x = (size_x - buffer_xy) / (length(n_list)) - (size_x - buffer_xy) / (length(n_list)+1)
	betweensquares_y = (size_y - buffer_xy) / (length(0:stdev_stepsize:1)) - (size_y - buffer_xy) / (length(0:stdev_stepsize:1)+1)

	boxPoints = Dict()
	for stdv in 0:stdev_stepsize:1, nindex in 1:length(n_list)
		n = n_list[nindex]
		newx = (nindex-1)/(length(n_list)-1) * length(n_list)/(length(n_list)+1) * (size_x - buffer_xy) - 0.5 * (size_x - buffer_xy)
		newy = (1-stdv) * length(0:stdev_stepsize:1)/(length(0:stdev_stepsize:1)+1) * (size_y - buffer_xy) - 0.5 * (size_y - buffer_xy)
		boxPoints[n,stdv] = Point((newx, newy))
	end

	#-------------------------------------------------------------------------#

	#Read the input file
	mapdata = CSV.read(inputfilename, DataFrame)
    filtereddata = filter(:n => n -> n in n_list, mapdata)
	mileslookup, denomlookup, countlookup = Dict(), Dict(), Dict()

	for stdv in 0:stdev_stepsize:1, n in n_list
		mileslookup[n,stdv] = "X"
        denomlookup[n,stdv] = 0
        countlookup[n,stdv] = 0
	end

	allmiles = []
	for row in 1:size(filtereddata)[1]
		n, stdv = abs(convert(Int, round(filtereddata[row,6], digits=0))), abs(round(stdev_stepsize * round(filtereddata[row,7] / stdev_stepsize),digits=2))
		if mileslookup[n,stdv] == "X"
            mileslookup[n,stdv] = 0
        end
        mileslookup[n,stdv] += (filtereddata[row,9] - filtereddata[row,8]) / filtereddata[row,9]
        denomlookup[n,stdv] += 1 #filtereddata[row,9] 
        countlookup[n,stdv] += 1
		#pctlookup[n,stdv] += 100 * (filtereddata[row,7] - filtereddata[row,6]) / filtereddata[row,7] 
		push!(allmiles, (filtereddata[row,9] - filtereddata[row,8])/filtereddata[row,9])
	end	

	maxmiles, minmiles = maximum(allmiles), minimum(allmiles)
    maxmiles = max(maxmiles, -1*minmiles)
    minmiles = min(-1*maxmiles/3, minmiles)

	#-------------------------------------------------------------------------#
	
	#Create a list of boxes and color types
	thickness = 0
	textfontsize = 10

	counter = 0

	xfont = 40
	numberfont = 30

	locationsquares, textsquares = [], []
	for stdv in 0:stdev_stepsize:1, n in n_list
		corner = boxPoints[n,stdv]
		center = boxPoints[n,stdv] + Point(boxwidth/2, boxheight/2) 
		if mileslookup[n,stdv] == "X"
			boxcolor = (200,200,200)
			textcolor = (150,150,150)
			actualtext = "X"
			push!(locationsquares, (corner, boxcolor, boxwidth, boxheight, thickness))
			push!(textsquares, (center, textcolor, actualtext, xfont))
		elseif mileslookup[n,stdv]/denomlookup[n,stdv] >= -0.005
			lamb = mileslookup[n,stdv] / (countlookup[n,stdv] * maxmiles)
			boxcolor = (98 * lamb + 255 * (1-lamb), 151 * lamb + 255 * (1-lamb), 236 * lamb + 255 * (1-lamb))
			textcolor = ((98 * lamb + 255 * (1-lamb)) / 2, (151 * lamb + 255 * (1-lamb)) / 2, (236 * lamb + 255 * (1-lamb)) / 2)
			actualtext = string(convert(Int, -1*round((100 * mileslookup[n,stdv])/denomlookup[n,stdv],digits=0)), "%")
			push!(locationsquares, (corner, boxcolor, boxwidth, boxheight, thickness))
			push!(textsquares, (center, textcolor, actualtext, numberfont))
			counter += 1 
		else
			lamb = mileslookup[n,stdv] / (countlookup[n,stdv] * minmiles)
			boxcolor = (247 * lamb + 255 * (1-lamb), 91 * lamb + 255 * (1-lamb), 95 * lamb + 255 * (1-lamb))
			textcolor = ((247 * lamb + 255 * (1-lamb)) / 2, (91 * lamb + 255 * (1-lamb)) / 2, (95 * lamb + 255 * (1-lamb)) / 2)
			actualtext = string("+",convert(Int,-1*round((100 * mileslookup[n,stdv])/denomlookup[n,stdv],digits=0)), "%")
			push!(locationsquares, (corner, boxcolor, boxwidth, boxheight, thickness))
			push!(textsquares, (center, textcolor, actualtext, numberfont))
			counter += 1
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
	fontsize(30)
	for lbl in textsquares
		fontsize(lbl[4])
		r_val, g_val, b_val = lbl[2]
		setcolor(convert(Colors.HSV, Colors.RGB(r_val/255, g_val/255, b_val/255)))
		Luxor.text(lbl[3], lbl[1], halign=:center,   valign = :middle)
		#label(lbl[3], :0, lbl[1])
	end

	#Chart Labels
	setcolor("black")
	setdash("solid")
	setline(2)

	#Box
	Luxor.rect(Point(-(size_x - buffer_xy) / 2, - (size_y - buffer_xy) / 2), size_x - buffer_xy, size_y - buffer_xy, :stroke)	

	#Ticks and labels
	for nindex in 1:length(n_list)
		n = n_list[nindex]
		xshift = nindex == 1 ? 0 : betweensquares_x / 2
		center = Point((nindex-1)/(length(n_list)-1) * length(n_list)/(length(n_list)+1) * (size_x - buffer_xy) - 0.5 * (size_x - buffer_xy) - xshift, 0.5 * (size_y - buffer_xy))
		Luxor.line(center, center + Point(0, 10), :stroke)
	end
	fontsize(46)
	for nindex in 1:length(n_list)
		n = n_list[nindex]
		center = Point((nindex-1)/(length(n_list)-1) *  length(n_list)/(length(n_list)+1) * (size_x - buffer_xy) - 0.5 * (size_x - buffer_xy),  0.5 * (size_y - buffer_xy))
		label(string(n), :S, center + Point(boxwidth/2, 10))
	end

	for stdev in 0:stdev_stepsize:1
		yshift = stdev == 0 ? 0 : betweensquares_y / 2
		center = Point(- 0.5 * (size_x - buffer_xy), stdev * length(0:stdev_stepsize:1)/(length(0:stdev_stepsize:1)+1) * (size_y - buffer_xy) - 0.5 * (size_y - buffer_xy) - yshift)
		Luxor.line(center - Point(10, 0), center , :stroke)
	end
	fontsize(46)
	for stdev in 0:stdev_stepsize:1
		center = Point(- 0.5 * (size_x - buffer_xy), (1-stdev) * length(0:stdev_stepsize:1)/(length(0:stdev_stepsize:1)+1) * (size_y - buffer_xy) - 0.5 * (size_y - buffer_xy))
		label(string(round(stdev,digits=1)), :W, center + Point(-10, boxheight/2))
	end

	fontsize(60)
	Luxor.text("Number of locations", Point(0.5 * length(0:stdev_stepsize:1)/(length(0:stdev_stepsize:1)+1) * (size_x - buffer_xy) - 0.5 * (size_x - buffer_xy),  0.5 * (size_y - buffer_xy)) + buffer_xy/2 - 20, halign=:center, valign = :bottom)
	Luxor.text("σ / d", Point(- 0.5 * (size_x - buffer_xy) - buffer_xy/2 + 35, 0.5 * length(0:stdev_stepsize:1)/(length(0:stdev_stepsize:1)+1) * (size_y - buffer_xy) - 0.5 * (size_y - buffer_xy)), halign=:center,   valign = :middle, angle = -pi/2)

	finish()
	preview()

end

#---------------------------------------------------------------------------------------#

generateheatmap(inputfilename, outputfilename, size_x, size_y)
