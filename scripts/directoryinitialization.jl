
#Create output folders
if !(isdir("outputs"))
	mkdir("outputs")
	mkdir("outputs/static")
	#mkdir("outputs/static/figures")
elseif !(isdir("outputs/static/"))
	mkdir("outputs/static/")
	#mkdir("outputs/static/figures")
elseif !(isdir("outputs/static/figures"))
	mkdir("outputs/static/figures")
	#mkdir("outputs/static/figures")
end

#if !(isdir(string(csvfoldername, runid)))
#	mkdir(string(csvfoldername, runid))
#end

#Create visualization folders
if maketimespacevizfiles + makespatialvizfiles + makeadvancedvizfiles + vizflag >= 1
	if !(isdir("visualizations"))
		mkdir("visualizations")
		mkdir("visualizations/static")
		mkdir(vizfoldername)
	elseif !(isdir("visualizations/static/"))
		mkdir("visualizations/static/")
		mkdir(vizfoldername)
	elseif !(isdir(vizfoldername))
		mkdir(vizfoldername)
	end

	#Create subfolders
	subfoldernames = []
	if maketimespacevizfiles == 1
		push!(subfoldernames, "TimeSpaceNetworkMaps")
	end
	if makespatialvizfiles == 1
		push!(subfoldernames, "SpatialNetworkMaps")
	end
	if makeadvancedvizfiles + vizflag >= 1
		push!(subfoldernames, "DriverMaps")
		push!(subfoldernames, "SpatialOrderMaps")
		push!(subfoldernames, "OrderMaps")
	end

	for sf in subfoldernames
		fullfoldername = string(vizfoldername, "/", sf)
		if !(isdir(fullfoldername))
			mkdir(fullfoldername)
		end
	end
end