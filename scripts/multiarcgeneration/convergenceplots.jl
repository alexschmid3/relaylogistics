

function write_cg_conv(filename, cg_iter, maximprove, totalorderarcs, totalorderpaths, rmp_obj)

		df = DataFrame(ID = [runid],
			example = [ex],
			week = [weekstart],
			method = [solutionmethod],
			timehorizon = [horizon], 
			iteration = [cg_iter],
			lowerbound = [rmp_obj + maximprove],
			upperbound = [rmp_obj], 
			path_count = [totalorderpaths],
			arc_count = [totalorderarcs]
           )

	if cg_iter == 1
		CSV.write(filename, df)
	else
		CSV.write(filename, df, append=true)
	end

end

#-------------------------------------------------------------------------------#

function findallpaths(A_plus_i, i)

	allpaths = 0
	for n in Origin[i]
		origintotalpaths = countpaths(A_plus_i, i, n)
		allpaths += origintotalpaths
	end

	return allpaths

end

#-------------------------------------------------------------------------------#

function countpaths(A_plus_i, i, n)

	exploredpaths = 0
	for nextarc in setdiff(A_plus_i[i,n], dummyarc)		
		newnode = arcLookup[nextarc][2]
		if newnode in Destination[i]
			exploredpaths += 1
		else
			exploredpaths += countpaths(A_plus_i, i, newnode)
		end
	end

	return exploredpaths

end
