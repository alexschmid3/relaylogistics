
#Removing item from list
function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end

#Print a description of arc a, format: "(startloc, starttime) ==> (endloc, endtime)"
function arcDesc(a)
	println(nodesLookup[arcLookup[a][1]], " ==> ", nodesLookup[arcLookup[a][2]])
end

#Runs arcDesc for all arcs in a fragment
function fragDesc(l,s,f)
	arclist = []
	for a in 1:numarcs
		if f in fragmentscontaining[l,s,a]
			push!(arclist, a)
		end
	end
	sort!(arclist, by = x -> nodesLookup[arcLookup[x][1]][2])
	for a in arclist
		arcDesc(a)
	end
end