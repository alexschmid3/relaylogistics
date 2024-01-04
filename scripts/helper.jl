
#Removing item from list
function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end

#Print a description of arc a, format: "(startloc, starttime) ==> (endloc, endtime)"
function arcDesc(a)
	println(nodesLookup[arcLookup[a][1]], " ==> ", nodesLookup[arcLookup[a][2]])
end