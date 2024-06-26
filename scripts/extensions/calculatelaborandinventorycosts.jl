
function calculatelaborandinventorycosts()

    journeylaborcost, journeyinventorycost = Dict(), Dict()

    for l in 1:numlocs, s in 1:numshifts, f in 1:numfragments[l,s]
        labor, inventory = 0, 0
       
        for a in fragmentarcs[l,s,f]

            l1,t1 = nodesLookup[arcLookup[a][1]]
            l2,t2 = nodesLookup[arcLookup[a][2]]
            #Check if labor (during working hours, away from home)
            if (a <= numarcs) & !(t1 in T_off[s]) & ((l != l1) || (l != l2))
                labor += (t2-t1) 
            end

            #Check if inventory (not moving, away from home)
            if (a <= numarcs) & (l1 == l2) & (l1 != l)
                inventory += (t2-t1)
            end 
        end

        journeylaborcost[l,s,f] = labor
        journeyinventorycost[l,s,f] = inventory
    end

    return journeylaborcost, journeyinventorycost

end

