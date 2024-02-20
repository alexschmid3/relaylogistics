
function calcdriveraveragemetrics(z)

    totaldist, totalshifts, totaltime = 0,0,0
    for d in drivers, f in workingfragments[driverHomeLocs[d], drivershift[d]]
        if value(z[d,f]) > 1e-4
            totalshifts += value(z[d,f]) * fragworkinghours[driverHomeLocs[d], drivershift[d],f]

            arclist = []
            for a in primaryarcs.A_space
                if f in fragmentscontaining[driverHomeLocs[d], drivershift[d],a]
                    push!(arclist, a)
                end
            end
            for a in arclist
                l1, l2 = nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1]
                totaldist += value(z[d,f]) * c[a]
                totaltime += value(z[d,f]) * arcLength_raw[l1,l2]
            end
        end
    end

    println("Average distance per shift = ", 12*totaldist / totalshifts)
    println("Average drive time per shift = ", 12*totaltime / totalshifts)

end