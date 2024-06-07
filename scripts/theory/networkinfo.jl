
function networkinfo()

    W = [i for i in 1:n]
    E = [i for i in n+1:2*n]

    demandlocs = union(W,E)
    allpairs = []
    for i in W, j in E
        push!(allpairs, (i,j))
    end
    for i in E, j in W
        push!(allpairs, (i,j))
    end

    coordinates = zeros(length(E) + length(W) + 2, 2)
    if length(W) == 1
        curr_y = h
    else
        curr_y = 2*h
    end
    for i in W
        coordinates[i,1] = 0
        coordinates[i,2] = curr_y
        curr_y -= 1 / (length(W) - 1)
    end

    if length(W) == 1
        curr_y = h
    else
        curr_y = 2*h
    end
    for i in E
        coordinates[i,1] = 3*w
        coordinates[i,2] = curr_y
        curr_y -= 1 / (length(W) - 1)
    end

    ps1 = length(E) + length(W) + 1
    ps2 = length(E) + length(W) + 2

    coordinates[ps1,1] = w
    coordinates[ps1,2] = h
    coordinates[ps2,1] = 2*w
    coordinates[ps2,2] = h

    pitstops = [ps1,ps2]

    p1,p2 = P[1], P[2]
    alllocs = union(demandlocs, P)

    return coordinates, pitstops, W, E, demandlocs, allpairs, p1, p2, alllocs

end

#---------------------------------------------------------------------------#

function networkinfo_configurable(K, m, q, n)

    coordinates = zeros(length(E) + length(W) + 2, 2)
    if length(W) == 1
        curr_y = h
    else
        curr_y = 2*h
    end
    for i in W
        coordinates[i,1] = 0
        coordinates[i,2] = curr_y
        curr_y -= 1 / (length(W) - 1)
    end

    if length(W) == 1
        curr_y = h
    else
        curr_y = 2*h
    end
    for i in E
        coordinates[i,1] = 3*w
        coordinates[i,2] = curr_y
        curr_y -= 1 / (length(W) - 1)
    end

    ps1 = length(E) + length(W) + 1
    ps2 = length(E) + length(W) + 2

    coordinates[ps1,1] = w
    coordinates[ps1,2] = h
    coordinates[ps2,1] = 2*w
    coordinates[ps2,2] = h

    pitstops = [ps1,ps2]

    return coordinates, pitstops

end