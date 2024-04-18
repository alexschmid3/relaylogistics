
function networkinfo()

    coordinates = zeros(length(E) + length(W) + 2, 2)
    curr_y = 2*h
    for i in W
        coordinates[i,1] = 0
        coordinates[i,2] = curr_y
        curr_y -= 1 / (length(W) - 1)
    end

    curr_y = 2*h 
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