
function networkinfo()

    coordinates = zeros(length(E) + length(W) + 2, 2)
    curr_y = 1.0 
    for i in W
        coordinates[i,2] = curr_y
        curr_y -= 1 / (length(W) - 1)
    end

    curr_y = 1.0 
    for i in E
        coordinates[i,1] = 3
        coordinates[i,2] = curr_y
        curr_y -= 1 / (length(W) - 1)
    end

    ps1 = length(E) + length(W) + 1
    ps2 = length(E) + length(W) + 2

    coordinates[ps1,1] = 1.0
    coordinates[ps1,2] = 0.5
    coordinates[ps2,1] = 2.0
    coordinates[ps2,2] = 0.5

    pitstops = [ps1,ps2]

    return coordinates, pitstops

end