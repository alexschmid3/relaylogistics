
include("scripts/visualizations/timespacenetwork.jl")

for i in [19]
    #magiparcs = [a for a in magarcs.A[i] if value(x_magip[i,a]) > 1e-4]
    iparcs = [a for a in magarcs.A[i] if value(x_magip[i,a]) > 1e-4]
    basisarcs = lpbasisarcs.A[i]
    maggenarcs = magarcs.A[i]
    removedarcs = setdiff(arcsbeforecolmgmt.A[i], magarcs.A[i])
    basisiparcs = [a for a in lpbasisarcs.A[i] if value(x_bip[i,a]) > 1e-4]

    fractlist = [[],Dict(),[]]

    arclistlist = [maggenarcs, basisarcs, basisiparcs, iparcs]
    colorlist = [(200,200,200), (255,204,129), (255,153,0), (0,0,0)] 
    thicknesslist = [6, 6, 12, 12]
    dashlist = ["solid", "solid", "solid", "solid"]

    timespacenetwork(string("outputs/viz/fig_order", i,".png"), arclistlist, colorlist, thicknesslist, dashlist, fractlist, 2400, 1800)

    arclistlist = [basisarcs, basisiparcs]
    colorlist = [(200,200,200), (0,0,0)] 
    thicknesslist = [6, 12]
    dashlist = ["solid","solid"]

    timespacenetwork(string("outputs/viz/fig_order", i,"_1.png"), arclistlist, colorlist, thicknesslist, dashlist, fractlist, 2400, 1800)

    arclistlist = [maggenarcs,  iparcs]
    colorlist = [(200,200,200), (0,0,0)] 
    thicknesslist = [6, 12]
    dashlist = ["solid", "solid"]

    timespacenetwork(string("outputs/viz/fig_order", i,"_2.png"), arclistlist, colorlist, thicknesslist, dashlist, fractlist, 2400, 1800)

    #Find the drivers
    basisotherorders, basisemptytrips = [], []
    for a in intersect(basisiparcs, primaryarcs.A_space)
        println("----- ARC $a -----")
        arcDesc(a)
        for l in 1:numlocs, s in 1:numshifts, f in fragmentscontaining[l,s,a], d in driversets[l,s]
            if value(z_bip[d,f]) > 1e-4
                println("    Driver $d")
                arclist = []
                for a2 in 1:numarcs
                    if f in fragmentscontaining[l,s,a2]
                        push!(arclist, a2)
                    end
                end
                sort!(arclist, by = x -> nodesLookup[arcLookup[x][1]][2])
                
                for a2 in setdiff(intersect(arclist, primaryarcs.A_space), a)
                    println("    ", a2)
                    ordertrips = 0
                    for j in orders
                        if (a2 in lpbasisarcs.A[j]) && (value(x_bip[j,a2]) > 1e-4)
                            println("     --> Order $j")
                            ordertrips += 1
                        end
                    end
                    if ordertrips == 0
                        push!(basisemptytrips, a2)
                    else    
                        push!(basisotherorders, a2)
                    end
                end
            end
        end
    end

    magotherorders, magemptytrips = [], []
    for a in intersect(iparcs, primaryarcs.A_space)
        println("----- ARC $a -----")
        arcDesc(a)
        for l in 1:numlocs, s in 1:numshifts, f in fragmentscontaining[l,s,a], d in driversets[l,s]
            if value(z_magip[d,f]) > 1e-4
                println("    Driver $d")
                arclist = []
                for a2 in 1:numarcs
                    if f in fragmentscontaining[l,s,a2]
                        push!(arclist, a2)
                    end
                end
                sort!(arclist, by = x -> nodesLookup[arcLookup[x][1]][2])
                
                for a2 in setdiff(intersect(arclist, primaryarcs.A_space), a)
                    println("    ", a2)
                    print("    ")
                    arcDesc(a2)
                    ordertrips = 0
                    for j in orders
                        if (a2 in magarcs.A[j]) && (value(x_magip[j,a2]) > 1e-4)
                            println("     --> Order $j")
                            ordertrips += 1
                        end
                    end
                    if ordertrips == 0
                        push!(magemptytrips, a2)
                    else    
                        push!(magotherorders, a2)
                    end
                end
            end
        end
    end
    
    push!(basisemptytrips, arcs[nodes[19,18], nodes[19,24]])
    push!(basisemptytrips, arcs[nodes[19,24], nodes[19,30]])
    push!(basisemptytrips, arcs[nodes[19,30], nodes[19,36]])
    push!(basisemptytrips, arcs[nodes[10,48], nodes[10,54]])
    push!(basisemptytrips, arcs[nodes[10,54], nodes[10,60]])
    push!(basisemptytrips, arcs[nodes[10,60], nodes[10,66]])
    push!(basisemptytrips, arcs[nodes[10,66], nodes[10,72]])

    arclistlist = [basisiparcs, iparcs, magemptytrips, basisemptytrips, magotherorders, basisotherorders]
    colorlist = [(255,153,0), (0,0,0), (190,190,190), (0,0,200), (190,190,190), (0,0,200)] 
    thicknesslist = [12, 12, 12, 12, 12, 12]
    dashlist = ["solid", "solid", "dashed", "dashed", "solid", "solid"]

    timespacenetwork(string("outputs/viz/fig_order", i,"_drivers.png"), arclistlist, colorlist, thicknesslist, dashlist, fractlist, 2400, 1800)

    arclistlist = [basisiparcs, basisemptytrips, basisotherorders]
    colorlist = [(0,0,0), (254,97,0), (220,38,127), (254,97,0)] 
    thicknesslist = [14, 12, 12]
    dashlist = ["solid", "dashed", "solid"]

    timespacenetwork(string("outputs/viz/fig_order", i,"_drivers1.png"), arclistlist, colorlist, thicknesslist, dashlist, fractlist, 2400, 1800)

    arclistlist = [iparcs, magemptytrips, magotherorders]
    colorlist = [(0,0,0), (254,97,0), (220,38,127)] 
    thicknesslist = [14, 12, 12]
    dashlist = ["solid", "dashed", "solid"]

    timespacenetwork(string("outputs/viz/fig_order", i,"_drivers2.png"), arclistlist, colorlist, thicknesslist, dashlist, fractlist, 2400, 1800)

end

#Blues
#(31,92,229)
#(139,171,241)

#Oranges
#(255,153,0)
#(255,204,129)