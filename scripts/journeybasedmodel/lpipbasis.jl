
function getnonzeroarcs(x, orderArcSet, A_plus_i, A_minus_i)

    orderArcSet_red, orderArcSet_space_red, A_plus_i_red, A_minus_i_red = Dict(), Dict(), Dict(), Dict()

    for i in orders
        orderArcSet_red[i] = [dummyarc]
        orderArcSet_space_red[i] = []
    end
    for i in orders, n in 1:extendednumnodes
        if n == Origin[i][1]
            A_plus_i_red[i,n] = [dummyarc]
            A_minus_i_red[i,n] = []
        elseif n == last(Destination[i])
            A_plus_i_red[i,n] = []
            A_minus_i_red[i,n] = [dummyarc]
        else 
            A_plus_i_red[i,n] = []
            A_minus_i_red[i,n] = []
        end
    end

    for i in orders
        destinationlocation = nodesLookup[Destination[i][1]][1]
        #push!(Destination[i], extendednodes[destinationlocation, dummyendtime])
        A_minus_i_red[i, extendednodes[destinationlocation, dummyendtime]] = []
        for n2 in N_end
            arc_ext = extendedarcs[n2, extendednodes[destinationlocation, dummyendtime]]
            push!(orderArcSet_red[i], arc_ext)
            push!(A_plus_i_red[i, n2], arc_ext)
            push!(A_minus_i_red[i, extendednodes[destinationlocation, dummyendtime]], arc_ext)
        end
    end

    for i in orders, a in orderArcSet[i]
        if x[i,a] > 1e-4
            orderArcSet_red[i] = union(orderArcSet_red[i], a)
            if a in A_space
                orderArcSet_space_red[i] = union(orderArcSet_space_red[i], a)
            end
        end
    end
    
    #Create A_plus and A_minus lists
    for i in orders, n in 1:numnodes, a in A_plus[n]
        if (a in orderArcSet_red[i]) & !(a in A_plus_i_red[i,n])
            push!(A_plus_i_red[i,n], a)
        end
    end
    for i in orders, n in 1:numnodes, a in A_minus[n]
        if (a in orderArcSet_red[i]) & !(a in A_minus_i_red[i,n])
            push!(A_minus_i_red[i,n], a)
        end
    end

    return orderArcSet_red, orderArcSet_space_red, A_plus_i_red, A_minus_i_red

end

#---------------------------------------------------------------------------------------#

#=
function getBasisHead(model::JuMP.Model)
    grb_model = model.internalModel
    bhead = zeros(Cint, Gurobi.num_constrs(grb_model.inner))
    ret = Gurobi.@grb_ccall(
        getBasisHead,
        Cint,
        (Ptr{Cvoid}, Ptr{Cint}),
        grb_model.inner, bhead
    )
    if ret != 0
        throw(Gurobi.GurobiError(grb_model.env.inner, ret))
    end
    return bhead .+ 1
end

function getrow(j, m)
    ci = m.internalModel.inner
    row = zeros(10)
    ccall((:CPXbinvarow,CPLEX.libcplex),Cint,(Ptr{Void},Ptr{Void},Cint,Ptr{Cdouble}),ci.env.ptr, ci.lp, j,row)
    return row
end

tableau = zeros(3,10)
for j in 0:2
    tableau[j+1,:] = getrow(j)
end

lp = Model(Gurobi.Optimizer)
@variable(lp, x[1:10] >= 0)
@objective(lp, Min, - x[10])
@constraint(lp, x[1] + x[2] + 3*x[5] + x[7] + 3*x[10] <= 20)
@constraint(lp, [1:9], x[1] >= x[2])
optimize!(lp)

ci = unsafe_backend(lp)

MOI.Utilities.attach_optimizer(ci)

bhead = zeros(10,10)
Gurobi.GRBgetBasisHead(ci, bhead)


model = direct_model(Gurobi.Optimizer())
@variable(model, x[1:10] >= 0)
@objective(model, Min, - x[10])
@constraint(model, x[1] + x[2] + 3*x[5] + x[7] + 3*x[10] <= 20)
@constraint(model, [1:9], x[1] >= x[2])
optimize!(model)
grb = unsafe_backend(model)
num_constraints = 10  # Number of rows in constraint matrix


basis = Vector{Cint}(undef, num_constraints)
GRBgetBasisHead(grb, basis)

Gurobi.column(grb, x[1])
=#