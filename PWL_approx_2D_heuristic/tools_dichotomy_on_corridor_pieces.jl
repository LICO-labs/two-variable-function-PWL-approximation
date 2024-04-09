using Gurobi

function select_useful_constraints(pieces, n_test)
    # renvoie un tableau (n,4) avec seulement les contraintes "utiles" pour le PL que va résoudre feasibility_model
    # c'est-à-dire que les contraintes intervenant aux mêmes positions (x,y) sont détectées et seule la plus forte est gardée

    # L'idée de la grille avec Set va marcher pour ce que je veux sauf que les erreurs d'arrondis dûes au float vont gêner.
    # il faudrait rajouter une ligne dans launch_bounding qui arrondit les valeurs de pieces[i,j,1:2] à suffisamment de
    # décimales près (faire une formule proprement pour être sûr qu'on perd pas d'infos)
    # et après je pourrais finir select_useful_constraints

    # préparation du tri avec des sets en x et sets en y
    X_list = sort!(collect(Set([pieces[i,j,1] for i in 1:n_test for j in 1:4])))
    Y_list = sort!(collect(Set([pieces[i,j,2] for i in 1:n_test for j in 1:4])))

    # déclaration de la matrice M de taille (length(X_list),length(Y_list),2) qui va garder en mémoire seulement les contraintes les plus fortes
    # associées à la position (X_list[i],Y_list[j]) (3e coordonnée : 1 pour l, 2 pour u)
    M = zeros(length(X_list),length(Y_list),2)
    ### ABANDON PARCE QUE AUTRE IDEE (NAIVE VERSION)
end

function select_useful_constraints_naive_version(pieces, eps = 1e-12)
    # version naïve avec construction pas à pas des contraintes fortes
    # réfléchir à un tri de points 2D par une liste croissante en x (permettant recherche dichotomique) et normale en y

    # Base.sort.insorted sera peut-être utile

    cons = []
    n = size(pieces)[1]
    for i in 1:n
        for j in 1:4
            pos = pieces[i,j,1:2]
            end_of_loop = false
            for k in 1:length(cons)
                if sqrt((pos[1]-cons[k][1])^2+(pos[2]-cons[k][2])^2) < eps
                    # test pour l
                    if cons[k][3] < pieces[i,j,3]
                        cons[k][3] = pieces[i,j,3]
                    end
                    if cons[k][4] > pieces[i,j,4]
                        cons[k][4] = pieces[i,j,4]
                    end
                    end_of_loop = true
                    break
                end
            end
            if !end_of_loop
                push!(cons,pieces[i,j,:])
            end
        end
    end
    return cons
end

function dichotomic_search_for_cons(cons, val)
    # recherche dichotomique de l'indice du premier élément de la liste croissante cons à être plus grand que val

    minpos = 1
    maxpos = length(cons)
    pos = Int(floor(minpos/2+maxpos/2))
    while maxpos > minpos
        if cons[pos][1] < val
            # test montée
            minpos = pos+1
            pos = Int(floor(minpos/2+maxpos/2))
        elseif cons[pos][1] >= val
            # test descente
            maxpos = pos
            pos = Int(floor(minpos/2+maxpos/2))
        end
    end
    return pos
end

function select_useful_constraints_dichotomy_version(pieces, eps = 1e-12)
    # même base que naive_version, mais les éléments de cons sont en permanence par ordre croissant de x
    # teste les sommets un par un en cherchant par dichotomie le plus petit indice i plus grand que x-eps
    # cherche si match en position (x,y) jusqu'à x+eps

    n = size(pieces)[1]

    cons = [pieces[1,1,:]]
    n = size(pieces)[1]
    for i in 1:n
        for j in 1:4
            # pour chaque sommet de chaque pièce, on cherche s'il y a une contrainte plus forte en u et en l
            end_of_loop = false
            # recherche de la position en x minimum du sommet (i,j) dans cons par dichotomie
            pos = dichotomic_search_for_cons(cons, pieces[i,j,1]-eps)
            starting_pos = pos
            # tant que la position en x n'a pas dépassé de plus de eps l'abscisse des contraintes de cons
            while cons[pos][1] < pieces[i,j,1]+eps
                # si les positions de pieces[i,j,1:2] et cons[pos][1:2] sont à eps ou moins en norme 2, on teste si une contrainte domine l'autre
                if sqrt((pieces[i,j,1]-cons[pos][1])^2+(pieces[i,j,2]-cons[pos][2])^2) < eps
                    # test pour la dominance en l
                    if cons[pos][3] < pieces[i,j,3]
                        cons[pos][3] = pieces[i,j,3]
                    end
                    # test pour la dominance en u
                    if cons[pos][4] > pieces[i,j,4]
                        cons[pos][4] = pieces[i,j,4]
                    end
                    end_of_loop = true
                    break
                end
                pos += 1
                if length(cons) < pos
                    break
                end
            end
            if !end_of_loop
                newpos = dichotomic_search_for_cons(cons[starting_pos:end], pieces[i,j,1])+starting_pos-1
                if newpos != length(cons)
                    insert!(cons, newpos, pieces[i,j,:])
                elseif pieces[i,j,1] > cons[end][1]
                    push!(cons,pieces[i,j,:])
                else
                    insert!(cons, newpos, pieces[i,j,:])
                end
            end
        end
    end
    return cons
end

function select_useful_constraints_dichotomy_list_version(list_pieces, eps = 1e-12)
    # même base que naive_version, mais les éléments de cons sont en permanence par ordre croissant de x
    # teste les sommets un par un en cherchant par dichotomie le plus petit indice i plus grand que x-eps
    # cherche si match en position (x,y) jusqu'à x+eps

    n = length(list_pieces)

    cons = [list_pieces[1][1]]
    n = length(list_pieces)
    for i in 1:n
        for j in 1:length(list_pieces[i])
            # pour chaque sommet de chaque pièce, on cherche s'il y a une contrainte plus forte en u et en l
            end_of_loop = false
            # recherche de la position en x minimum du sommet (i,j) dans cons par dichotomie
            pos = dichotomic_search_for_cons(cons, list_pieces[i][j][1]-eps)
            starting_pos = pos
            # tant que la position en x n'a pas dépassé de plus de eps l'abscisse des contraintes de cons
            while cons[pos][1] < list_pieces[i][j][1]+eps
                # si les positions de list_pieces[i][j][1:2] et cons[pos][1:2] sont à eps ou moins en norme 2, on teste si une contrainte domine l'autre
                if sqrt((list_pieces[i][j][1]-cons[pos][1])^2+(list_pieces[i][j][2]-cons[pos][2])^2) < eps
                    # test pour la dominance en l
                    if cons[pos][3] < list_pieces[i][j][3]
                        cons[pos][3] = list_pieces[i][j][3]
                    end
                    # test pour la dominance en u
                    if cons[pos][4] > list_pieces[i][j][4]
                        cons[pos][4] = list_pieces[i][j][4]
                    end
                    end_of_loop = true
                    break
                end
                pos += 1
                if length(cons) < pos
                    break
                end
            end
            if !end_of_loop
                newpos = dichotomic_search_for_cons(cons[starting_pos:end], list_pieces[i][j][1])+starting_pos-1
                if newpos != length(cons)
                    insert!(cons, newpos, list_pieces[i][j][:])
                elseif list_pieces[i][j][1] > cons[end][1]
                    push!(cons,list_pieces[i][j][:])
                else
                    insert!(cons, newpos, list_pieces[i][j][:])
                end
            end
        end
    end
    return cons
end

function check_useful_constraint_selection(pieces, cons, eps = 1e-12)
    # test si cons est bien ce que l'on veut
    # stratégie : pour chaque contrainte de pieces, on regarde si une contrainte au moins aussi forte existe
    # il manque le cas où la contrainte entrée dans cons est strictement plus forte que toutes celles de pieces
    # mais ça a pas de raison d'arriver puisque pendant l'écriture de cons on utilise que les valeurs de pieces

    n_pieces = size(pieces)[1]
    n_cons = length(cons)
    for i in 1:n_pieces
        for j in 1:4
            bool_test = false
            for k in 1:n_cons
                if dist2(cons[k][1:2],pieces[i,j,1:2]) < eps
                    if cons[k][3]+eps < pieces[i,j,3] || cons[k][4]-eps > pieces[i,j,4]
                        error("cons pas bon en (i,j,k) = ($i,$j,$k)\ncons : $(cons[k]) et pieces : $(pieces[i,j,:])")
                    end
                    bool_test = true
                    break
                end
            end
            if !bool_test
                error("aucune contrainte trouvée pour (i,j) = ($i,$j) avec pieces = $(pieces[i,j,:])")
            end
        end
    end
    return 0
end

function check_equality_list(l1,l2, eps = 1e-12)
    # retourne true si tous les éléments de la liste l1 ont un élément égal à précision eps près (car vecteur) dans l2
    # sachant que ce sont des vecteurs donc à précision eps près
    # retourne false sinon

    n = length(l1)
    n2 = length(l2)
    if n != n2
        return false
    end
    for i in 1:n
        bool_test = false
        for j in 1:n
            if dist2(l1[i],l2[j]) < eps
                bool_test = true
            end
        end
        if !bool_test
            return false
        end
    end
    return true
end

function cut_pieces_with_separating_line(pieces, separating_line, polygon)
    # return a list of pieces (not always 4 vertices per piece) that are the intersection of pieces[i,:,:],
    # polygon and the halfplane on the lower side of separating_line (with separating_line = [a,b,c] meaning points satisfying ax+by >= c)

    # build halfplane with separating_line
    a,b,c = separating_line # ax+by >= c
    hs = Polyhedra.HalfSpace([-a,-b], -c) # HalfSpace([d,e], f) is the set of (x,y) satisfying dx+ey <= f

    # intersect polygon in halfspaces and halfplane with separating_line
    halfplanes = [hs]
    for i in 1:length(polygon)
        x1 = polygon[i]
        x2 = polygon[mod(i,length(polygon))+1]
        norm_vec,val = line_characteristics(x1,x2) # right halfplane is norm_vec.X >= val
        hs_polygon = Polyhedra.HalfSpace(-1 .* norm_vec, -val)
        push!(halfplanes, hs_polygon)
    end

    # build vertex description of new_polygon

    list_pieces = []
    n = size(pieces)[1]
    for i in 1:n
        # uses Polyhedra.intersect to intersect pieces[i,:,:] and vpolygon
        new_piece = [pieces[i,j,1:2] for j in 1:4]

        # check if at least one vertex of new piece is not inside
        included = true
        for j in 1:4
            if !in(new_piece[j], halfplanes)
                included = false
                break
            end
        end
        # if not included, intersect with Polyhedra.intersect
        if !included
            # save linear function on piece
            Al,Bl,Cl = get_planar_equation([pieces[i,j,1:3] for j in 1:3])
            Au,Bu,Cu = get_planar_equation([pieces[i,j,[1,2,4]] for j in 1:3])
            # intersect
            intersec = copy(polyhedron(vrep(new_piece)))
            for k in 1:length(halfplanes)
                intersec = Polyhedra.intersect(intersec, halfplanes[k])
            end
            # rebuild vrep of new_piece
            new_piece = reorganise_vrep(intersec)
            # rebuild evaluation of l and u on vertices of new_pieces
            for j in 1:length(new_piece)
                push!(new_piece[j], new_piece[j][1]*Al + new_piece[j][2]*Bl + Cl)
                push!(new_piece[j], new_piece[j][1]*Au + new_piece[j][2]*Bu + Cu)
            end
        else # no intersection
            # get evaluation of l and u on vertices of new_pieces
            new_piece = [pieces[i,j,:] for j in 1:4]
        end

        # add to list_pieces
        if length(new_piece) > 0
            push!(list_pieces, new_piece)
        end
    end

    return list_pieces
end

function cut_pieces_with_separating_line_list_version(pieces, separating_line, polygon, eps = 1e-12)
    # return a list of pieces (not always 4 vertices per piece) that are the intersection of pieces[i,:,:],
    # polygon and the halfplane on the lower side of separating_line (with separating_line = [a,b,c] meaning points satisfying ax+by >= c)

    # build halfplane with separating_line
    a,b,c = separating_line # ax+by >= c
    hs = Polyhedra.HalfSpace([-a,-b], -c) # HalfSpace([d,e], f) is the set of (x,y) satisfying dx+ey <= f

    # intersect polygon in halfspaces and halfplane with separating_line
    halfplanes = [hs]

    list_pieces = []
    n = length(pieces)
    for i in 1:n
        # uses Polyhedra.intersect to intersect pieces[i,:,:] and vpolygon
        new_piece = [pieces[i][j][1:2] for j in 1:length(pieces[i])]

        # check if at least one vertex of new piece is not inside
        included = true
        for j in 1:length(pieces[i])
            if new_piece[j][1]*a + new_piece[j][2]*b < c - eps
                included = false
                break
            end
        end
        # if not included, intersect with Polyhedra.intersect
        if !included
            # save linear function on piece
            Al,Bl,Cl = get_planar_equation([pieces[i][j][1:3] for j in 1:3])
            Au,Bu,Cu = get_planar_equation([pieces[i][j][[1,2,4]] for j in 1:3])
            # intersect
            intersec = copy(polyhedron(vrep(new_piece)))
            for k in 1:length(halfplanes)
                intersec = Polyhedra.intersect(intersec, halfplanes[k])
            end
            # rebuild vrep of new_piece
            new_piece = reorganise_vrep(intersec)
            # rebuild evaluation of l and u on vertices of new_pieces
            for j in 1:length(new_piece)
                push!(new_piece[j], new_piece[j][1]*Al + new_piece[j][2]*Bl + Cl)
                push!(new_piece[j], new_piece[j][1]*Au + new_piece[j][2]*Bu + Cu)
            end
        else # no intersection
            # get evaluation of l and u on vertices of new_pieces
            new_piece = [pieces[i][j][:] for j in 1:length(pieces[i])]
        end

        # add to list_pieces
        if length(new_piece) > 0
            push!(list_pieces, new_piece)
        end
    end

    return list_pieces
end

function cut_pieces_and_select_useful_constraints(pieces, n_test, polygon, forward_direction)
    # prepare pieces for feasibility_model2 with cutting parts of pieces unnecessary and selecting useful constraints

    # find separating line for the n_test first pieces of pieces
    if n_test != length(pieces)
        separating_line = splitting_line(pieces, n_test, polygon, forward_direction)
    else
        separating_line = [1,0,-Inf]
    end

    # cut the pieces on the wrong side of the separating line
    list_pieces = cut_pieces_with_separating_line_list_version(pieces[1:n_test,:,:], separating_line, polygon)

    # select useful constraints
    cons = select_useful_constraints_dichotomy_list_version(list_pieces)

    return cons
end

function cut_pieces_and_select_useful_constraints_real_dichotomy(pieces, threshold_c, polygon, forward_direction)
    # prepare pieces for feasibility_model2 with cutting parts of pieces unnecessary and selecting useful constraints

    # find separating line for the n_test first pieces of pieces.
    # [a,b,c] represents: all [x,y] satisfying ax + by <= c
    separating_line = [-forward_direction[1], -forward_direction[2], -threshold_c]

    # cut the pieces on the wrong side of the separating line
    list_pieces = cut_pieces_with_separating_line_list_version(pieces, separating_line, polygon)

    # select useful constraints
    cons = select_useful_constraints_dichotomy_list_version(list_pieces)

    return cons
end

function feasibility_model2(pieces, n_test, time_limit, polygon, forward_direction, additional_constraints = []; arguments = Dict("LP_SOLVER"=>"Gurobi"))
    # construit et résout un modèle LP pour vérifier s'il existe un plan passant à l'intérieur du corridor définit par les n_test premières pièces de pieces
    #global to
    if arguments["num_func"] != "brekelmans"
        cut_wanted = true # change here to deactivate cutting pieces as well as in tools_compute_and_select_corridor_bounding
        if cut_wanted
            #@timeit to "select_useful_constraints" cons = cut_pieces_and_select_useful_constraints(pieces, n_test, polygon, forward_direction)
            cons = cut_pieces_and_select_useful_constraints(pieces, n_test, polygon, forward_direction)
        else
            #@timeit to "select_useful_constraints" cons = select_useful_constraints_dichotomy_list_version(pieces[1:n_test])
            cons = select_useful_constraints_dichotomy_list_version(pieces[1:n_test])
        end
    else
        #@timeit to "no_select_useful_constraints" cons = [pieces[i,j,:] for i in 1:n_test for j in 1:4] # to adapt to list version of pieces
        cons = [pieces[i,j,:] for i in 1:n_test for j in 1:4] # to adapt to list version of pieces
    end

    n = length(cons)

    # définition du modèle
    LP_SOLVER = get(arguments, "LP_SOLVER", "GLPK")
	if LP_SOLVER == "GLPK"
		model = Model(GLPK.Optimizer)
	elseif LP_SOLVER == "Gurobi"
    	model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "LogToConsole", 0)
        set_optimizer_attribute(model, "Threads", 1) # to remove when you don't need to be deterministic
	elseif LP_SOLVER == "Cbc"
		model = Model(Cbc.Optimizer)
	else
		error("a value of $LP_SOLVER for LP_SOLVER is not supported. Only GLPK, Gurobi and Cbc are supported.")
	end
    set_time_limit_sec(model, time_limit)
    @variable(model, a)
    @variable(model, b)
    @variable(model, c)
    @constraint(model, [i=1:n], cons[i][3] <= a*cons[i][1] + b*cons[i][2] + c <= cons[i][4])

    #=# ajout contraintes additionnelles
    for i in 1:length(additional_constraints)
        println("entrée dans l'ajout de contraintes additionnelles avec i = $i")
        l = additional_constraints[i]
        @constraint(model, l[1]*a + l[2]*b + l[3]*c <= l[4])
    end=#

    # résolution
    status = JuMP.optimize!(model)

    # enregistrement de la solution sous la forme [[a,b,c],1] si le modèle est faisable, et [-1,0] sinon
    term_status = JuMP.termination_status(model)
    if term_status == MOI.OPTIMAL
        val_a = JuMP.value.(a)
        val_b = JuMP.value.(b)
        val_c = JuMP.value.(c)
        return [val_a,val_b,val_c], 1
    elseif term_status == MOI.TIME_LIMIT
        error("time_limit atteint dans le modèle LP")
    else
        return -1, 0
    end
end

function feasibility_model2_real_dichotomy(pieces, threshold_c, time_limit, polygon, forward_direction, additional_constraints = []; arguments = Dict("LP_SOLVER"=>"Gurobi"))
    # construit et résout un modèle LP pour vérifier s'il existe un plan passant à l'intérieur du corridor définit par les n_test premières pièces de pieces
    #global to
    if arguments["num_func"] != "brekelmans"
        cut_wanted = true # change here to deactivate cutting pieces as well as in tools_compute_and_select_corridor_bounding
        if cut_wanted
            #@timeit to "select_useful_constraints" cons = cut_pieces_and_select_useful_constraints_real_dichotomy(pieces, threshold_c, polygon, forward_direction)
            cons = cut_pieces_and_select_useful_constraints_real_dichotomy(pieces, threshold_c, polygon, forward_direction)
        else
            #@timeit to "select_useful_constraints" cons = select_useful_constraints_dichotomy_list_version(pieces[1:n_test])
            cons = select_useful_constraints_dichotomy_list_version(pieces[1:n_test])
        end
    else
        #@timeit to "no_select_useful_constraints" cons = [pieces[i,j,:] for i in 1:n_test for j in 1:4] # to adapt to list version of pieces
        cons = [pieces[i,j,:] for i in 1:n_test for j in 1:4]
    end
    n = length(cons)

    # définition du modèle
    LP_SOLVER = get(arguments, "LP_SOLVER", "GLPK")
	if LP_SOLVER == "GLPK"
		model = Model(GLPK.Optimizer)
	elseif LP_SOLVER == "Gurobi"
    	model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "LogToConsole", 0)
        set_optimizer_attribute(model, "Threads", 1) # to remove when you don't need to be deterministic
	elseif LP_SOLVER == "Cbc"
		model = Model(Cbc.Optimizer)
	else
		error("a value of $LP_SOLVER for LP_SOLVER is not supported. Only GLPK, Gurobi and Cbc are supported.")
	end
    set_time_limit_sec(model, time_limit)
    @variable(model, a)
    @variable(model, b)
    @variable(model, c)
    @constraint(model, [i=1:n], cons[i][3] <= a*cons[i][1] + b*cons[i][2] + c <= cons[i][4])

    # résolution
    status = JuMP.optimize!(model)

    # enregistrement de la solution sous la forme [[a,b,c],1] si le modèle est faisable, et [-1,0] sinon
    term_status = JuMP.termination_status(model)
    if term_status == MOI.OPTIMAL
        val_a = JuMP.value.(a)
        val_b = JuMP.value.(b)
        val_c = JuMP.value.(c)
        return [val_a,val_b,val_c], 1
    elseif term_status == MOI.TIME_LIMIT
        error("time_limit atteint dans le modèle LP")
    else
        return -1, 0
    end
end

function find_in_increasing_list_by_dichotomy(l, val)
    # l is a list of floats with increasing value
    # return the max index with l[index] <= val

    # special case: l[1] > val
    if l[1] > val
        error("Error in find_in_increasing_list_by_dichotomy: all values inside list are greater than $val")
    end

    return searchsortedlast(l,val)
end

function feasibility_model(pieces, n_test, additional_constraints = []; arguments = Dict("LP_SOLVER"=>"Gurobi"))
    # construit et résout un modèle LP pour vérifier s'il existe un plan passant à l'intérieur du corridor définit par les n_test premières pièces de pieces

    # définition du modèle
    LP_SOLVER = get(arguments, "LP_SOLVER", "GLPK")
	if LP_SOLVER == "GLPK"
		model = Model(GLPK.Optimizer)
	elseif LP_SOLVER == "Gurobi"
    	model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "LogToConsole", 0)
        set_optimizer_attribute(model, "Threads", 1) # to remove when you don't need to be deterministic
	elseif LP_SOLVER == "Cbc"
		model = Model(Cbc.Optimizer)
	else
		error("a value of $LP_SOLVER for LP_SOLVER is not supported. Only GLPK, Gurobi and Cbc are supported.")
	end
    @variable(model, a)
    @variable(model, b)
    @variable(model, c)
    @constraint(model, [i=1:n_test,j=1:4], pieces[i,j,3] <= a*pieces[i,j,1] + b*pieces[i,j,2] + c <= pieces[i,j,4])

    # ajout contraintes additionnelles
    for i in 1:length(additional_constraints)
        println("entrée dans l'ajout de contraintes additionnelles avec i = $i")
        l = additional_constraints[i]
        @constraint(model, l[1]*a + l[2]*b + l[3]*c <= l[4])
    end

    # résolution
    status = JuMP.optimize!(model)

    # enregistrement de la solution sous la forme [[a,b,c],1] si le modèle est faisable, et [-1,0] sinon
    term_status = JuMP.termination_status(model)
    if term_status == MOI.OPTIMAL
        val_a = JuMP.value.(a)
        val_b = JuMP.value.(b)
        val_c = JuMP.value.(c)
        return [val_a,val_b,val_c], 1
    else
        return -1, 0
    end
end
