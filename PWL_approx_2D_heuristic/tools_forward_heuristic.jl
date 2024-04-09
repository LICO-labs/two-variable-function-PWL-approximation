include("tools_compute_and_select_corridor_bounding.jl")
include("tools_sort_pieces_by_direction.jl")
include("tools_dichotomy_on_corridor_pieces.jl")

function compute_and_select_corridor_bounding(polygon, first_vertex, forward_direction, arguments)
    # calcul de l'approximation intérieure du corridor sur polygon

    # calcul de l'aire atteignable par une pièce, partant dans la direction forward_direction
    #global to
    #@timeit to "compute_reachable_polygon" reachable_polygon = compute_reachable_polygon_improved(polygon, first_vertex, forward_direction, arguments)
    reachable_polygon = compute_reachable_polygon_improved(polygon, first_vertex, forward_direction, arguments)

    # calcul de l'approximation intérieure
    #global to
    #@timeit to "get_corridor_by_inner_approx" approximate_corridor = get_corridor_by_inner_approximation(reachable_polygon,forward_direction,arguments)
    approximate_corridor = get_corridor_by_inner_approximation(reachable_polygon,forward_direction,arguments)

    # enlevage des morceaux de l'approximation intérieure qui ne sont pas dans polygon
    #global to
    #@timeit to "remove_pieces_outside_polygon" approximate_corridor = remove_pieces_outside_polygon_list_version(approximate_corridor, reachable_polygon)
    approximate_corridor = remove_pieces_outside_polygon_list_version(approximate_corridor, reachable_polygon)

    return reachable_polygon,approximate_corridor
end

function sort_pieces_by_direction(approximate_corridor, forward_direction)
    # retourne les pièces de approximate_corridor triées dans l'ordre croissant de produit scalaire avec le vecteur forward_direction
    # c'est avec le premier sommet atteint avec une avancée dans la direction forward_direction que le produit scalaire est calculé

    # préparation du tri d'approximate_corridor
    n = size(approximate_corridor)[1]
    sort_values = [get_min_scalar_product(approximate_corridor[i,:,:], forward_direction) for i in 1:n]

    # extraction de la permutation pour le tri
    permutation = sortperm(sort_values)

    # application de la permutation du tri
    return approximate_corridor[permutation,:,:]
end

function sort_pieces_by_direction_list_version(approximate_corridor, forward_direction)
    # retourne les pièces de approximate_corridor triées dans l'ordre croissant de produit scalaire avec le vecteur forward_direction
    # c'est avec le premier sommet atteint avec une avancée dans la direction forward_direction que le produit scalaire est calculé

    # préparation du tri d'approximate_corridor
    n = length(approximate_corridor)
    sort_values = [get_min_scalar_product_list_version(approximate_corridor[i], forward_direction) for i in 1:n]

    # extraction de la permutation pour le tri
    permutation = sortperm(sort_values)

    # application de la permutation du tri
    return approximate_corridor[permutation]
end # ok

function dichotomy_on_corridor_pieces(polygon, pieces, forward_direction, arguments)
    # retourne un nombre num voulant dire que les num premières pièces de pieces rentre dans une seule pièce, mais pas num+1
    #println("we are in dichotomy_on_corridor_pieces for cutting_corridor version")

    # récupération des arguments nécessaires à la dichotomie
    ratio = get(arguments, "DSR", 0.125) # ratio of pieces of pieces to start with

    if typeof(pieces) != Array{Float64,3}
        n_piece = length(pieces)
    else
        n_piece = size(pieces)[1]
    end
    # n_test est le nombre de pièces qui va être testé
    n_test = max(1,Int(ceil(ratio*n_piece)))
    n_min = 0
    n_max = n_piece

    # récupère le dernier plan solution de problème de faisabilité
    last_valid_sol = []
    # booléen pour fin de while
    not_terminated = true
    while not_terminated

        #global to
        #@timeit to "2D_LP" sol,b = feasibility_model2(pieces, n_test, get(arguments, "time_limit", 3600), polygon, forward_direction, arguments = arguments)
        sol,b = feasibility_model2(pieces, n_test, get(arguments, "time_limit", 3600), polygon, forward_direction, arguments = arguments)

        if b == 0
            # problème aucun morceau valide
            if n_test == 1
                error("aucun morceau n'est valide dans la dichotomie\npremière pièce : $(pieces[1,:,:])\npour le polygone atteignable $polygon")
            end

            n_max = n_test-1
            n_test = Int(floor((n_min+n_max+1)/2))

        elseif n_test != n_max
            # while pas fini
            n_min = n_test
            n_test = max(Int(floor((n_min+n_max)/2)), n_test+1)
            # cut pieces n_min+1:n_test # it does not work because the min scalar product depends on this cutting
            last_valid_sol = sol
        else
            # n_test == n_max : while fini
            n_min = n_test
            not_terminated = false
            # enregistrement de la solution
            last_valid_sol = sol
        end

        # si dernière itération est un échec, arrêt de l'algo grâce à ce if
        if n_min == n_max
            not_terminated = false
        end
    end
    #println("$n_test valid pieces out of $(length(pieces))")

    return n_test, last_valid_sol
end

function dichotomy_on_corridor_pieces_real_dichotomy(polygon, pieces, forward_direction, arguments, eps = 1e-4) # first expes with eps = 1e-4. It does change the computation time sensibly, but a greater eps would probably increase the number of pieces
    # retourne un float seuil de produit scalaire délimitant polygon dans la direction forward_direction
    #println("we are in dichotomy_on_corridor_pieces_real_dichotomy")

    # récupération des arguments nécessaires à la dichotomie
    ratio = get(arguments, "DSR", 0.125) # ratio of pieces of pieces to start with

    if typeof(pieces) != Array{Float64,3}
        n_piece = length(pieces)
    else
        n_piece = size(pieces)[1]
    end

    # compute max scalar product and min scalar product
    min_sc = minimum([scalar_product(polygon[i], forward_direction) for i in 1:length(polygon)])
    max_sc = maximum([scalar_product(polygon[i], forward_direction) for i in 1:length(polygon)])
    current_min_sc = min_sc
    current_max_sc = max_sc
    current_sc = (min_sc+max_sc)/2

    # compute min scalar product for all pieces
    min_sc_per_piece = [minimum([scalar_product(pieces[i][j][1:2], forward_direction) for j in 1:length(pieces[i])]) for i in 1:length(pieces)]

    #println("démarrage de la dichotomie entre $min_sc et $max_sc de seuil de produit scalaire")

    # récupère le dernier plan solution de problème de faisabilité
    last_valid_sol = []
    # retrieve the index of the last piece valid
    last_valid_last_piece = 0
    # booléen pour fin de while
    not_terminated = true

    while not_terminated
        #global to
        last_piece = find_in_increasing_list_by_dichotomy(min_sc_per_piece, current_sc)
        #@timeit to "2D_LP" sol,b = feasibility_model2_real_dichotomy(pieces[1:last_piece], current_sc, get(arguments, "time_limit", 3600), polygon, forward_direction, arguments = arguments)
        sol,b = feasibility_model2_real_dichotomy(pieces[1:last_piece], current_sc, get(arguments, "time_limit", 3600), polygon, forward_direction, arguments = arguments)

        if b == 0
            # problème aucun morceau valide
            if current_sc-eps < min_sc
                error("the maximal piece seems too small : current value+$eps < minimum value = $min_sc")
                error("aucun morceau n'est valide dans la dichotomie\npremière pièce : $(pieces[1,:,:])\npour le polygone atteignable $polygon")
            end

            current_max_sc = current_sc
            current_sc = (current_min_sc+current_max_sc)/2
            if current_sc == current_min_sc
                # terminate while loop
                current_max_sc = current_sc
                not_terminated = false
            elseif current_sc-eps < current_min_sc
                # terminate while loop next iteration
                current_sc = current_min_sc
                current_max_sc = current_min_sc
            end
        else
            last_valid_sol = sol
            last_valid_last_piece = last_piece
            if current_sc == current_max_sc
                # terminate while loop
                current_min_sc = current_sc
                not_terminated = false
                # enregistrement de la solution
            elseif current_sc+eps < current_max_sc
                # while not finished
                current_min_sc = current_sc
                current_sc = (current_min_sc+current_max_sc)/2
            elseif current_sc+eps >= current_max_sc
                # current_sc really close to current_max_sc, terminate while loop next iteration
                current_min_sc = current_sc
                current_sc = current_max_sc
            else
                # this case should not happen
                error("dichotomy terminated by a non coded case")
            end
        end
    end
    #println("threshold for scalar product $current_sc when it should be between $min_sc and $max_sc")

    # compute the splitting line
    line_equation = [-forward_direction[1], -forward_direction[2], -current_sc]
    return line_equation, last_valid_sol, last_valid_last_piece
end

function splitting_line(approximate_corridor, number_of_valid_pieces, reachable_polygon, forward_direction)
    # return the line equation that limits the valid pieces of approximate_corridor
    # the [a,b,c] returned should represent all [x,y] satisfying ax + by == c

    # calcul du premier point d'une pièce non valide atteint en allant dans la direction forward_direction
    if typeof(approximate_corridor) != Array{Float64,3}
        first_non_valid_piece = approximate_corridor[number_of_valid_pieces+1]
        first_invalid_point = get_argmin_scalar_product_list_version(first_non_valid_piece, forward_direction)
    else
        first_non_valid_piece = approximate_corridor[number_of_valid_pieces+1,:,:]
        first_invalid_point = get_argmin_scalar_product(first_non_valid_piece, forward_direction)
    end

    # calcul de la droite orthogonale à forward_direction passant par first_invalid_point
    # et de vecteur directeur pointant vers le sommet de départ
    # équation de la droite de la forme ax + by = c (ATTENTION : ce n'est pas ax + by + c = 0)
    a,b = -forward_direction
    c = scalar_product(first_invalid_point, [a,b]) # c computed as: all [x,y] satisfying ax + by == c

    return [a,b,c]
end

function evaluate_polygon_in_3D(piece,plane_equation)
    # retourne un polygone en 3D avec les sommets de piece et les hauteurs calculées
    # grâce à l'équation de plan dans plane_equation

    n = length(piece)
    a,b,c = plane_equation
    for i in 1:n
        # a,b,c décrit z=ax+by+c
        vertex = piece[i][1:2]
        push!(vertex,a*piece[i][1]+b*piece[i][2]+c)
        piece[i] = vertex
    end

    return piece
end
