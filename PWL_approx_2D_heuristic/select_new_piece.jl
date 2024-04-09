include("forward_heuristic.jl")
include("tools_select_new_piece.jl")

function check_already_computed_polygon(reusable_polygons, first_vertex, current_region, piece_trials, values, arguments)
    # complète piece_trials et values s'il existe un polygone de reusable_polygons dont le premier
    # sommet est first_vertex
    # retourne already_computed à true si c'est le cas, à false sinon

    already_computed = false
    for j in 1:length(reusable_polygons)
        if reusable_polygons[j][1][1:2] == first_vertex[1:2]
            push!(piece_trials, reusable_polygons[j])
            push!(values, evaluate_piece(current_region, reusable_polygons[j], arguments))
            already_computed = true
            break
        end
    end

    return already_computed, piece_trials, values
end

function compute_and_evaluate_piece(current_region, first_vertex, forward_direction, arguments)
    compute_piece = get(arguments, "compute_piece", compute_piece_by_forward_heuristic_cutting_corridor)
    try
        # construit une nouvelle pièce en partant d'un sommet
        #global to
        #@timeit to "compute_piece" candidate_piece = compute_piece(current_region, first_vertex, forward_direction, arguments)
        candidate_piece = compute_piece(current_region, first_vertex, forward_direction, arguments)

        # évalue la nouvelle pièce
        #global to
        #@timeit to "evaluate_piece" value = evaluate_piece(current_region, candidate_piece, arguments)
        value = evaluate_piece(current_region, candidate_piece, arguments)
        return candidate_piece, value
    catch e
        filename = "errors_compute_and_evaluate_piece.txt"
        file = open(filename, "a")
        println(file,"instance : fonction $(get(arguments, "f", "UNK")) err $(get(arguments, "err","UNK")) valeur de inclined_bounding $(get(arguments, "inclined_bounding","UNK"))")
        println(file,"la pièce partant du sommet $first_vertex du polygon $current_region n'a pas pu être construite")
        println(file,"à cause de l'erreur $e")
        println("la pièce partant du sommet $first_vertex du polygon $current_region n'a pas pu être construite")
        println("à cause de l'erreur $e")
        close(file)
        if e == ErrorException("time_limit atteint dans le modèle LP")
            error("time_limit atteint dans le modèle LP")
        end
        return [[]], -1
    end
end

function near60_heuristic(P)
    # retourne l'indice du sommet le plus proche de 60°
    # parcours du premier au dernier donc favorise le premier si égalité
    n = length(P)
    diff_to_60_angles = []
    for i in 1:n
        s = P[i][1:2]
        s_after = P[mod(i,n)+1][1:2]
        s_before = P[mod(i-2,n)+1][1:2]
        e_after = s_after.-s
        e_before = s_before.-s
        push!(diff_to_60_angles, abs(eval_angle(e_after[1],e_after[2],e_before[1],e_before[2])-pi/3))
    end
    return argmin(diff_to_60_angles)
end

function use_inclined_bounding_option(inclined_bounding, current_region, first_vertex, forward_direction, arguments)
    # compute piece according to inclined_bounding value

    if inclined_bounding == "hybrid"
        # computation and evaluation of a piece for inclined_bounding = true and false, then keep the best

        # first evaluation with inclined_bounding true
        arguments["inclined_bounding"] = true
        candidate_piece1,value1 = compute_and_evaluate_piece(current_region, first_vertex, forward_direction, arguments)

        # second evaluation with inclined_bounding false
        arguments["inclined_bounding"] = false
        candidate_piece2,value2 = compute_and_evaluate_piece(current_region, first_vertex, forward_direction, arguments)

        # keep the best of the two
        if value1 >= value2
            candidate_piece,value = candidate_piece1,value1
        else
            candidate_piece,value = candidate_piece2,value2
        end

        # put back "hybrid" in "inclined_bounding" argument of arguments
        arguments["inclined_bounding"] = "hybrid"
    else
        candidate_piece,value = compute_and_evaluate_piece(current_region, first_vertex, forward_direction, arguments)
    end
    return candidate_piece,value
end

function select_new_piece_component_forward_direction(current_region, reusable_polygons, arguments)
    # calcul d'un certain nombre de nouvelles pièces et choix parmi celles-ci
    n = length(current_region)

    # récupération des arguments nécessaires
    inclined_bounding = get(arguments, "inclined_bounding", "false")
    first_vertex_heuristic = get(arguments, "firstv", "all")

    if first_vertex_heuristic == "near60"
        first_vertex = near60_heuristic(current_region)
        #global to
        #@timeit to "get_args_compute_piece" first_vertex,forward_direction = get_args_for_compute_piece(current_region, first_vertex, arguments)
        first_vertex,forward_direction = get_args_for_compute_piece(current_region, first_vertex, arguments)
        candidate_piece,value = use_inclined_bounding_option(inclined_bounding, current_region, first_vertex, forward_direction, arguments)
        return candidate_piece,[]
    elseif first_vertex_heuristic == "all"
        # rangement des tentatives dans piece_trials
        piece_trials = []
        values = []
        for i in 1:n

            # vérification que la pièce pour le sommet i n'a pas déjà été calculé en regardant dans reusable_polygons
            # si le premier sommet d'un polygone de reusable_polygons est le sommet i, alors la pièce est déjà calculée!
            already_computed = false
            #global to
            #@timeit to "check_already_computed_polygon" already_computed, piece_trials, values = check_already_computed_polygon(reusable_polygons, current_region[i], current_region, piece_trials, values, arguments)
            already_computed, piece_trials, values = check_already_computed_polygon(reusable_polygons, current_region[i], current_region, piece_trials, values, arguments)

            if !already_computed
                # calcul des arguments
                #global to
                #@timeit to "get_args_compute_piece" first_vertex,forward_direction = get_args_for_compute_piece(current_region, i, arguments)
                first_vertex,forward_direction = get_args_for_compute_piece(current_region, i, arguments)

                candidate_piece,value = use_inclined_bounding_option(inclined_bounding, current_region, first_vertex, forward_direction, arguments)

                # stocke le nouvel essai dans piece_trials
                push!(piece_trials, candidate_piece)
                push!(values, value)
                # to remove after problem with real_dichotomy solved:
                if candidate_piece == Any[[0.5, 3.5, -13.07225953230288], [0.5, 3.202270507813, -11.45317545554829], [2.062123341603556, 1.851135253906452, 2.415066863151544], [3.968468468468469, 3.5, 1.405849142371888]]
                    println("stopping because end of computation of wanted candidate piece:\n$candidate_piece")
                    2*eee
                end

                # test si la nouvelle pièce couvre tout le polygone
                if length(candidate_piece) >= 3
                    if abs(polygon_surface(current_region)-polygon_surface(candidate_piece)) <= 1e-12
                        return candidate_piece,[]
                    end
                end
            else
                #println("candidate piece $i already computed: $(piece_trials[end])")
            end
        end

        argmaximum = argmax(values)

        # calcul des pièces non choisies
        non_chosen_polygons = copy(piece_trials)
        filter!(e->e≠piece_trials[argmaximum], non_chosen_polygons)
        if values[argmaximum] == -1
            error("aucune pièce retournée à cette itération, l'option 'lipschitz2' peut être la responsable")
        end

        #println("candidate piece $argmaximum added out of $(length(piece_trials)) candidate pieces")
        return piece_trials[argmaximum], non_chosen_polygons
    elseif first_vertex_heuristic == "random"
        first_vertex = rand(1:length(current_region))
        #global to
        #@timeit to "get_args_compute_piece" first_vertex,forward_direction = get_args_for_compute_piece(current_region, first_vertex, arguments)
        first_vertex,forward_direction = get_args_for_compute_piece(current_region, first_vertex, arguments)
        candidate_piece,value = use_inclined_bounding_option(inclined_bounding, current_region, first_vertex, forward_direction, arguments)
        return candidate_piece,[]
    end
end

function get_specific_forward_direction(FDH)
    # extract the specific forward direction inside FDH
    # FDH is in the form "specific_direction_{a}_{b}" where [a,b] is the forward direction
    spl = split(FDH, "_")
    d = [parse(Float64,spl[3]), parse(Float64,spl[4])]
    return d
end

function select_new_piece_specific_forward_direction(current_region, arguments)
    # compute a piece with the specific forward_direction in argument
    FDH = get(arguments, "FDH", "inconnu")
    # extract the forward direction from the string FDH
    forward_direction = get_specific_forward_direction(FDH)

    # find first vertex as vertex of current_region with minimum scalar product with forward_direction
    # specific case for brekelmans: first_vertex = argmin([current_region[i][1]+current_region[i][2] for i in 1:length(current_region)])
    first_vertex = argmin([scalar_product(forward_direction, current_region[i][1:2]) for i in 1:length(current_region)])
    println("first vertex: $(current_region[first_vertex]) and forward_direction $forward_direction")

    # launch compute_piece
    inclined_bounding = get(arguments, "inclined_bounding", "true")
    compute_piece = get(arguments, "compute_piece", compute_piece_by_forward_heuristic_cutting_corridor)
    #@timeit to "compute_piece" candidate_piece = compute_piece(current_region, first_vertex, forward_direction, arguments)
    candidate_piece = compute_piece(current_region, first_vertex, forward_direction, arguments)
    return candidate_piece,[]
end

function select_new_piece(current_region, reusable_polygons, arguments)
    # redirige vers la bonne heuristique grâce à la valeur de l'argument "FDH" de arguments
    FDH = get(arguments, "FDH", "inconnu")
    if FDH == "bisector_direction" || FDH == "max_curvature" || FDH == "mean_along_edge"
        return select_new_piece_component_forward_direction(current_region, reusable_polygons, arguments)
    elseif FDH == "equally_spaced_forward_direction"
        return select_new_piece_random_forward_direction(current_region, arguments)
    elseif occursin("specific_direction", FDH)
        return select_new_piece_specific_forward_direction(current_region, arguments)
    else
        error("argument pour choisir l'heuristique de choix de la nouvelle pièce est inconnu")
    end
end
