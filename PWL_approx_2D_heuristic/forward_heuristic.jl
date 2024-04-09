include("../supporting_julia_files/functions.jl")
include("tools_forward_heuristic.jl")
include("debug_corridor.jl")

function compute_piece_by_forward_heuristic_cutting_corridor(polygon, first_vertex, forward_direction, arguments)
    # calcul d'une pièce dans le corridor partant du sommet first_vertex de polygon et s'étendant dans la direction forward_direction

    # calcul des pièces approximant le corridor sur polygon et enlevage des pièces qui ne s'intersectent pas avec polygon
    #global to
    #@timeit to "corridor_bounding" reachable_polygon,approximate_corridor = compute_and_select_corridor_bounding(polygon, first_vertex, forward_direction, arguments)
    reachable_polygon,approximate_corridor = compute_and_select_corridor_bounding(polygon, first_vertex, forward_direction, arguments)
    # trie les pièces de manière croissante en fonction du produit scalaire entre forward_direction et le sommet le plus proche de cette direction

    #global to
    if typeof(approximate_corridor) != Array{Float64,3}
        n_tot_piece = length(approximate_corridor)
        #@timeit to "sort_pieces" approximate_corridor = sort_pieces_by_direction_list_version(approximate_corridor, forward_direction) # ok
        approximate_corridor = sort_pieces_by_direction_list_version(approximate_corridor, forward_direction)
    else
        n_tot_piece = size(approximate_corridor)[1]
        #@timeit to "sort_pieces" approximate_corridor = sort_pieces_by_direction(approximate_corridor, forward_direction) # changed here
        approximate_corridor = sort_pieces_by_direction(approximate_corridor, forward_direction)
    end
    if n_tot_piece == 0
        println("problème 0 pièces dans approximate_corridor :\nforward_direction $forward_direction\nfirst_vertex $first_vertex")
        println("polygon\n$polygon\nreachable_polygon\n$reachable_polygon")
    end

    # dichotomie sur le nombre de pièces qui peuvent être incluses dans la pièce que l'on construit, et récupération du plan permettant de valider ces pièces
    #global to
    #@timeit to "dichotomy_2D" number_of_valid_pieces, solution_plane = dichotomy_on_corridor_pieces(reachable_polygon, approximate_corridor, forward_direction, arguments)
    number_of_valid_pieces, solution_plane = dichotomy_on_corridor_pieces(reachable_polygon, approximate_corridor, forward_direction, arguments)

    if number_of_valid_pieces != n_tot_piece
        # calcul le demi-plan séparant la partie de polygon qui constitue la nouvelle pièce et le reste de polygon
        line_equation = splitting_line(approximate_corridor, number_of_valid_pieces, reachable_polygon, forward_direction)

        # retourne les deux morceaux
        #global to
        #@timeit to "intersect_polygon_and_half_plane" new_piece = intersect_polygon_and_half_plane(polygon, line_equation, first_vertex)
        new_piece = intersect_polygon_and_half_plane(polygon, line_equation, first_vertex)
    else
        # tout le polygone atteignable est valide
        new_piece = reachable_polygon
    end

    # reconstruction de la troisième coordonnée de new_piece grâce à solution_plane
    new_piece = evaluate_polygon_in_3D(new_piece,solution_plane)

    # check la validité du polygone new_piece
    if number_of_valid_pieces != 1
        #global to
        if typeof(approximate_corridor) == Array{Float64,3}
            #@timeit to "check_piece_match_approx_corridor" count_errors_l,count_errors_u = check_piece_match_approx_corridor(new_piece, approximate_corridor[1:number_of_valid_pieces,:,:])
            count_errors_l,count_errors_u = check_piece_match_approx_corridor(new_piece, approximate_corridor[1:number_of_valid_pieces,:,:])
            if count_errors_l+count_errors_u > 0
                file = open("debug_erreur_contrainte.txt", "w")
                global BAD_CORRIDOR = deepcopy(approximate_corridor)
                println(file,approximate_corridor)
                global BAD_POLYGON = deepcopy(new_piece)
                println(file,new_piece)
                println(new_piece)
                global BAD_NUMBER = number_of_valid_pieces
                println(file,number_of_valid_pieces)
                println(number_of_valid_pieces)
            end
        end
    end

    return new_piece
end

function compute_piece_by_forward_heuristic_real_dichotomy(polygon, first_vertex, forward_direction, arguments)
    # calcul d'une pièce dans le corridor partant du sommet first_vertex de polygon et s'étendant dans la direction forward_direction

    # calcul des pièces approximant le corridor sur polygon et enlevage des pièces qui ne s'intersectent pas avec polygon
    #global to
    #@timeit to "corridor_bounding" reachable_polygon,approximate_corridor = compute_and_select_corridor_bounding(polygon, first_vertex, forward_direction, arguments)
    reachable_polygon,approximate_corridor = compute_and_select_corridor_bounding(polygon, first_vertex, forward_direction, arguments)
    # trie les pièces de manière croissante en fonction du produit scalaire entre forward_direction et le sommet le plus proche de cette direction

    #global to
    if typeof(approximate_corridor) != Array{Float64,3}
        n_tot_piece = length(approximate_corridor)
        #@timeit to "sort_pieces" approximate_corridor = sort_pieces_by_direction_list_version(approximate_corridor, forward_direction) # ok
        approximate_corridor = sort_pieces_by_direction_list_version(approximate_corridor, forward_direction)
    else
        n_tot_piece = size(approximate_corridor)[1]
        #@timeit to "sort_pieces" approximate_corridor = sort_pieces_by_direction(approximate_corridor, forward_direction) # changed here
        approximate_corridor = sort_pieces_by_direction(approximate_corridor, forward_direction)
    end
    if n_tot_piece == 0
        println("problème 0 pièces dans approximate_corridor :\nforward_direction $forward_direction\nfirst_vertex $first_vertex")
        println("polygon\n$polygon\nreachable_polygon\n$reachable_polygon")
    end

    # dichotomie sur le nombre de pièces qui peuvent être incluses dans la pièce que l'on construit, et récupération du plan permettant de valider ces pièces
    #global to
    #@timeit to "dichotomy_2D" line_equation, solution_plane, number_of_valid_pieces = dichotomy_on_corridor_pieces_real_dichotomy(reachable_polygon, approximate_corridor, forward_direction, arguments)
    line_equation, solution_plane, number_of_valid_pieces = dichotomy_on_corridor_pieces_real_dichotomy(reachable_polygon, approximate_corridor, forward_direction, arguments)

    # line_equation is the halfplane separating the part of polygon that is the new piece and the remaining part

    # retourne les deux morceaux
    #global to
    #@timeit to "intersect_polygon_and_half_plane" new_piece = intersect_polygon_and_half_plane(polygon, line_equation, first_vertex)
    new_piece = intersect_polygon_and_half_plane(polygon, line_equation, first_vertex)

    # reconstruction de la troisième coordonnée de new_piece grâce à solution_plane
    new_piece = evaluate_polygon_in_3D(new_piece,solution_plane)

    # check la validité du polygone new_piece
    if number_of_valid_pieces != 1
        #global to
        if typeof(approximate_corridor) == Array{Float64,3}
            #@timeit to "check_piece_match_approx_corridor" count_errors_l,count_errors_u = check_piece_match_approx_corridor(new_piece, approximate_corridor[1:number_of_valid_pieces,:,:])
            count_errors_l,count_errors_u = check_piece_match_approx_corridor(new_piece, approximate_corridor[1:number_of_valid_pieces,:,:])
            if count_errors_l+count_errors_u > 0
                file = open("debug_erreur_contrainte.txt", "w")
                global BAD_CORRIDOR = deepcopy(approximate_corridor)
                println(file,approximate_corridor)
                global BAD_POLYGON = deepcopy(new_piece)
                println(file,new_piece)
                println(new_piece)
                global BAD_NUMBER = number_of_valid_pieces
                println(file,number_of_valid_pieces)
                println(number_of_valid_pieces)
                error("Erreur détectée dans forward_heuristic")
                close(file)
            end
        end
    end

    return new_piece
end
