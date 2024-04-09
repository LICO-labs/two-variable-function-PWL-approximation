using LinA

include("../supporting_julia_files/basic_functions.jl")
include("../supporting_julia_files/pwl_constant_bounding.jl")
include("../supporting_julia_files/estimate_longest_1D_piece.jl")
#include("../heuristique triangle+rectangle/polygon_difference.jl")
include("tools_compute_reachable_polygon.jl")
include("tools_inclined_rectangle_bounding.jl")
include("tools_select_new_piece.jl")
include("debug_corridor.jl")

function compute_reachable_polygon(polygon, first_vertex, forward_direction, arguments)
    # retourne un polygone inclus dans polygon, correspondant à une surestimation de la zone qui peut être dans le corridor
    #  en partant du sommet first_vertex de polygon et en se déplaçant direction forward_direction

    # calcul de l'avancée max dans 3 directions :
    # forward_direction, en suivant l'arète de droite, et en suivant l'arète de gauche
    # puis calcule de la droite orthogonale à forward_direction la plus proche du sommet de départ :
    # c'est cette ligne qui séparera la partie atteignable de polygon du reste (c'est pas optimal mais on s'en contente)

    # récupération des arguments nécessaires
    starting_vertex = polygon[first_vertex]
    err = get(arguments, "err", Absolute(1))
    str_exprf = get(arguments, "str_exprf", "x*x+y*y")
    domain = get(arguments, "domain", get_rectangular_domain_from_polygon(polygon))
    starting_polygon = init_paving_from_domain(domain)

    # calcul des 3 angles
    println("\n----- starting building reachable polygon -----")
    angle_forward_direction = oriented_angle([1.0,0.0],forward_direction)
    println("angle_forward_direction: $angle_forward_direction")
    n = length(polygon)
    s_right = polygon[mod(first_vertex,n)+1]
    s_left = polygon[mod(first_vertex-2,n)+1]
    s_first = starting_vertex
    angle_right_edge = oriented_angle([1.0,0.0], [s_right[1]-s_first[1], s_right[2]-s_first[2]])
    angle_left_edge = oriented_angle([1.0,0.0], [s_left[1]-s_first[1], s_left[2]-s_first[2]])
    alphas = [angle_forward_direction, angle_right_edge, angle_left_edge]

    # calcul des 3 points limites
    limit_points_min = []
    limit_points_max = []
    #global to
    #@timeit to "1D_corridor" p1,pwl1 = end_of_first_piece(starting_vertex,alphas[1],starting_polygon,err,str_exprf,domain,arguments)
    p1,pwl1 = end_of_first_piece(starting_vertex,alphas[1],starting_polygon,err,str_exprf,domain,arguments)
    #println("on avance de $starting_vertex à $p1 avec l'angle $(alphas[1]) et longueur de pwl1 = $(length(pwl1))")
    if length(pwl1) > 1
        push!(limit_points_min,p1)
        println("adding p1 to min")
    else
        push!(limit_points_max,p1)
    end
    #global to
    #@timeit to "1D_corridor" p2,pwl2 = end_of_first_piece(starting_vertex,alphas[2],starting_polygon,err,str_exprf,domain,arguments)
    p2,pwl2 = end_of_first_piece(starting_vertex,alphas[2],starting_polygon,err,str_exprf,domain,arguments)
    if length(pwl2) > 1
        push!(limit_points_min,p2)
        println("adding p2 to min")
    else
        push!(limit_points_max,p2)
    end
    #global to
    #@timeit to "1D_corridor" p3,pwl3 = end_of_first_piece(starting_vertex,alphas[3],starting_polygon,err,str_exprf,domain,arguments)
    p3,pwl3 = end_of_first_piece(starting_vertex,alphas[3],starting_polygon,err,str_exprf,domain,arguments)
    if length(pwl3) > 1
        push!(limit_points_min,p3)
        println("adding p3 to min")
    else
        push!(limit_points_max,p3)
    end

    # calcul de la plus grande avancée possible avec les 3 utilisations d'end_of_first_piece
    if length(limit_points_min) > 0
        # il y a au moins une pwl parmi les 3 qui était en plus qu'un morceau, choix du demi-plan à partir de ceux-là
        # calcul du point donnant le plus petit produit scalaire entre les points de limit_points_min et forward_direction
        max_forward_point = sort(limit_points_min, by = x -> scalar_product(x,forward_direction))[1]
    else
        # les 3 pwl font un seul morceau, choix du demi-plan le plus loin
        # calcul du point donnant le plus grand produit scalaire entre les points de limit_points_max et forward_direction

        # CORRECTION TO THE CODE USED FOR ARTICLE ISCO2023:
        # if no limit_points_min were found, the function should return the whole polygon because this intersection with the biggest half plane can miss some parts of polygon
        ## return polygon

        max_forward_point = sort(limit_points_max, by = x -> scalar_product(x,forward_direction))[end]
    end

    # calcul de la droite orthogonale à forward_direction la plus proche du sommet de départ
    # comme l'équation d'une droite ax + by = c (ATTENTION : ce n'est pas ax + by + c = 0)
    a,b = -forward_direction
    c = scalar_product(max_forward_point, [a,b])
    line_equation = [a,b,c]

    # calcul du polygone atteignable
    reachable_polygon = intersect_polygon_and_half_plane(polygon,line_equation,first_vertex)
    return reachable_polygon
end

function compute_reachable_polygon_improved(polygon, first_vertex, forward_direction, arguments, eps=1e-12)
    # retourne un polygone inclus dans polygon, correspondant à une surestimation de la zone qui peut être dans le corridor
    #  en partant du sommet first_vertex de polygon et en se déplaçant direction forward_direction

    # calcul de l'avancée max dans 3 directions :
    # forward_direction, en suivant l'arète de droite, et en suivant l'arète de gauche
    # puis calcule de la droite orthogonale à forward_direction la plus proche du sommet de départ :
    # c'est cette ligne qui séparera la partie atteignable de polygon du reste (c'est pas optimal mais on s'en contente)

    # récupération des arguments nécessaires
    starting_vertex = polygon[first_vertex]
    err = get(arguments, "err", Absolute(1))
    str_exprf = get(arguments, "str_exprf", "x*x+y*y")
    domain = get(arguments, "domain", get_rectangular_domain_from_polygon(polygon))
    starting_polygon = init_paving_from_domain(domain)

    # calcul des 3 angles
    # find the vertex of polygon that is the furthest away in the direction of forward_direction
    points = copy(polygon)
    furthest_point = sort(points, by = x -> scalar_product(x,forward_direction))[end]
    angle_furthest_vector = oriented_angle([1.0,0.0],furthest_point.-polygon[first_vertex])
    n = length(polygon)
    s_right = polygon[mod(first_vertex,n)+1]
    s_left = polygon[mod(first_vertex-2,n)+1]
    s_first = starting_vertex
    angle_right_edge = oriented_angle([1.0,0.0], [s_right[1]-s_first[1], s_right[2]-s_first[2]])
    angle_left_edge = oriented_angle([1.0,0.0], [s_left[1]-s_first[1], s_left[2]-s_first[2]])
    alphas = [angle_furthest_vector, angle_right_edge, angle_left_edge]

    # calcul des 3 points limites
    limit_points_min = []
    limit_points_max = []
    #global to
    # pwl1 is not a pwl, but a list just made to contain one or two elements depending on if the endpoint is reached in one piece (one element) or not (two elements)
    #@timeit to "1D_corridor" p1,pwl1 = end_of_first_piece(starting_vertex,alphas[1],starting_polygon,err,str_exprf,domain,arguments)
    p1,pwl1 = end_of_first_piece(starting_vertex,alphas[1],starting_polygon,err,str_exprf,domain,arguments)
    # CORRECTION TO THE CODE USED FOR ARTICLE ISCO2023: (changing tests of the 3 if-else statements)
    #if test_point_inside_polygon(polygon, p1, true)
    if scalar_product(p1, forward_direction)+eps < scalar_product(furthest_point, forward_direction)
        push!(limit_points_min,p1)
    else
        push!(limit_points_max,p1)
    end
    #global to
    #@timeit to "1D_corridor" p2,pwl2 = end_of_first_piece(starting_vertex,alphas[2],starting_polygon,err,str_exprf,domain,arguments)
    p2,pwl2 = end_of_first_piece(starting_vertex,alphas[2],starting_polygon,err,str_exprf,domain,arguments)
    if scalar_product(p2, forward_direction)+eps < scalar_product(s_right, forward_direction)
        push!(limit_points_min,p2)
    else
        push!(limit_points_max,p2)
    end
    #global to
    #@timeit to "1D_corridor" p3,pwl3 = end_of_first_piece(starting_vertex,alphas[3],starting_polygon,err,str_exprf,domain,arguments)
    p3,pwl3 = end_of_first_piece(starting_vertex,alphas[3],starting_polygon,err,str_exprf,domain,arguments)
    if scalar_product(p3, forward_direction)+eps < scalar_product(s_left, forward_direction)
        push!(limit_points_min,p3)
    else
        push!(limit_points_max,p3)
    end

    # calcul de la plus grande avancée possible avec les 3 utilisations d'end_of_first_piece
    if length(limit_points_min) > 0
        # il y a au moins une pwl parmi les 3 qui était en plus qu'un morceau, choix du demi-plan à partir de ceux-là
        # calcul du point donnant le plus petit produit scalaire entre les points de limit_points_min et forward_direction
        max_forward_point = sort(limit_points_min, by = x -> scalar_product(x,forward_direction))[1]
    else
        # les 3 pwl font un seul morceau, choix du demi-plan le plus loin
        # calcul du point donnant le plus grand produit scalaire entre les points de limit_points_max et forward_direction

        # CORRECTION TO THE CODE USED FOR ARTICLE ISCO2023:
        # if no limit_points_min were found, the function should return the whole polygon because this intersection with the biggest half plane can miss some parts of polygon
        return polygon

        #max_forward_point = sort(limit_points_max, by = x -> scalar_product(x,forward_direction))[end]
    end

    # calcul de la droite orthogonale à forward_direction la plus proche du sommet de départ
    # comme l'équation d'une droite ax + by = c (ATTENTION : ce n'est pas ax + by + c = 0)
    a,b = -forward_direction
    c = scalar_product(max_forward_point, [a,b])
    line_equation = [a,b,c]

    # calcul du polygone atteignable
    reachable_polygon = intersect_polygon_and_half_plane(polygon,line_equation,first_vertex)

    return reachable_polygon
end

function get_inclined_bounding(polygon, forward_direction, f, nx, ny, corridor_type, bounding_method, arguments)
    # calcul de u et l sur polygon pour que les pièces soient inclinées selon forward_direction

    # changement de f en g et de A,B,C,D en A,B,C,D différent pour que ça aille avec un rectangle
    # contenant polygon le plus petit possible orienté selon la direction d'avancée
    # 1) calcul du rectangle incliné et de son angle theta d'inclinaison
    rectangle,theta = get_inclined_rectangle(polygon,forward_direction)
    # 2) calcul de g et A,B,C,D
    g,domain = get_new_args_for_bounding(f,rectangle,theta)
    global fincl = g
    A,B,C,D = domain
    # 3) calculer approximate_corridor
    approximate_corridor = launch_bounding(g,A,B,C,D,nx,ny,corridor_type,bounding_method,arguments)
    # 4) changer les positions des sommets des pièces en les bonnes (les sous et surestimations sont bonnes déjà)
    approximate_corridor = get_inclined_positions(approximate_corridor,theta)

    return approximate_corridor
end

function get_corridor_by_inner_approximation(polygon,forward_direction,arguments)
    # prepare function launch_bounding with arguments as set containing necessary arguments

    # récupération des arguments nécessaires
    f = get(arguments, "f", x -> x[1]^2+x[2]^2)
    A,B,C,D = get_rectangular_domain_from_polygon(polygon)
    nx,ny = get(arguments, "n_corridor", [1,1])
    corridor_type = get(arguments, "corridor_type", "double")
    bounding_method = get(arguments, "BM", "efficiency_refinement")
    inclined_bounding = get(arguments, "inclined_bounding", false)

    # construction de l'approximation intérieure du corridor

    if inclined_bounding
        #global to
        #@timeit to "pwl_inner_bounding_construction" approximate_corridor = get_inclined_bounding(polygon, forward_direction, f, nx, ny, corridor_type, bounding_method, arguments)
        approximate_corridor = get_inclined_bounding(polygon, forward_direction, f, nx, ny, corridor_type, bounding_method, arguments)
    elseif !inclined_bounding
        # rien à faire en plus
        #global to
        #@timeit to "pwl_inner_bounding_construction" approximate_corridor = launch_bounding(f,A,B,C,D,nx,ny,corridor_type,bounding_method,arguments)
        approximate_corridor = launch_bounding(f,A,B,C,D,nx,ny,corridor_type,bounding_method,arguments)
    else
        error("inclined_bounding ne vaut ni true ni false : $inclined_bounding")
    end
    return approximate_corridor
end

function remove_pieces_outside_polygon(approximate_corridor, polygon)
    # enlève les pièces de approximate_corridor qui ne s'intersectent pas avec polygon

    # attention aux types de polygon et pieces :
    # polygon est une liste de sommets
    # approximate_corridor est une matrice 3D comme la sortie de lipschitz_bounding

    n = size(approximate_corridor)[1]
    pieces_to_keep = []

    # boucle sur toutes les pièces
    for i in 1:n
        # conversion des pièces de approximate_corridor en polygones
        piece = approximate_corridor[i,:,1:2]
        piece = [piece[1,:],piece[2,:],piece[3,:],piece[4,:]]
        if test_polygon_intersection(polygon,piece)
            push!(pieces_to_keep, i)
        end
    end

    # élimination des pièces qui ne s'intersectent pas
    approximate_corridor = approximate_corridor[pieces_to_keep,:,:]

    return approximate_corridor
end

function remove_pieces_outside_polygon_list_version(approximate_corridor, polygon)
    # enlève les pièces de approximate_corridor qui ne s'intersectent pas avec polygon
    # et coupe les parties à l'extérieur de polygon des autres pièces

    # attention aux types de polygon et pieces :
    # polygon est une liste de sommets
    # approximate_corridor est une matrice 3D comme la sortie de lipschitz_bounding

    n = size(approximate_corridor)[1]
    pieces_to_keep = []

    # boucle sur toutes les pièces
    for i in 1:n
        # conversion des pièces de approximate_corridor en polygones
        piece = approximate_corridor[i,:,1:2]
        piece = [piece[1,:],piece[2,:],piece[3,:],piece[4,:]]
        if test_polygon_intersection(polygon,piece)
            push!(pieces_to_keep, i)
        end
    end

    # élimination des pièces qui ne s'intersectent pas
    approximate_corridor = approximate_corridor[pieces_to_keep,:,:]

    # transform to list version
    list_pieces = []
    for i in 1:size(approximate_corridor)[1]
        new_piece = [approximate_corridor[i,j,:] for j in 1:4]
        push!(list_pieces, new_piece)
    end

    cut_wanted = true # change here to deactivate cutting pieces as well as in tools_dichotomy_on_corridor_pieces
    if cut_wanted
        return cut_pieces_outside_polygon(list_pieces, polygon)
    else
        return list_pieces
    end
end

function cut_pieces_outside_polygon(approximate_corridor, polygon)
    # coupage des pièces restantes qui dépassent du polygon

    # intersect polygon in half spaces
    halfplanes = []
    for i in 1:length(polygon)
        x1 = polygon[i]
        x2 = polygon[mod(i,length(polygon))+1]
        norm_vec,val = line_characteristics(x1,x2) # right halfplane is norm_vec.X >= val
        hs_polygon = Polyhedra.HalfSpace(-1 .* norm_vec, -val)
        push!(halfplanes, hs_polygon)
    end

    list_pieces = []
    n = length(approximate_corridor)
    for i in 1:n
        # uses Polyhedra.intersect to intersect pieces[i,:,:] and vpolygon
        new_piece = [approximate_corridor[i][j][1:2] for j in 1:4]

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
            Al,Bl,Cl = get_planar_equation([approximate_corridor[i][j][1:3] for j in 1:3])
            Au,Bu,Cu = get_planar_equation([approximate_corridor[i][j][[1,2,4]] for j in 1:3])
            # intersect
            intersec = polyhedron(vrep(new_piece))
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
            new_piece = [approximate_corridor[i][j] for j in 1:4]
        end

        # add to list_pieces
        if length(new_piece) > 0
            push!(list_pieces, new_piece)
        end
    end

    return list_pieces
end
