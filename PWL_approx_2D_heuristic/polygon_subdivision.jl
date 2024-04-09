#include("../heuristique triangle+rectangle/polygon_difference.jl")

function emptying_polygon_list(current_region,missing_regions,polygon_list)
    # vide polygon_list dans current_region et missing_regions
    # vérification que polygon_list n'est pas vide
    if length(polygon_list) == 0
        return current_region, missing_regions
    end
    # current_region reçoit le premier polygone
    current_region = pop!(polygon_list)
    # ajout des polygones polygon_list dans missing_regions
    for i in 1:length(polygon_list)
        push!(missing_regions, polygon_list[i])
    end
    return current_region, missing_regions
end

function polygon_subdivision_test(current_region, missing_regions, arguments, threshold_angle = pi/2, last_iter_test = false)
    # subdivision du polygone current_polygon si certaines conditions sont validées et met à jour current_region et missing_regions en conséquence
    # condition : tous les angles sont plus grand que threshold_angle radian (coupe entre les deux plus grand angles)

    # test subdivision
    bonus_plot_infos = arguments["bonus_plot_infos"]
    missing_regions,current_region,last_iteration_candidate = subdivide_polygon(missing_regions, current_region, arguments, threshold_angle, last_iter_test, bonus_plot_infos)
    #println("\naprès subdivide_polygon :\ncurrent_region $current_region\nmissing regions $missing_regions\nlast_iteration_candidate $last_iteration_candidate")

    # vidage de last_iteration_candidate, sans tester ses polygones parce que
    # la subdivision a été faite à la main et pas avec polygon_difference
    current_region,missing_regions = emptying_polygon_list(current_region,missing_regions,last_iteration_candidate)

    return current_region, missing_regions
end

function repare_vertex_on_edge(polygon, eps = 1e-12)
    # repare a kind of numerical error during general_polygon_difference
    # where 3 consecutive vertices (x1,x2,x3) are on the same line at less than 1e-12 say (n.xi-b <= 1e-12)
    # and x3 \in [x1,x2]
    # in this case, cut the middle point x2.
    k = length(polygon)
    for i in 1:k
        i1 = i
        i2 = mod(i,k)+1
        i3 = mod(i+1,k)+1
        x1 = polygon[i1]
        x2 = polygon[i2]
        x3 = polygon[i3]
        n,b = line_characteristics(x1,x2)
        # test if x3 on segment[x1,x2]
        #println(scalar_product(n,x3)-b)
        #println(x1[1], " ", x2[1], " ",x3[1])
        #println(x1[1] <= x2[1] <= x3[1])
        if abs(scalar_product(n,x3)-b) <= eps && min(x1[1],x2[1]) <= x3[1] <= max(x1[1],x2[1])
            #println("repare_vertex_on_edge detected a reparation with\n$polygon\nand vertices\n$x1\n$x2\n$x3")
            # project i3 on [i1,i2]
            proj = polygon[i3] - (scalar_product(n,polygon[i3])-b).*n
            #println("projection: $proj\nof scalar product $(scalar_product(n,proj)-b)")
            polygon[i3] = proj
            # remove i2 from polygon
            P = []
            for i in 1:k
                if i != i2
                    push!(P,polygon[i])
                end
            end
            return P
        end
    end
    return polygon
end

function general_polygon_difference(base_polygon, polygon_to_remove)
    # retourne le polygone intersection de base_polygon et du complémentaire de polygon_to_remove de type liste
    # base_polygon et polygon_to_remove sont de type liste aussi
    magnitude1,precision1 = find_magnitude_and_precision(base_polygon)
    magnitude2,precision2 = find_magnitude_and_precision(polygon_to_remove)
    magnitude = max(magnitude1,magnitude2)
    precision = max(precision1,precision2)
    base_poly = create_polygon(base_polygon, magnitude, precision)
    poly_to_remove = create_polygon(polygon_to_remove, magnitude, precision)
    c = Clip()
    add_path!(c, base_poly, PolyTypeSubject, true)
    add_path!(c, poly_to_remove, PolyTypeClip, true)
    global result, polys = execute(c, ClipTypeDifference, PolyFillTypeEvenOdd, PolyFillTypeEvenOdd)
    #println("Dans general_polygon_difference :\nresult : $result\npolys :\n$polys\nlongueur de polys : $(length(polys))")
    if !result
        #println("result dans polygon_difference = $result\npolys :\n$polys\nbase_polygon :\n$base_polygon\npolygon_to_remove :\n$polygon_to_remove")
        if repare_false_polygon_difference(polys,base_polygon,polygon_to_remove)
            # on retourne res == [] pour dire que pas d'intersection
            return []
        else
            # on retourne une erreur
            error("résultat faux en sortie de polygon_différence")
        end
    end

    corrected = false
    #println("here in the middle of general_polygon_difference")
    #file_diff = open("../heuristique triangle+rectangle/log_diff.txt","a")
    #println(file_diff,"polygon_difference :\nP = $base_polygon\nT = $polygon_to_remove")
    #close(file_diff)

    # met test_fd à true si on est bien dans le problème identifié de Clipper qui est :
    # renvoie 2 polygones pour l'intersection de base_polygon et polygon_to_remove, l'un est base_polygon et l'autre polygon_to_remove à l'envers
    # en plus ça n'arrive que quand polygon_to_remove est un triangle
    global test_fd = false
    if length(polys) == 2
        try
            global test_fd = test_failed_difference(polys,polygon_to_remove)
        catch e
            global polys = [polys[2],polys[1]]
            global test_fd = test_failed_difference(polys,polygon_to_remove)
        end
    end
    #println("test_fd dans general_polygon_difference vaut $test_fd")
    if test_fd
        #println("polygon_difference a fail\nréparation")
        polys = failed_difference_correction(base_polygon,polygon_to_remove)
        corrected = true
    end
    # dernière étape : on reconvertit polys en une liste de polygone (eux-mêmes de type liste) en prenant en compte magnitude et precision
    res = []
    for i in 1:length(polys)
        poly = polys[i]
        pol = []
        for j in 1:length(poly)
            if !corrected
                try
                    push!(pol, [Float64(poly[j].X)/10^(precision),Float64(poly[j].Y)/10^(precision)])
                catch e
                    push!(pol, [Float64(poly[j][1])/10^(precision),Float64(poly[j][2])/10^(precision)])
                end
            else
                try
                    push!(pol, [Float64(poly[j].X),Float64(poly[j].Y)])
                catch e
                    push!(pol, [Float64(poly[j][1]),Float64(poly[j][2])])
                end
            end
        end
        push!(res,pol)
    end

    # try to repare a vertex on an edge
    for i in 1:length(res)
        res[i] = repare_vertex_on_edge(res[i])
    end

    return res
end
