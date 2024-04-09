using ForwardDiff, Plots

function check_inclined_bounding(f, approx_corridor, delta)
    # vérifie que u et l définis dans approx_corridor sont bien respectivement des sous-estimations et surestimations de g+-delta
    global ERR = false
    global TAB_ERRORS = []
    for i in 1:size(approx_corridor)[1]
        for j in 1:4
            if f(approx_corridor[i,j,1:2])+delta < approx_corridor[i,j,4]
                println("erreur échantillon ($i,$j) en $(approx_corridor[i,j,1:2]) : f+delta = $(f(approx_corridor[i,j,1:2])+delta) et u = $(approx_corridor[i,j,4])")
                push!(TAB_ERRORS,[i,j,f(approx_corridor[i,j,1:2])+delta,approx_corridor[i,j,4]])
                global ERR = true
            end
            if f(approx_corridor[i,j,1:2])-delta > approx_corridor[i,j,3]
                println("erreur échantillon ($i,$j) en $(approx_corridor[i,j,1:2]) : g-delta = $(f(approx_corridor[i,j,1:2])-delta) et u = $(approx_corridor[i,j,3])")
                push!(TAB_ERRORS,[i,j,f(approx_corridor[i,j,1:2])-delta,approx_corridor[i,j,3]])
                global ERR = true
            end
        end
    end
    if ERR
        global APPROX_CORRIDOR_FALSE = approx_corridor
        error("erreur détectée dans la construction de u et l avec inclined_bouding true")
    end
    return 0
end

function check_per_piece(f,approx_corridor,i,delta)
    for j in 1:4
        println("$(approx_corridor[i,j,1:2])\t$(approx_corridor[i,j,3])\t< $(f(approx_corridor[i,j,1:2]))\t< $(approx_corridor[i,j,4])")
    end
end

function point_in_piece(piece, X, Y, eps = 1e-6)
    # renvoie true si [X,Y] est dans piece et false sinon
    # eps est la marge avec laquelle on accepte des points légèrement à l'extérieur de piece

    # si piece est une double liste, on le cast en array 2D
    #println("type de piece : $(typeof(piece))")
    if typeof(piece) == Vector{Vector{Float64}} || typeof(piece) == Vector{Any}
        #println("ancienne piece : $piece")
        piece = list_to_array(piece)
        #println("nouvelle piece : $piece")
    end

    # n est le nombre de sommet dans piece
    n = size(piece)[1]

    for j in 1:size(piece)[1]
        v1 = piece[j,1:2]
        v2 = piece[mod(j,n)+1,1:2]
        edge = [v2[1]-v1[1], v2[2]-v1[2]]
        normal_vec = [-edge[2], edge[1]]
        normal_vec = normal_vec ./ norm(normal_vec)
        c = scalar_product(normal_vec,v1)
        if !(scalar_product(normal_vec,[X,Y]) >= c-eps)
            return false
        end
    end

    return true
end

function find_piece(approximate_corridor, X, Y)
    # retourne le numéro de pièce de approximate_corridor contenant [X,Y]
    eps = 1e-6
    global LIST_NEAR_OK = []
    global LIST_3 = []
    for i in 1:size(approximate_corridor)[1]
        test = point_in_piece(approximate_corridor[i,:,:],X,Y,1e-6)
        if test
            return i
        end
    end
    return 0
end

function check_piece_match_approx_corridor(T, approx_corridor, eps=1e-5)
    # vérifie que T est bien entre les sommets de u et l définis dans approx_corridor
    #println("start of check_piece_match")
    #println(typeof(T))
    if typeof(T) == Vector{Any}
        #println("taille de T avant changement en liste : $(size(T))")
        #println(T)
        T = array_to_list2(T)
        #println("T changed to list :\n$T")
    end
    a,b,c = get_planar_equation(T)
    #println("end of get_planar_equation")
    count_errors_u = 0
    count_errors_l = 0
    total_count_u = 0
    total_count_l = 0
    list_false_u = []
    list_false_l = []
    n = size(approx_corridor)[1]
    #println("start of the main for loop")
    file = open("debug_wrong_corridor.txt","w")
    for i in 1:n
        # est-ce que la pièce i est partiellement dans T? oui si au moins un de ses sommets est dedans
        test_in_piece = [false,false,false,false]
        for j in 1:4
            test_in_piece[j] = point_in_piece(T, approx_corridor[i,j,1], approx_corridor[i,j,2])
        end
        if sum(test_in_piece[k]*1 for k in 1:4) >= 1
            # vérification de la satisfaction des contraintes discrétisées
            total_count_u += 4
            total_count_l += 4
            for j in 1:4
                l,u = approx_corridor[i,j,3:4]
                valf = a*approx_corridor[i,j,1] + b*approx_corridor[i,j,2] + c
                if valf < l-eps
                    count_errors_l += 1
                    push!(list_false_l, approx_corridor[i,j,1:2])
                    println("en $(approx_corridor[i,j,1:2]) pièce $i, l > f : $l > $valf écart $(l-valf)")
                    println(file,"en $(approx_corridor[i,j,1:2]) pièce $i, l > f : $l > $valf écart $(l-valf)")
                end
                if valf > u+eps
                    count_errors_u += 1
                    push!(list_false_u, approx_corridor[i,j,1:2])
                    println("en $(approx_corridor[i,j,1:2]) pièce $i, u < f : $u < $valf écart $(valf-u)")
                    println(file,"en $(approx_corridor[i,j,1:2]) pièce $i, u < f : $u < $valf écart $(valf-u)")
                end
            end
        end
    end
    close(file)
    if count_errors_u+count_errors_l > 0
        p = plot(title="position des erreurs, rouge u, bleu l",legend=false)
        p = plot_polygon(p,T,:black)
        for i in 1:length(list_false_u)
            errpos = list_false_u[i]
            p = scatter!(p,[errpos[1]],[errpos[2]],color=:red)
        end
        for i in 1:length(list_false_l)
            errpos = list_false_l[i]
            p = scatter!(p,[errpos[1]],[errpos[2]],color=:blue)
        end
        display(p)
        global BAP = deepcopy(approx_corridor)
        error("le corridor n'est pas respecté par u et l, aller voir erreur fichier debug_wrong_corridor.txt\nerreurs en u :$count_errors_u sur $(total_count_u)\nerreurs en l : $count_errors_l sur $(total_count_l)")
    end
    #println("erreurs en l : $count_errors_l sur $(total_count_l)")
    #println("erreurs en u : $count_errors_u sur $(total_count_u)")
    return count_errors_l,count_errors_u
end

# does not work :
function check_piece_match_approx_corridor_inside_corridor(T, approx_corridor, polygon, separating_line, eps=1e-5) # does not work
    # vérifie que T est bien entre les sommets de u et l définis dans approx_corridor

    return "not working"

    #println("start of check_piece_match")
    #println(typeof(T))
    if typeof(T) == Vector{Any}
        #println("taille de T avant changement en liste : $(size(T))")
        #println(T)
        T = array_to_list2(T)
        #println("T changed to list :\n$T")
    end
    a,b,c = get_planar_equation(T)
    #println("end of get_planar_equation")
    count_errors_u = 0
    count_errors_l = 0
    total_count_u = 0
    total_count_l = 0
    list_false_u = []
    list_false_l = []
    n = size(approx_corridor)[1]
    #println("start of the main for loop")
    file = open("debug_wrong_corridor.txt","w")
    for i in 1:n
        # est-ce que la pièce i est partiellement dans T? oui si au moins un de ses sommets est dedans
        test_in_piece = [false,false,false,false]
        for j in 1:4
            test_in_piece[j] = point_in_piece(T, approx_corridor[i,j,1], approx_corridor[i,j,2])
        end
        if sum(test_in_piece[k]*1 for k in 1:4) >= 1
            # vérification de la satisfaction des contraintes discrétisées
            total_count_u += 4
            total_count_l += 4
            for j in 1:4
                point = approx_corridor[i,j,1:2]
                if point in vrep(polyhedron(polygon)) && point[1]*a + point[2]*b >= c
                    l,u = approx_corridor[i,j,3:4]
                    valf = a*approx_corridor[i,j,1] + b*approx_corridor[i,j,2] + c
                    if valf < l-eps
                        count_errors_l += 1
                        push!(list_false_l, approx_corridor[i,j,1:2])
                        println("en $(approx_corridor[i,j,1:2]) pièce $i, l > f : $l > $valf écart $(l-valf)")
                        println(file,"en $(approx_corridor[i,j,1:2]) pièce $i, l > f : $l > $valf écart $(l-valf)")
                    end
                    if valf > u+eps
                        count_errors_u += 1
                        push!(list_false_u, approx_corridor[i,j,1:2])
                        println("en $(approx_corridor[i,j,1:2]) pièce $i, u < f : $u < $valf écart $(valf-u)")
                        println(file,"en $(approx_corridor[i,j,1:2]) pièce $i, u < f : $u < $valf écart $(valf-u)")
                    end
                end
            end
        end
    end
    close(file)
    if count_errors_u+count_errors_l > 0
        p = plot(title="position des erreurs, rouge u, bleu l",legend=false)
        p = plot_polygon(p,T,:black)
        for i in 1:length(list_false_u)
            errpos = list_false_u[i]
            p = scatter!(p,[errpos[1]],[errpos[2]],color=:red)
        end
        for i in 1:length(list_false_l)
            errpos = list_false_l[i]
            p = scatter!(p,[errpos[1]],[errpos[2]],color=:blue)
        end
        display(p)
        global BAP = deepcopy(approx_corridor)
        error("le corridor n'est pas respecté par u et l, aller voir erreur fichier debug_wrong_corridor.txt\nerreurs en u :$count_errors_u sur $(total_count_u)\nerreurs en l : $count_errors_l sur $(total_count_l)")
    end
    #println("erreurs en l : $count_errors_l sur $(total_count_l)")
    #println("erreurs en u : $count_errors_u sur $(total_count_u)")
    return count_errors_l,count_errors_u
end

function get_hessian(f, X, Y)
    # renvoie la fonction hessienne de f, et la hessienne de f en [X,Y]

    H = x -> ForwardDiff.hessian(f,x)
    return H, H([X,Y])
end
