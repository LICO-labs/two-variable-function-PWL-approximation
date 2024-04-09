using LinA
using Plots
using AngleBetweenVectors

include("../supporting_julia_files/basic_functions.jl")
include("../supporting_julia_files/plot_polygons.jl")
include("../supporting_julia_files/tools_save_triangulation.jl")
include("../supporting_julia_files/utils_intersection.jl")

function constraint_translator(con)
    # print en clair la contrainte con
    i11, i12, val1, i21, i22, val2, c, str = con
    if val1 < 0
        val1 = -val1
        val2 = -val2
        c = -c
        if str == "<="
            str = ">="
        elseif str == ">="
            str = "<="
        end
    end
    letters = ["x","y"]
    numbers = ["1","2","3"]
    s = ""
    if val1 != 0
        val2 = val2/val1
        c = c/val1
        val1 = 1
        if val1 < 0
            if str == ">="
                str = "<="
            elseif str == "<="
                str = ">="
            end
        end
        s = string(s, "$(val1)$(letters[Int(i11)])$(numbers[Int(i12)]) + ")
    end
    if val2 != 0
        s = string(s, "$(val2)$(letters[Int(i21)])$(numbers[Int(i22)]) ")
    else
        s = s[1:end-2]
    end
    s = string(s, str)
    if c != 0
        s = string(s, " $(-c)")
    else
        s = string(s, " 0")
    end
    return s
    if val1 != 0 && val2 != 0
        return "$(val1)$(letters[int(i12)])$(numbers[int(i11)]) + $(val2)$(letters[int(i22)])$(numbers[int(i21)]) + $c = 0"
    elseif val1 != 0
        return "$(val1)$(letters[int(i12)])$(numbers[int(i11)]) + $c = 0"
    end
end

function print_constraints(cons)
    # affiche les contraintes cons en lisible
    for i in 1:length(cons)
        println(constraint_translator(cons[i]))
    end
end

function expression_on_line2(str_expr2D, v, alpha)
    # retourne une string str_expr1D de str_expr2D pour tout (x,y)
    # f(x,y) est remplacée par f(x(t),y(t)) où x(t) et y(t) sont définies pour t positif
    # avec (x(0),y(0)) = v et (x(t)-x(0),y(t)-y(0)) = (t*cos(alpha),t*(sin(alpha)))
    # ou de manière équivalente x(t) = v[1] + t*cos(alpha); y(t) = v[2] + t*sin(alpha)
    # on remplace les x par des t parce qu'à la fin il ne faut que des x
    # et on modifie les "x" en "X" et les "y" en "Y", parce qu'on ne veut que des majuscules
    if occursin("x",str_expr2D)
        # gestion des 'exp' en les remplaçant temporairement par 'EP'
        str_expr2D = replace(str_expr2D, "exp" => "EP")
        str_expr2D = replace(str_expr2D, "x" => "X")
        str_expr2D = replace(str_expr2D, "EP" => "exp")
    end
    if occursin("y",str_expr2D)
        str_expr2D = replace(str_expr2D, "y" => "Y")
    end
    str_expr1D = replace(str_expr2D, "X" => "t")
    cos_alpha = round(cos(alpha),digits=8)
    if abs(cos_alpha) < 1e-15
        cos_alpha = 0
    end
    str_expr1D = replace(str_expr1D, "t" => "($(v[1])+x*$cos_alpha)")
    sin_alpha = round(sin(alpha),digits=8)
    if abs(sin_alpha) < 1e-15
        sin_alpha = 0
    end
    str_expr1D = replace(str_expr1D, "Y" => "($(v[2])+x*$sin_alpha)")

    # remplacement de 'exp' par '2.718281828' pour test
    #str_expr1D = replace(str_expr1D, "exp" => "2.718281828^")
    return str_expr1D
end

function float_magnitude(x)
    # retourne le nombre de chiffres avant la virgule utilisé dans x
    # 1 est le minimum
    x = abs(x)
    if x < 1
        return 1
    end
    return Int(ceil(log(10,x)))
end

function float_precision(x, seuil = 12)
    # retourne le nombre de chiffres après la virgule utilisé dans x
    val = 2
    x = round(x,digits=seuil)
    while round(x,digits=val) != x
        val += 1
    end
    return val
end

function find_magnitude_and_precision(polygon)
    # retourne la magnitude et la precision à appliquer à tous les sommets pour coder toute la fraction
    magnitude = 0
    precision = 0
    for i in 1:length(polygon)
        #println("dans l'itération $i de find_magnitude_and_precision")
        #println("polygon[i][1] $(polygon[i][1])")
        #println(polygon)
        val1 = Float64(polygon[i][1])
        val2 = Float64(polygon[i][2])
        magnitude = max(float_magnitude(val1),float_magnitude(val2),magnitude)
        precision = max(float_precision(val1),float_precision(val2),precision)
        #println("magnitude : $magnitude")
        #println("précision : $precision")
    end
    return magnitude, precision
end

function create_polygon(polygon, magnitude = -1, precision = -1)
    # retourne un polygone de type IntPoint pour utilisation par la librairie Clipper
    # polygon est une liste de sommets fermée (premier=dernier) ou non, ça n'a pas l'air de changer quelque chose
    poly = Vector{IntPoint}()
    if magnitude == -1 && precision == -1
        magnitude,precision = find_magnitude_and_precision(polygon)
        for i in 1:length(polygon)
            val1 = Float64(polygon[i][1])
            val2 = Float64(polygon[i][2])
            push!(poly, IntPoint(val1,val2,magnitude,precision+magnitude))
        end
    elseif magnitude == -1
        false_magnitude,precision = find_magnitude_and_precision(polygon)
        for i in 1:length(polygon)
            val1 = Float64(polygon[i][1])
            val2 = Float64(polygon[i][2])
            push!(poly, IntPoint(val1,val2,magnitude,precision+magnitude))
        end
    elseif precision == -1
        magnitude,false_precision = find_magnitude_and_precision(polygon)
        for i in 1:length(polygon)
            val1 = Float64(polygon[i][1])
            val2 = Float64(polygon[i][2])
            push!(poly, IntPoint(val1,val2,magnitude,precision+magnitude))
        end
    else
        false_magnitude,false_precision = find_magnitude_and_precision(polygon)
        for i in 1:length(polygon)
            val1 = Float64(polygon[i][1])
            val2 = Float64(polygon[i][2])
            push!(poly, IntPoint(val1,val2,magnitude,precision+magnitude))
        end
    end
    return poly
end

function test_points_in_polygon(polygon, points, prec_points = 12)
    # retourne une liste contenant true aux indices i où points[i] n'est pas à l'extérieur de polygon et false sinon
    mag,prec = find_magnitude_and_precision(polygon)
    # prétraitement de points pour que les points ne soient définis qu'avec prec_points décimales maximum
    list_points = []
    for i in 1:size(points)[1]
        p1 = Int(floor(round(points[i,1],digits=prec_points)*10^(prec_points)))
        p2 = Int(floor(round(points[i,2],digits=prec_points)*10^(prec_points)))
        push!(list_points, [p1,p2])
    end
    poly = create_polygon(polygon,mag,prec_points)
    println(poly)
    println(list_points)
    # remplissage de la liste de test, qui valent 1 si à l'intérieur, -1 si sur les arètes, et 0 si à l'extérieur
    test_inside = []
    for i in 1:length(list_points)
        p = list_points[i]
        val = pointinpolygon(IntPoint(p[1],p[2]), poly)
        if val != 0
            push!(test_inside, true)
        else
            push!(test_inside, false)
        end
    end
    return test_inside
end

function test_failed_difference(diff,T,eps = 1e-2)
    # retourne true si polygon_difference a fail en inversant l'ordre des sommets du polygon_to_remove en [1,3,2]
    #println("longueur de la différence de polygones dans test_failed_difference $(length(diff))")
    if length(diff) != 2
        # si diff n'est pas deux polygones, on est pas dans le problème identifié de Clipper
        return false
    end
    #file_diff = open("log_diff.txt","a")
    #println(file_diff,"test_failed_difference en cours\ndiff $diff\nT $T")
    #println("test_failed_difference en cours\ndiff $diff\nT $T")
    #close(file_diff)
    poly = diff[2]
    n = length(poly)

    # calcul d'un sommet commun à poly (diff[2]) et T, indiqué par cpt, l'indice du sommet dans poly
    cpt = 1
    global power = -1
    global intpower = -1
    for i in 1:n
        # les deux sommets sont les mêmes si la diff en x ET en y en log en base 10 est entière à cause de la sortie en Int64 de poly
        global power = (round((log(poly[i].X)-log(T[1][1]))/log(10),digits = 8) + round((log(poly[i].Y)-log(T[1][2]))/log(10),digits = 8))/2 ## POSSIBLE PROBLEME SI ARETE VERTICALE!!! IL FAUT CALCULER EN Y AUSSI
        #println("power inside try : $power")
        if abs(power-round(power,digits=0)) < 1e-8
            global intpower = Int(power)
            #println("calcul de intpower réussi!\npower $power\nintpower $intpower")
            break
        else
            #println("calcul de intpower raté... power $power")
        end
        cpt += 1
    end

    #println("en sortie de boucle for :\ncpt $cpt\npower $power intpower $intpower")
    #file_diff = open("log_diff.txt","a")
    #println(file_diff,"power $power intpower $intpower poly[cpt] $(poly[cpt]) T11 $(T[1][1]) T11*10^power $(T[1][1]*10^power)")
    #println("poly[cpt] $(poly[cpt]) T11 $(T[1][1]) T11*10^power $(T[1][1]*10^power)")
    #close(file_diff)
    #val1 = norm2([poly[1].X-round(T[1][1],digits=intpower)*10^power, poly[1].Y-round(T[1][2],digits=intpower)*10^power])
    #val2 = norm2([poly[2].X-round(T[3][1],digits=intpower)*10^power, poly[2].Y-round(T[3][2],digits=intpower)*10^power])
    #val3 = norm2([poly[3].X-round(T[2][1],digits=intpower)*10^power, poly[3].Y-round(T[2][2],digits=intpower)*10^power])

    # vérification que poly et T sont bien le même triangle, mais avec poly qui tourne dans le mauvais sens
    i = cpt
    ip1 = mod(i,n)+1
    ip2 = mod(i+1,n)+1
    val1 = norm2([poly[i].X-round(T[1][1],digits=intpower)*10^power, poly[i].Y-round(T[1][2],digits=intpower)*10^power])
    val2 = norm2([poly[ip1].X-round(T[3][1],digits=intpower)*10^power, poly[ip1].Y-round(T[3][2],digits=intpower)*10^power])
    val3 = norm2([poly[ip2].X-round(T[2][1],digits=intpower)*10^power, poly[ip2].Y-round(T[2][2],digits=intpower)*10^power])
    #file_diff = open("log_diff.txt","a")
    #println(file_diff,"val1 $val1 val2 $val2 val3 $val3")
    #println("val1 $val1 val2 $val2 val3 $val3")
    #close(file_diff)
    b1 = (val1 < eps)
    b2 = (val2 < eps)
    b3 = (val3 < eps)
    #file_diff = open("log_diff.txt","a")
    #println(file_diff,"b1 $b1 b2 $b2 b3 $b3")
    #close(file_diff)
    # b1, b2 et b3 sont true si c'est bien le même triangle à l'envers
    return (b1 && b2 && b3)
end

function failed_difference_correction(P, T, eps=1e-10)
    # renvoie la différence corrigée entre P et T parce que Clipper n'y arrive pas
    # étape 1 : trouver le sommet commun s1P / s1T
    nP = length(P)
    nT = length(T)
    #println()
    #println("nP $nP nT $nT")
    s1P = 0
    s1T = 0
    finished = false
    for i in 1:nP
        for j in 1:nT
            if norm2(P[i]-T[j][1:2]) < eps
                finished = true
                s1P = i
                s1T = j
                break
            end
        end
        if finished
            break
        end
    end

    # étape 2 : trouver sur quel arète e1 de P est le sommet s1T+1 de T
    s1Tp1 = mod(s1T,nT)+1
    s = T[s1Tp1][1:2]
    e1 = 0
    #println("s1P $s1P s1T $s1T")
    #println("P[s1P+1] = $(P[mod(s1P,nP)+1])\tT[s1T+1] = $(T[mod(s1T,nT)+1])")
    #println("calcul de e1 avec s $s et P :\n$P")
    for i in 1:nP
        edge = [P[i],P[mod(i,nP)+1]]
        #println("itération $i : edge $edge")
        #println("vertex_on_edge(s,edge) = $(vertex_on_edge(s,edge))")
        if vertex_on_edge(s,edge)[1]
            e1 = i
            break
        end
    end

    # étape 3 : trouver sur quel arète e2 de P est le sommet s1T+2 de T
    s1Tp2 = mod(s1T+1,nT)+1
    s = T[s1Tp2][1:2]
    e2 = 0
    e2_found = false
    for i in 1:nP
        edge = [P[i],P[mod(i,nP)+1]]
        if vertex_on_edge(s,edge)[1]
            e2 = i
            e2_found = true
            break
        end
    end
    # cas particulier : si e2_found == false, alors T[s1Tp2] n'est pas sur une arète de P,
    # et donc il suffit d'ajouter T[s1Tp2] dans P juste après P[s1P]!
    if !e2_found
        new_P = P[1:s1P]
        push!(new_P, s)
        append!(new_P, P[(s1P+1):end])
        return [new_P]
    end

    # étape 4 : construire la différence de polygone à partir de s1P,s1T,e1,e2 et P et T
    # on vérifie que e1 et e2 soit bien adjacent dans P, sinon argh
    if abs(e1-e2) == 1 || (min(e1,e2) == 1 && max(e1,e2) == nP)
        #println("j'ai codé ce cas")
    else
        error("j'ai pas codé ce cas :\nP $P\nT $T\ns1P $s1P s1T $s1T e1 $e1 e2 $e2 s1Tp1 $s1Tp1 s1Tp2 $s1Tp2")
    end
    # sommets de T qui seront dans la différence de polygone
    Tp1 = T[s1Tp1][1:2]
    Tp1 = vertex_on_edge(Tp1,[P[e1],P[mod(e1,nP)+1]])[2]
    Tp2 = T[s1Tp2][1:2]
    Tp2 = vertex_on_edge(Tp2,[P[e2],P[mod(e2,nP)+1]])[2]
    #file_diff = open("log_diff.txt","a")
    #println(file_diff,"Tp1 $Tp1 Tp2 $Tp2")
    #close(file_diff)
    # on réorganise P en PP pour que P[s1P] corresponde à PP[1] en faisant une rotation sur les sommets
    PP = deepcopy(P[s1P:nP])
    for i in 1:(s1P-1)
        push!(PP,P[i])
    end
    #file_diff = open("log_diff.txt","a")
    #println(file_diff,"P $P")
    #println(file_diff,"PP $PP")
    #close(file_diff)
    #rot = s1P-1
    #e1 = mod(e1-rot-1,nP)+1
    #e2 = mod(e2-rot-1,nP)+1
    # et on calcule la différence de polygone :
    # 1) on ajoute Tp2 et Tp1 à poly
    poly = [Tp2,Tp1]
    # 2) on ajoute le sommet en fin de e1 s'il n'est pas confondu avec Tp1
    #println(PP[2][1])
    #println([Tp1[1]-PP[2][1], Tp1[2]-PP[2][2]])
    if norm2([Tp1[1]-PP[2][1], Tp1[2]-PP[2][2]]) > eps
        # Tp1 est différent du sommet suivant s1P alors on ajoute le sommet suivant s1P
        push!(poly,PP[2])
    end
    # 3) on ajoute le gros des sommets (entre 3 et nP-1)
    for i in 3:(nP-1)
        push!(poly,PP[i])
    end
    # 4) on ajoute le sommet en début de e2 s'il n'est pas confondu avec Tp2
    if norm2([Tp2[1]-P[end][1], Tp2[2]-P[end][2]]) > eps
        push!(poly,PP[end])
    end
    #file_diff = open("log_diff.txt","a")
    #println(file_diff,"poly :\n$poly")
    #close(file_diff)
    return [poly]
end

function polygon_rotation(P,rot = 1)
    # fait une rotation de rot éléments sur le polygone P
    for i in 1:rot
        temp = P[2:end]
        push!(temp,P[1])
        P = temp
    end
    return P
end

function repare_false_polygon_difference(polys,P,T,eps = 1e-8)
    # regarde si P et T coïncident et retourne true si oui, false si non
    if length(P) == length(T)
        n = length(P)
        test = zeros(n)
        for i in 1:n
            # pour chaque sommet de P, on rotationne P d'un cran, et on teste
            # si c'est égal à T à eps près par coordonnée
            bool = test_equal_polygon(P,T)
            if bool
                return true
            end
            P = polygon_rotation(P)
        end
    end
    return false
end

function edge_intersection(edge1, edge2, eps = 1e-8, eps2 = 1e-12, verbose = false)
    # retourne un booléen valant true si intersection et false sinon, plus un point valant l'intersection si intersection
    # retourne false si les deux arètes sont parallèles parce que ça fait plus d'une intersection
    s11 = edge1[1]
    s12 = edge1[2]
    s21 = edge2[1]
    s22 = edge2[2]
    n1, b1 = line_characteristics(s11,s12)
    n2, b2 = line_characteristics(s21,s22)
    if verbose
        println()
        println("charac line 1 : $s11\t$s12\t$n1\t$b1")
        println("charac line 1 : $s21\t$s22\t$n2\t$b2")
    end
    # est-ce que les deux sommets de edge2 sont de part et d'autres de edge1?
    scalprod21 = n1[1]*s21[1] + n1[2]*s21[2]
    scalprod22 = n1[1]*s22[1] + n1[2]*s22[2]
    if (scalprod21-b1)*(scalprod22-b1) <= eps2
        bool2 = true
    else
        bool2 = false
    end
    # est-ce que les deux sommets de edge1 sont de part et d'autres de edge2?
    scalprod11 = n2[1]*s11[1] + n2[2]*s11[2]
    scalprod12 = n2[1]*s12[1] + n2[2]*s12[2]
    if (scalprod11-b2)*(scalprod12-b2) <= eps2
        bool1 = true
    else
        bool1 = false
    end
    #println()
    if verbose
        println("($scalprod11,$scalprod12,$scalprod21,$scalprod22)\t($(scalprod11-b2),$(scalprod12-b2),$(scalprod21-b1),$(scalprod22-b1))\t($(((scalprod11-b2)*(scalprod12-b2))),$((scalprod21-b1)*(scalprod22-b1)))")
        println("($bool1,$bool2)")
    end
    # intersection si bool1 et bool2 valent true
    if bool1 && bool2
        # détection du cas particulier arètes parallèles et retourne false si c'est le cas
        if abs(n1[1]*(-n2[2]) + n1[2]*n2[1]) < eps
            return false, [Inf,Inf]
        end
        # calcul de l'intersection
        A = [n1[1]-n2[1] n1[2]-n2[2];n1[1] n1[2]]
        # formule pour B de signe opposé à ce qu'il y a sur la feuille parce que la formule devrait être ax+by-c=0 avec notre calcul de c (b dans le code)
        B = [b1-b2, b1]
        if verbose
            println("taille de A : $(size(A))\tA[1,2] = $(A[1,2])\tn1[2]-n2[2] = $(n1[2]-n2[2])\tn1[2] = $(n1[2])")
            println("n1 $n1\tn2 $n2")
            println("A et B :\n$A\n$B")
        end
        sol = A\B
        return (true, sol)
    else
        # fin de la fonction
        return false, zeros(2)
    end
end

function end_of_first_piece(v1::Vector{Float64},alpha::Float64,starting_polygon::Vector{Vector{Float64}},err::LinA.ErrorType,str_expr_f::String,limits, args, eps::Float64 = 1e-8)
    # DO NOT CHANGE IT for computing reachable polygon because other functions depend on it (get_mean_along_edge for example)

    #println("\n----- entrée end_of_first_piece -----")
    # tire un rayon et renvoie la fin du premier segment point
    A,B,C,D = limits
    # vecteur normal vec à la droite
    vec = [-sin(alpha),cos(alpha)]
    # coordonnées de la droite passant par v1 et d'angle avec l'axe des abscisses de alpha = [vec[1],vec[2],v1.vec]
    line = [vec[1], vec[2], -vec[1]*v1[1]-vec[2]*v1[2]]
    # calcul de l'expression expr_f1D de la fonction f sur la droite line
    str_expr_f1D = expression_on_line2(str_expr_f, v1, alpha)
    #println()
    #println("expression à la sortie de expression_on_line2 dans end_of_first_piece : $str_expr_f1D")
    # t1 et t2 donnent les bornes sur lesquelles calculer la première pièce
    t1 = 0.0
    # t2 demande de calculer l'intersection entre la demi-droite partant de v1 et d'angle alpha
    # avec les limites du domaine qui sont dans starting_polygon
    bigM = sqrt((B-A)^2+(D-C)^2)
    edge1 = [v1,[v1[1]+cos(alpha)*bigM, v1[2]+sin(alpha)*bigM]]
    #println("alpha : $alpha")
    #println("edge1 : $edge1")
    for i in 1:4
        edge2 = [starting_polygon[i], starting_polygon[mod(i,4)+1]]
        bool, intersec = edge_intersection(edge1,edge2)
        #println("itération $i : edge2 $edge2\t(bool,intersec) = ($bool,$intersec)")
        # test si edge1 et edge2 sont parallèles parce que ça compte pas comme une intersection
        vec1 = [edge1[2][1]-edge1[1][1], edge1[2][2]-edge1[1][2]]
        vec2 = [edge2[2][1]-edge2[1][1], edge2[2][2]-edge2[1][2]]
        #println("vec1 $vec1\nvec2 $vec2\nedge1 $edge1\nedge2 $edge2\nbool $bool\nintersec $intersec")
        if abs(vec1[1]*vec2[2]-vec1[2]*vec2[1])/norm2(vec1)/norm2(vec2) < eps
            bool = false
        end
        if bool
            #println("entrée dans le if parce que bool = $bool")
            # il faut vérifier si l'intersection est bien avec la demi-droite partant de v1 et d'angle alpha
            if intersec[1]*cos(alpha)+intersec[2]*sin(alpha) > v1[1]*cos(alpha)+v1[2]*sin(alpha)
                #println("entrée dans le deuxième if car intersec[1]*cos(alpha)+intersec[2]*sin(alpha) > v1[1]*cos(alpha)+v1[2]*sin(alpha) ($(intersec[1]*cos(alpha)+intersec[2]*sin(alpha)) > $(v1[1]*cos(alpha)+v1[2]*sin(alpha)))")
                global t2 = norm2([v1[1]-intersec[1], v1[2]-intersec[2]])
                #println("t2 = $t2 avec v1 = $v1 et intersec = $intersec")
                #println("définition de t2!")
                if abs(t2) >= eps
                    #println("break")
                    break
                end
            end
        end
    end
    # récupération de l'expression 1D de f
    expr_f1D::Expr = Meta.parse(str_expr_f1D)
    # calcul de la pwl grâce à LinA
    # On peut gagner du temps en ne calculant que la première pièce de la pwl,
    # mais il faudrait regarder un peu en détail l'implémentation
    #println("avant le lancement de LinA :\nexpr_f1D = $expr_f1D\nt1 = $t1\nt2 = $t2\nerr = $err\nalpha = $alpha")
    global time_exactLin
    #println("expr_f1D $expr_f1D t1 $t1 t2 $t2 err $err")
    #println("là on va rentrer dans le try de end_of_first_piece")

    # surestimation du point le plus loin atteignable dans le corridor
    X1 = v1[1:2]
    X2 = X1 .+ t2.*[cos(alpha),sin(alpha)]
    X = estimate_longest_1D_corridor(X1,X2,args)

    # calcul du retour en fonction de si X == X2 ou non
    if dist2(X,X2) < eps
        return X,[[]]
    else
        return X,[[],[]]
    end

    #=pwl = launch_exactLin(expr_f1D,t1,t2,err,alpha,starting_polygon,limits)

    println("longueur de pwl : $(length(pwl))")
    #pwl = LinA.ExactLin(expr_f1D,t1,t2,err,HeuristicLin())
    # extraction de la fin de la première pièce
    tMax = pwl[1].xMax
    point = [v1[1]+tMax*cos(alpha), v1[2]+tMax*sin(alpha)]

    # enregistre les deux avancées en norme dans un fichier pour comparaison
    val_me = dist2(X1,X)
    val_lina = dist2(X1,point)
    s = "moi $val_me LinA $val_lina ratio $(val_me/val_lina)"
    println(s)
    file_spec = open("compare_1D_piece.txt", "a")
    println(file_spec, s)
    close(file_spec)

    return point, pwl=#
end

function check_polygon(P, eps = 1e-6)
    # vérifie si le polygone P ne possède pas deux segments qui se touchent
    # sur plus qu'un point et si ce n'est pas le cas, supprime cette portion
    # retourne le polygone corrigé
    n = length(P)
    i = 1
    while i <= length(P)
        e1 = [P[i],P[mod(i,n)+1]]
        e2 = [P[mod(i,n)+1],P[mod(i+1,n)+1]]
        n1,b1 = line_characteristics(e1[1],e1[2])
        n2,b2 = line_characteristics(e2[1],e2[2])
        scal1 = scalar_product(e1[1],n2)
        scal2 = scalar_product(e1[2],n2)
        #println("itération $i : $(abs(scal1-b2)), $(abs(scal2-b2))")
        if abs(scal1-b2) < eps && abs(scal2-b2) < eps
            #=println("détection de deux segments partiellement confondus entre les segments des sommets $i $(i+1) et $(i+1) $(i+2)")
            println("...réparation...")
            file = open("debug_suppression_ratee_sommet_double.txt", "w")
            println(file,"détection de deux segments partiellement confondus entre les segments des sommets $i $(i+1) et $(i+1) $(i+2)")
            println(file,"...réparation...")
            close(file)=#
            # étape 1 : trouver les deux sommets qui forme la portion commune aux deux segments i et i+1
            # premier sommet facile : c'est une extrémité des deux segments, donc c'est e1[2]=e2[1]
            s1 = e1[2]
            # deuxième sommet: c'est le plus proche de s1 parmi e1[1] et e2[2] (et troisième sommet l'autre)
            v1 = [s1[1]-e1[1][1], s1[2]-e1[1][2]]
            v2 = [s1[1]-e2[2][1], s1[2]-e2[2][2]]
            if norm2(v1) <= norm2(v2)
                s2 = e1[1]
                s3 = e2[2]
                del_e1 = true
            else
                s2 = e2[2]
                s3 = e1[1]
                del_e1 = false
            end
            # étape 2 : éliminer le segment qui est la partie commune et le remplacement
            # de l'autre segment par [s2,s3] ou [s3,s2] se fait tout seul
            deleteat!(P, mod(i,n)+1)
            # étape 3 : changer le compteur i pour recommencer la vérification
            i = 0
            n = length(P)
        end
        i += 1
    end
    # cas particulier où P est maintenant un seul point
    # (parce que c'était un quadrilatère avec deux côtés de longueur ridicule)
    if length(P) == 1
        P = []
        file_diff = open("log_diff.txt", "a")
        println(file_diff,"remplacement de P par [] car c'est un point seul")
        close(file_diff)
    end
    return P
end

function subdivide_polygon(P_remainder, P, args, angle_seuil = pi/2, last_iter_test = true, bonus_plot_infos = [])
    # subdivise P en deux si :
    # 1) soit tous les angles sont plus grand que 90° (coupe entre les deux plus grand angles)
    # (not implemented; for now it errors out) 2) soit P n'est pas convexe (coupe entre le sommet qui n'est pas un point extrême et le plus grand angle non voisin)
    n = length(P)
    angles = get_angles(P)
    # liste qui reste vide si pas de subdivision, et qui prend les deux polygones subidivisé pour test last_iteration
    last_iteration_candidate = []
    # si P == [], on quitte maintenant
    if P == []
        P = pop!(P_remainder)
        file_diff = open("log_diff.txt","a")
        println(file_diff,"changement manuel de P car il est vide dans subdivide_polygon :\nP = $P\nP_remainder = $P_remainder")
        close(file_diff)
        return P_remainder,P,[]
    end
    if maximum(angles) > pi
        non_convex = true
    else
        non_convex = false
    end
    file_diff = open("log_diff.txt","a")
    println(file_diff,"non_convex = $non_convex\tangles = $angles")
    close(file_diff)
    if !non_convex
        #println("convex")
        #println("angles : ", angles)
        #println(minimum(angles)-angle_seuil)
        if minimum(angles) > angle_seuil
            # cas 1 : coupe entre le sommet s1 de plus petit angle et le sommet s2 qui donne
            # deux angles en s1 de max le plus petit possible
            #=file_diff = open("log_diff.txt","a")
            println(file_diff,"subdivision")
            close(file_diff)=#
            s1 = argmin(angles)
            future_angles = []
            s1m1 = mod(s1-2,n)+1
            s1p1 = mod(s1,n)+1
            # b pour before, a pour after
            xb2 = P[s1m1][1]-P[s1][1]
            yb2 = P[s1m1][2]-P[s1][2]
            xa1 = P[s1p1][1]-P[s1][1]
            ya1 = P[s1p1][2]-P[s1][2]
            for i in 1:(n-3)
                num = mod(s1+i,n)+1
                xb1 = P[num][1]-P[s1][1]
                yb1 = P[num][2]-P[s1][2]
                xa2 = xb1
                ya2 = yb1
                angle_before = eval_angle(xb1,yb1,xb2,yb2)
                angle_after = eval_angle(xa1,ya1,xa2,ya2)
                push!(future_angles,max(angle_before,angle_after))
            end
            s2 = mod(s1+argmin(future_angles),n)+1
            #println("s1 : $s1\ts2 : $s2\tpolygone $P")
            # ajout d'une ligne en rouge dans le plot algorithm_plot pour visualiser les subdivisions
            X1 = P[s1]
            X2 = P[s2]
            push!(bonus_plot_infos, [[X1[1],X2[1]], [X1[2],X2[2]]])
            args["bonus_plot_infos"] = bonus_plot_infos
            # construction des deux polygones P1 et P2 à partir de P
            if s1 < s2
                P1 = P[1:s1]
                for i in s2:n
                    push!(P1,P[i])
                end
                P2 = P[s1:s2]
            else
                P1 = P[1:s2]
                for i in s1:n
                    push!(P1,P[i])
                end
                P2 = P[s2:s1]
            end
            #println("P1 :\n$P1\nP2 :\n$P2")
            ##### ajout de P2 dans P_remainder et remplacement de P par P1
            # ajout de P1 et P2 dans last_iteration_candidate
            if last_iter_test
                last_iteration_candidate = [P1,P2]
                P = []
            else
                push!(P_remainder,P2)
                P = P1
            end
        end
    else
        error("non convex: $angles\n$P")
        # cas 2 : coupe entre le sommet qui n'est pas un point extrême et le plus grand angle non voisin
    end
    return P_remainder,P,last_iteration_candidate
end

function polygon_intersection(polygon1, polygon2)
    # retourne l'intersection de polygon1 et polygon2 qui sont convexes, sinon l'intersection ne marche pas je crois
    magnitude1,precision1 = find_magnitude_and_precision(polygon1)
    magnitude2,precision2 = find_magnitude_and_precision(polygon2)
    magnitude = max(magnitude1,magnitude2)
    precision = max(precision1,precision2)
    poly1 = create_polygon(polygon1, magnitude, precision)
    poly2 = create_polygon(polygon2, magnitude, precision)
    c = Clip()
    add_path!(c, poly1, PolyTypeSubject, true)
    add_path!(c, poly2, PolyTypeClip, true)
    result, polys = execute(c, ClipTypeIntersection, PolyFillTypeEvenOdd, PolyFillTypeEvenOdd)
    # dernière étape : on reconvertit polys en une liste de polygone (eux-mêmes de type liste) en prenant en compte magnitude et precision
    res = []
    for i in 1:length(polys)
        poly = polys[i]
        pol = []
        for j in 1:length(poly)
            push!(pol, [Float64(poly[j].X)/10^(precision),Float64(poly[j].Y)/10^(precision)])
        end
        push!(res,pol)
    end
    return res
end

function test_polygon_intersection(polygon1,polygon2)
    # retourne true si intersection entre polygon1 et polygon2, false sinon
    intersection = polygon_intersection(polygon1,polygon2)
    if intersection == Any[]
        return false
    end
    return true
end
