using StringEncodings

include("print_triangulation.jl")
include("plot_polygons.jl")

SEP1 = "::"
SEP2 = " "
SEP3 = "_"
SEP_MAT = ["___","__","_"]

function write_triangulation_to_txt(T_set, sep_mat)
    # retourne un string contenant T_set au format txt avec sep_mat[1] pour
    # délimiteur interpolygone, sep_mat[2] pour délimiteur de sommet d'un même polygone
    # et sep_mat[3] pour délimiteur de coordonnées d'un sommet
    sep1 = sep_mat[1]
    sep2 = sep_mat[2]
    sep3 = sep_mat[3]
    s = ""
    for i in 1:length(T_set)
        T = T_set[i]
        for j in 1:length(T)
            S = T[j]
            for k in 1:length(S)
                s = string(s,"$(S[k])")
                if k < length(S)
                    s = string(s,sep3)
                end
            end
            if j < length(T)
                s = string(s,sep2)
            end
        end
        if i < length(T_set)
            s = string(s,sep1)
        end
    end
    return s
end

function parse_triangulation_from_txt(s,sep_mat)
    # parse une triangulation mise en txt grâce à write_triangulation_to_txt
    T_set = Vector{Vector{Float64}}[]
    infos_triangles = split(s,sep_mat[1])
    for i in 1:length(infos_triangles)
        s_triangle = infos_triangles[i]
        T = Vector{Float64}[]
        infos_vertices = split(s_triangle,sep_mat[2])
        for j in 1:length(infos_vertices)
            s_vertex = infos_vertices[j]
            S = Float64[]
            infos_vertex = split(s_vertex,sep_mat[3])
            for k in 1:length(infos_vertex)
                push!(S,parse(Float64,infos_vertex[k]))
            end
            push!(T,S)
        end
        push!(T_set,T)
    end
    return T_set
end

function save_triangle_plus_rectangle_result(filename,str_exprf,limits,err,size_Z1,size_Z2,choices,t,T_set,name_algo="triangle+rectangle")
    # écrit au format txt les informations utiles du run de l'algo triangle+rectangle
    s = ""
    # récupération des délimiteurs
    sep1 = SEP1
    sep2 = SEP2
    sep3 = SEP3
    sep_mat = SEP_MAT

    # début de l'écriture du string résultat dans s
    # ajout de la signature de l'algo
    s = string(s,name_algo,sep1)
    # ajout de l'expression de la fonction
    s = string(s,str_exprf,sep1)
    # ajout des limites du domaine
    A,B,C,D = limits
    s = string(s,"$A",sep2,"$B",sep2,"$C",sep2,"$D",sep1)
    # ajout de l'erreur
    if typeof(err) == Absolute
        delta = err.delta
        eps = 0
    elseif typeof(err) == Relative
        eps = err.percent/100
        delta = 0
    end
    s = string(s,"delta$(Float64(delta))",sep2,"epsilon$(Float64(eps))",sep1)
    # ajout du nombre de polygones dans la triangulation
    s = string(s,"$(length(T_set))",sep2,"polygones",sep1)
    # t
    s = string(s,t,sep2,"secondes",sep1)
    # ajout des paramètres restants
    # size_Z1
    s = string(s,"size_Z1",sep2,"$(size_Z1[1])",sep3,"$(size_Z1[2])",sep2)
    # size_Z2
    s = string(s,"size_Z2",sep2,"$(size_Z2[1])",sep3,"$(size_Z2[2])",sep2)
    # choices
    for i in 1:length(choices)
        println(choices[i])
        s = string(s,"$(choices[i][1])",sep2,"$(choices[i][2])")
        if i < length(choices)
            s = string(s,sep2)
        end
    end
    s = string(s,sep1)
    # ajout de la triangulation
    s_triangulation = write_triangulation_to_txt(T_set,sep_mat)
    s = string(s,s_triangulation)

    # écriture de s dans le fichier filename
    file = open(filename,"a")
    println(file,s) # pas besoin de mettre de \n dans s parce que println le fait (et pas write)
    close(file)
    return s
end

function parse_all_heuristics_result(s)
    # parse une ligne résultat de triangle+rectangle

    # récupération des délimiteurs
    sep1 = SEP1
    sep2 = SEP2
    sep3 = SEP3
    sep_mat = SEP_MAT

    # parsing
    elements = split(s,sep1)

    method = elements[1]

    str_exprf = elements[2]

    s_limits = split(elements[3],sep2)
    limits = [parse(Float64,s_limits[1]),parse(Float64,s_limits[2]),parse(Float64,s_limits[3]),parse(Float64,s_limits[4])]

    s_err = split(elements[4],sep2)
    delta,epsilon = [parse(Float64,s_err[1][6:end]),parse(Float64,s_err[2][8:end])]
    if epsilon == 0.0
        err = Absolute(delta)
    elseif delta == 0.0
        err = Relative(epsilon/100)
    else
        error("parsage de l'erreur non géré : $(elements[4])")
    end

    s_n = split(elements[5],sep2)
    n_triangle = parse(Int64, s_n[1])

    s_t = split(elements[6],sep2)
    t = 0
    try
        t = parse(Float64, s_t[1])
    catch e
        t = s_t[1]
    end
    elements[7] = replace(elements[7], ", " => ",")
    s_choices = split(elements[7],sep2)
    choices = []
    #println(s_choices)
    if s_choices[end] == ""
        s_choices = s_choices[1:end-1]
    end
    for i in 1:Int(length(s_choices)/2)
        param = [s_choices[2*i-1], s_choices[2*i]]
        param[2] = replace(param[2], " " => "")
        push!(choices, param)
    end

    if elements[8] != ""
        T_set = parse_triangulation_from_txt(elements[8],sep_mat)
    else
        T_set = [[]]
    end

    return str_exprf, limits, err, n_triangle, t, choices, T_set
end

function get_all_parsed_results(filename = "/home/aduguet/Documents/doctorat/2dpwlb/codes/julia/heuristique triangle+rectangle/resultats/resultats.txt")
    # retourne une liste de string dont chacun décrit une sortie de l'algo triangle+rectangle
    file = open(filename,"r")
    strings = read(file)
    close(file)
    strings = decode(Vector{UInt8}(strings), "utf-8")
    s = split(strings,"\n")
    s = s[1:end]
    if s[end] == ""
        s = s[1:end-1] # ne pas toucher, supprime la dernière ligne vide
    end

    # parsing de tous les résultats
    parsed = []
    for i in 1:length(s)
        if s[i] != ""
            push!(parsed, parse_all_heuristics_result(s[i]))
        end
    end
    return parsed
end

function sum_up_results(parsed, filename = "../tableau_resultats/dernier_tableau.txt", list_exit = false)
    # prend la sortie de get_all_parsed_results et retourne un tableau :
    # (str_exprf, delta, solved,time,number_of_polygones) (expression de la fonction, erreur, booléen valant true si bien résolu)
    # le quatrième contient les temps d'exécutions (Inf si non résolu) à la seconde près
    # le cinquième le nombre de polygones dans la triangulation (Inf si non résolu)
    # list_exit == true renvoie les infos sous forme de list sinon, renvoie sous forme de string

    # écriture des infos pertinentes
    infos = []
    n = length(parsed)
    for i in 1:n
        current_parsed = parsed[i]
        info = []
        # entrée de la fonction et du delta
        push!(info, current_parsed[1])
        push!(info, current_parsed[3])
        # entrée de solved, avec temps = Inf et nombre de polygones = Inf si false
        # solved = true ssi le nombre de polygones est différent de 0 et sinon les 40 premiers caractères du message d'erreur
        if current_parsed[4] == 0.0
            push!(info, current_parsed[5][1:min(end,40)])
            push!(info, Inf)
            push!(info, Inf)
        else
            push!(info, true)
            push!(info, current_parsed[5])
            push!(info, current_parsed[4])
        end
        push!(infos,info)
    end

    if list_exit
        return infos
    end

    # tri de infos pour mettre les expés sur la même fonction à côté
    sort!(infos, by = x -> x[1])

    # affichage des infos pertinentes
    # nombres de caractères par colonne : 24 pour fonction, 18 pour delta, 9 pour # polygones, 6 pour le temps et 15 pour l'erreur éventuelle
    s = string("fonction"," "^16,"delta"," "^13,"# polys"," "^2,"temps"," ","solved","\n")
    for i in 1:n
        info = infos[i]
        n_expr = length(info[1])
        n_delta = length(string(info[2]))
        n_solved = length(string(info[3]))
        n_time = length(string(info[4]))
        n_poly = length(string(info[5]))
        println(n_expr," ",n_delta," ",n_solved," ",n_time," ",n_poly)
        line = string(info[1]," "^(24-n_expr),info[2]," "^(18-n_delta),info[5]," "^(9-n_poly),info[4]," "^(6-n_time),info[3],"\n")
        s = string(s,line)
    end


    file = open(filename, "a")
    println(file,s)
    close(file)
    return s
end

function print_parsed_result(res, aff = true, aff_poly = true) # pas adapté à flex_heuristic, mais à triangle_plus_rectangle
    # affiche de manière lisible les paramètres et la solution d'une instance
    str_exprf, limits, err, n_triangle, t, size_Z1, size_Z2, choices, T_set = res
    s = ""
    s = string(s, "Fonction : $str_exprf\n")
    A,B,C,D = limits
    s = string(s, "Domaine : [$A,$B]x[$C,$D]\n")
    s = string(s, "Erreur : $err\n")
    s = string(s, "---> $n_triangle polygones en $t secondes\n")
    s = string(s, "Informations complémentaires :\ngrille sur Z1 : [$(Int(size_Z1[1])),$(Int(size_Z1[2]))]\ngrille sur Z2 : [$(Int(size_Z2[1])),$(Int(size_Z2[2]))]\n")
    for i in 1:length(choices)
        name,val = choices[i]
        s = string(s, "$name : $val\n")
    end
    s = string(s, "Triangulation :\n")
    for i in 1:length(T_set)
        s_poly = print_polygon(T_set[i])
        s = string(s, "Triangle $i : $s_poly\n")
    end
    s = string(s,"\n")
    if aff
        println(s)
    end
    if aff_poly
        title = "$str_exprf d'erreur $(Float64(delta)) : $n_triangle polygones en $(round(t,digits=2))s"
        p = plot(title=title,legend=false)
        for i in 1:length(T_set)
            plot_polygon(p,T_set[i],:green,true)
        end
        display(p)
    end
    return s
end
