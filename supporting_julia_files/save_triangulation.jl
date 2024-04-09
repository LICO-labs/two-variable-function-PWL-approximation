include("tools_save_triangulation.jl")

function save_triangulation_result(filename,str_exprf,domain,err,arguments,t = "UNK",polygons = [],name_algo = "UNK")
    # écrit au format txt les informations utiles d'un run dans le fichier filename
    # pour construire une pwl dans un corridor défini par str_exprf et err, sur un domaine domain
    # avec les arguments complémentaires arguments, un temps d'exécution t, une triangulation polygons
    # et pour l'algo name_algo
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
    A,B,C,D = domain
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
    s = string(s,"$(length(polygons))",sep2,"polygones",sep1)
    # t
    s = string(s,t,sep2,"secondes",sep1)
    # ajout des paramètres restants
    # arguments
    for key in keys(arguments)
        if key != "bonus_plot_infos"
            name = key
            value = arguments[key]
            s = string(s,"$name",sep2,"$value")
            s = string(s,sep2)
        end
    end
    s = string(s,sep1)
    # ajout de la triangulation
    s_triangulation = write_triangulation_to_txt(polygons,sep_mat)
    s = string(s,s_triangulation)

    # écriture de s dans le fichier filename
    file = open(filename,"a")
    println(file,s) # pas besoin de mettre de \n dans s parce que println le fait (et pas write)
    close(file)
    return s
end
