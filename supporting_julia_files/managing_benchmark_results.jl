using Colors, LaTeXStrings

include("basic_functions.jl")
include("corridor.jl")
include("tools_save_triangulation.jl")
include("functions.jl")

solution_filename = "../resultats_numeriques/SDI/hybrid.txt"
method_name = "SDIEff95Hybrid"
filename_table = "../resultats_numeriques/tableau_recap_par_heuristique_comparaison_incl.txt"
filename_time_table = "../resultats_numeriques/tableau_recap_temps_par_heuristique_comparaison_incl.txt"

vals_reb2D = [16,20,48,80,224,24,28,84,121,351,4,12,20,59,94,2,6,10,31,350,5,8,16,44,85,2,4,9,23,40,6,6,21,96,274,6,9,12,40,87,2,4,6,84,86]
vals_reb1D = [6,8,15,32,60,6,8,15,32,60,6,8,18,45,98,4,9,12,30,270,12,16,24,72,125,3,10,28,27,48,20,42,80,168,340,6,7,18,40,83,3,4,7,18,34]
vals_kazda = [6,6,Inf,Inf,Inf,Inf,6,Inf,Inf,Inf,4,6,12,Inf,Inf,1,2,4,Inf,Inf,Inf,Inf,8,Inf,Inf,2,4,6,12,Inf,1,4,5,Inf,Inf,3,4,6,Inf,Inf,1,2,4,5,Inf]
vals_LinA1D = [6,8,15,21,60,6,8,15,21,60,6,8,18,45,91,4,9,12,30,238,12,16,22,68,120,3,8,21,27,44,20,42,72,168,340,9,9,16,49,81,4,4,9,16,36]
vals_KL23LP = [6,10,18,36,98,7,11,20,Inf,Inf,6,10,16,36,64,1,3,6,19,Inf,2,6,10,15,36,2,4,6,14,25,1,4,13,24,43,3,4,14,22,50,1,2,19,19,19]
vals_KL23relaxed = [6,8,16,32,83,6,7,14,25,Inf,4,10,14,30,57,1,4,4,18,181,3,3,8,19,27,2,4,6,14,26,1,4,6,21,44,3,4,10,16,38,1,2,12,8,14]
times_reb2D = [30.8,84.4,150.4,272.6,380.6,26.8,7.4,38.0,35.8,171.7,0.8,72.4,4.7,59.3,45.3,0.3,18.7,12.7,54.6,652.6,1.0,13.1,30.0,74.6,141.9,1.8,1.0,25.8,14.4,161.4,7.7,1.3,30.8,73.0,305.5,22.8,15.6,22.9,202.8,174.1,0.8,66.6,4.4,231.5,57.8]
times_reb1D = [0.5,0.4,0.5,1.2,1.6,0.5,0.4,0.5,1.2,1.6,0.4,0.5,0.6,1.0,1.3,0.4,0.6,0.7,0.9,2.7,0.6,1.2,1.3,1.7,2.4,0.5,0.6,1.4,1.0,2.6,1.2,1.4,1.3,2.6,3.4,0.5,0.5,0.7,0.8,1.0,0.3,0.3,0.4,0.6,0.7]
times_kazda = [34+53.1,88.84+60.03,Inf,Inf,Inf,Inf,1313.49+181.97,Inf,Inf,Inf,2.22+16.26,297.79+119.01,14565.71+196.79,Inf,Inf,0.76+5.20,1.25+12.16,6.02+32.98,Inf,Inf,Inf,Inf,166.63+175.93,Inf,Inf,0.6+12.7,1.2+18.39,11.11+55.79,6928.58+317.92,Inf,0.42+5.88,23.58+63.35,1880.45+210.8,Inf,Inf,3.53+19.73,14.66+112.58,57.96+247.53,Inf,Inf,0.29+2.62,0.66+8.84,5.47+30.02,222.11+155.57,Inf]
times_LinA1D = [0.014,0.004,0.004,0.004,0.006,0.003,0.003,0.003,0.004,0.004,0.003,0.004,0.004,0.006,0.006,0.004,0.006,0.008,0.005,0.008,0.004,0.007,0.007,0.008,0.013,0.004,0.007,0.007,0.006,0.006,0.006,0.007,0.007,0.01,0.012,0.267,0.004,0.005,0.006,0.006,0.003,0.003,0.004,0.004,0.004]
# the two following counts only MILP solving time
times_KL23LP = [0.97,17,51,2197,8454,6,392,11449,Inf,Inf,0.20,8.3,152,1912,30824,0.03,0.06,0.34,244,Inf,0.05,3.9,11,45,3514,0.03,0.08,2.8,106,267,0.02,0.46,37,2031,35999,0.22,1.2,62,115,9448,0.01,0.05,40,285,178]
times_KL23relaxed = [0.29,0.80,3,23,156,30,1.7,45,361,Inf,0.10,1.6,6.0,46,88,0.04,0.14,0.43,6.2,8418,0.05,0.11,0.31,8.0,37,0.04,0.12,0.67,3.4,12,0.03,0.52,1.1,44,130,0.29,0.43,11,17,80,0.004,0.07,2.9,1.8,8.0]

expr_functions = []
for i in 1:9
    push!(expr_functions, get_limits_and_func(i)[6])
end
deltas1 = [1.5,1.0,0.5,0.25,0.1]
deltas2 = deltas1
deltas3 = [1.0,0.5,0.25,0.1,0.05]
deltas4 = [0.1,0.05,0.03,0.01,0.001]
deltas5 = deltas3
deltas6 = [0.5,0.25,0.1,0.05,0.03]
deltas7 = deltas3
deltas8 = deltas3
deltas9 = deltas3
deltas = [deltas1,deltas2,deltas3,deltas4,deltas5,deltas6,deltas7,deltas8,deltas9]

function plot_benchmark_result(well_solved_mat, times_mat, number_of_polygons_mat, names, foldername, scatter = false)
    # sauvegarde un plot avec le nombre de polygones et un avec les temps pour les résultats donnés

    n_exp = length(names)
    colors = [:blue,:red,:green,:purple]

    # résultats de Rebennack
    reb_polygons = zeros(5,9)
    reb_times = zeros(5,9)
    reb_polygons[:,1] = [16,20,48,80,224]
    reb_polygons[:,2] = [24,28,84,121,351]
    reb_polygons[:,3] = [4,12,20,59,94]
    reb_polygons[:,4] = [2,6,10,31,350]
    reb_polygons[:,5] = [5,8,16,44,85]
    reb_polygons[:,6] = [2,4,9,23,40]
    reb_polygons[:,7] = [6,6,21,96,274]
    reb_polygons[:,8] = [6,9,12,40,87]
    reb_polygons[:,9] = [2,4,6,84,86]
    reb_times[:,1] = [31,84,150,273,381]
    reb_times[:,2] = [27,7,38,36,172]
    reb_times[:,3] = [1,72,5,59,45]
    reb_times[:,4] = [0,19,13,55,653]
    reb_times[:,5] = [1,13,30,75,142]
    reb_times[:,6] = [2,1,26,14,161]
    reb_times[:,7] = [8,1,31,73,306]
    reb_times[:,8] = [23,16,23,203,174]
    reb_times[:,9] = [1,67,4,232,58]

    # préparation des tableaux de la bonne taille
    well_solved = zeros(n_exp,5,9)
    times = zeros(n_exp,5,9)
    n_polygons = zeros(n_exp,5,9)

    # remplissage des tableaux
    for i in 1:2
        for j in 1:2
            well_solved[(i-1)*2+j,:,:] = well_solved_mat[i,j,:,:]
            times[(i-1)*2+j,:,:] = times_mat[i,j,:,:]
            n_polygons[(i-1)*2+j,:,:] = number_of_polygons_mat[i,j,:,:]
        end
    end

    # gestion des instances non résolues en mettant Inf en temps et en nombre de polygones
    times[well_solved.==false] .= Inf
    n_polygons[well_solved.==false] .= Inf

    # préparation des plots
    name_time = "Time comparison on Rebennack instances"
    name_number = "Comparison of the number of polygons on Rebennack instances"
    p_time = plot()
    p_number = plot()
    p_time = plot!(p_time,title=name_time)
    p_number = plot!(p_number,title=name_number)

    if !scatter
        # remplissage des plots
        for i in 1:4
            p_time = plot!(p_time,[k for k in 1:5], times[i,:,1], label="$(names[i])",color=colors[i])
            p_number = plot!(p_number,[k for k in 1:5], n_polygons[i,:,1], label="$(names[i])",color=colors[i])
            for j in 2:9
                p_time = plot!(p_time,[k+5*(j-1) for k in 1:5], times[i,:,j],label=false,color=colors[i])
                p_number = plot!(p_number,[k+5*(j-1) for k in 1:5], n_polygons[i,:,j],label=false,color=colors[i])
            end
        end
    else
        # remplissage des plots
        for i in 1:4
            p_time = scatter!(p_time,[k for k in 1:5], times[i,:,1], label="$(names[i])",color=colors[i])
            p_number = scatter!(p_number,[k for k in 1:5], n_polygons[i,:,1], label="$(names[i])",color=colors[i])
            for j in 2:9
                p_time = scatter!(p_time,[k+5*(j-1) for k in 1:5], times[i,:,j],label=false,color=colors[i])
                p_number = scatter!(p_number,[k+5*(j-1) for k in 1:5], n_polygons[i,:,j],label=false,color=colors[i])
            end
        end
    end

    # sauvegarde des plots
    plot_name_time = string(foldername,"/",name_time,"_scatter_",scatter)
    plot_name_number = string(foldername,"/",name_number,"_scatter_",scatter)
    for i in 1:length(names)
        plot_name_time = string(plot_name_time,"_",names[i])
        plot_name_number = string(plot_name_number,"_",names[i])
    end
    plot_name_time = string(plot_name_time, ".pdf")
    plot_name_number = string(plot_name_number, ".pdf")
    savefig(p_time, plot_name_time)
    savefig(p_number, plot_name_number)

    # affichage
    display(p_time)
    display(p_number)

    return well_solved,times,n_polygons
end

function add_column_to_comparison_table(col_name, vals, filename, char_per_info = 15)
    # ajoute la colonne de nom col_name et de valeurs vals dans le fichier filename en ajoutant juste l'info minimum
    # à la fin des length(vals)+1 lignes concernées, puis sauvegarde dans le même fichier

    # erreur si col_name ou vals[i] sont trop longs : il ne doivent pas dépasser char_per_info-2 caractères
    if length(col_name) > char_per_info-2 && length(col_name) <= char_per_info+1
        char_per_info += 2
    elseif length(col_name) > char_per_info+1
        error("col_name est plus long que char_per_info+1 qui vaut $(char_per_info+1), diminuer sa longueur puis relancer")
    end
    for i in 1:length(vals)
        if length(string(vals[i])) > char_per_info-2
            error("la valeur $i de vals est trop longue, le max autorisé est char_per_info-2 qui vaut $(char_per_info-2), diminuer sa longueur puis relancer")
        end
    end

    # récupération des lignes
    file = open(filename, "r")
    lines = readlines(file)
    close(file)

    println("lines :\n$lines")

    # ajout aux lignes
    # le nom de la colonne
    println("col_name $col_name\nlines[1] $(lines[1])")
    lines[1] = string(lines[1], col_name, " "^(char_per_info-length(col_name)))
    # chaque valeur
    for i in 2:(length(vals)+1)
        lines[i] = string(lines[i], vals[i-1], " "^(char_per_info-length(string(vals[i-1]))))
    end

    # enregistrement
    file = open(filename, "w")
    for i in 1:length(lines)
        println(file, lines[i])
    end
    close(file)

    return lines
end

function add_column_complete(solution_filename, method_name, filename_table = "../resultats_numeriques/tableau_recap_par_heuristique.txt", filename_time_table = "../resultats_numeriques/tableau_recap_temps_par_heuristique.txt")
    # rajoute une colonne à filename_table en lui donnant solution_filename le fichier txt contenant les lignes résultats
    # avec method_name le nom de la méthode

    parsed = get_all_parsed_results(solution_filename)
    infos = sum_up_results(parsed,"",true)

    # récupération du nombre de polygones dans vals
    vals = [infos[i][5] for i in 1:length(infos)]

    # récupération du temps dans vals_time
    vals_time = [infos[i][4] for i in 1:length(infos)]
    for i in 1:length(infos)
        if vals_time[i] != Inf
            vals_time[i] = Int(round(vals_time[i], digits=0))
        end
    end

    vals2 = []
    vals2_time = []
    for i in 1:9
        for j in 1:5
            push!(vals2, vals[(j-1)*9+i])
            push!(vals2_time, vals_time[(j-1)*9+i])
        end
    end
    vals = vals2
    vals_time = vals2_time

    s1 = add_column_to_comparison_table(method_name,vals,filename_table)
    s1 = add_column_to_comparison_table(method_name,vals_time,filename_time_table)

    return s1
end

function convert_name_exps_to_fr(name_exps)
    for i in 1:length(name_exps)
        name_exps[i] = replace(name_exps[i], "bd"=>"DB")
        name_exps[i] = replace(name_exps[i], "med"=>"ProMA")
        name_exps[i] = replace(name_exps[i], "Area"=>"Aire")
        name_exps[i] = replace(name_exps[i], "PaD"=>"VaT")
        name_exps[i] = replace(name_exps[i], "eff"=>"")
        name_exps[i] = replace(name_exps[i], "_true"=>"")
        name_exps[i] = replace(name_exps[i], "near60"=>"pp60")
    end
    return name_exps
end

function get_best_plus_x_bar_plots(results, name_exps, title = "comparison of algorithms", type_of_plot = "plot", save_figure = false, nice_leg = false, fr_name = true)
    # construit un histogramme best+x pour chacune des colonnes de résultats fournies dans results
    # name_exps indique le nom de l'expérience

    if fr_name
        name_exps = convert_name_exps_to_fr(name_exps)
    end

    n_exp = length(results)
    println(results)
    n_inst = length(results[1])

    # récupérer les minimums par instance
    best_results = []
    for i in 1:n_inst
        l_scores = []
        for j in 1:n_exp
            push!(l_scores, results[j][i])
        end
        push!(best_results, minimum(l_scores))
    end

    # calculer les pourcentages d'écarts au meilleur résultat
    diff_results = []
    for i in 1:n_exp
        l_diff = []
        for j in 1:n_inst
            push!(l_diff, (results[i][j]/best_results[j])-1)
        end
        push!(diff_results, l_diff)
    end

    # calculer le nombres d'instances par colonne de l'histogramme
    cols_threshold = [0,2,5,10,15,20,25,30,40,50,100,200]
    cols_threshold = [0,2,5,10,15,20,30,40,50,100,200]
    col_values = []
    for i in 1:n_exp
        l_values = []
        for threshold in cols_threshold
            push!(l_values, sum((diff_results[i][j] <= threshold/100) for j in 1:n_inst)/n_inst*100)
        end
        push!(col_values, l_values)
    end

    if type_of_plot == "bar"
        # préparer les plots
        for i in 1:n_exp
            p = bar(col_values[i], xticks = (1:length(cols_threshold), cols_threshold), xlabel = "% extra linear pieces in comparison to the best known solution", ylabel = "% of instances", title=name_exps[i], legend = false)
            p = plot!([0,length(cols_threshold)],[100,100], color=:black)
            display(p)
        end
    elseif type_of_plot == "plot"
        # préparer le plot
        p = plot(legend=:bottomright)
        new_plot = true
        if new_plot
            p = plot!([0,log(1+cols_threshold[end])],[100,100], color=:black, label = "")
            p = plot!([0,log(1+cols_threshold[end])],[0,0], color=:black, label = "")
            if !nice_leg
                for i in 1:n_exp
                    p = plot!(log.(1 .+ cols_threshold), col_values[i], xticks = (log.(1 .+ cols_threshold), 1 .+ cols_threshold./100), lw = 3, xlabel = L"\rho", ylabel = "% instances", label=name_exps[i])
                end
            else
                # suppose that the first half of the experiments are of one component and the other half of the other and differentiate them with :dash and :solid
                # while using different colors
                if length(n_exp) > 14
                    rng = MersenneTwister(1000)
                    colors = compute_list_of_colors()
                    colors = shuffle(colors)
                else
                    colors = ["blue","red","yellow","grey","orange","green","violet","pink","gold","black","aliceblue","sienna","coral","lawngreen"]
                end
                n_expd2 = Int(n_exp/2)
                for i in 1:n_expd2
                    # xlabel = "% extra linear pieces in comparison to the best known solution", ylabel = "% of instances"
                    p = plot!(log.(1 .+ cols_threshold), col_values[i], xticks = (log.(1 .+ cols_threshold), 1 .+ cols_threshold./100), lw = 3, xlabel = L"\rho", ylabel = "% instances", label=name_exps[i], color = colors[i], ls = :dash)
                end
                for i in n_expd2+1:n_exp
                    p = plot!(log.(1 .+ cols_threshold), col_values[i], xticks = (log.(1 .+ cols_threshold), 1 .+ cols_threshold./100), lw = 3, xlabel = L"\rho", ylabel = "% instances", label=name_exps[i], color = colors[i-n_expd2], ls = :solid)
                end
            end
        else
            # old plot
            #=p = plot!([0,log(1+cols_threshold[end])],[100,100], color=:black, label = "")
            p = plot!([0,log(1+cols_threshold[end])],[0,0], color=:black, label = "")
            if !nice_leg
                for i in 1:n_exp
                    p = plot!(log.(1 .+ cols_threshold), col_values[i], xticks = (log.(1 .+ cols_threshold), cols_threshold), lw = 3, xlabel = "% de pièces en plus en comparaison du plus petit nombre", ylabel = "% d'instances", label=name_exps[i])
                end
            else
                # suppose that the first half of the experiments are of one component and the other half of the other and differentiate them with :dash and :solid
                # while using different colors
                if length(n_exp) > 14
                    rng = MersenneTwister(1000)
                    colors = compute_list_of_colors()
                    colors = shuffle(colors)
                else
                    colors = ["blue","red","yellow","grey","orange","green","violet","pink","gold","black","aliceblue","sienna","coral","lawngreen"]
                end
                n_expd2 = Int(n_exp/2)
                for i in 1:n_expd2
                    # xlabel = "% extra linear pieces in comparison to the best known solution", ylabel = "% of instances"
                    p = plot!(log.(1 .+ cols_threshold), col_values[i], xticks = (log.(1 .+ cols_threshold), cols_threshold), lw = 3, xlabel = "% de pièces en plus en comparaison du plus petit nombre", ylabel = "% d'instances", label=name_exps[i], color = colors[i], ls = :dash)
                end
                for i in n_expd2+1:n_exp
                    p = plot!(log.(1 .+ cols_threshold), col_values[i], xticks = (log.(1 .+ cols_threshold), cols_threshold), lw = 3, xlabel = "% de pièces en plus en comparaison du plus petit nombre", ylabel = "% d'instances", label=name_exps[i], color = colors[i-n_expd2], ls = :solid)
                end
            end=#
        end
        if save_figure
            savefig(p,"../resultats_numeriques/unclassified analysis/"*title*".pdf")
        end
        display(p)
    end

    return col_values
end

function select_instances(results, functions, deltas)
    # results est la liste de benchmark de Rebennack et on ne conserve que les résultats
    # pour lesquels la fonction est dans functions et le delta est dans deltas
    # ([1,3,4] veut dire on garde les deltas 1, 3 et 4 dans l'ordre décroissant par fonction)

    # prétraitement de functions en funcs
    funcs = zeros(9)
    for i in 1:9
        if i in functions
            funcs[i] = 1
        end
    end

    # prétraitement de deltas en dels
    dels = zeros(5)
    for i in 1:5
        if i in deltas
            dels[i] = 1
        end
    end

    # traitement de res avec funcs et dels
    exps = []
    for res in results
        exp = []
        for i in 1:9
            for j in 1:5
                if funcs[i]+dels[j] == 2
                    push!(exp, res[(i-1)*5+j])
                end
            end
        end
        push!(exps, exp)
    end

    return exps
end

function create_title_comparison(functions, deltas)
    # retourne un titre pour un get_best_plus_x_bar_plots avec les numéros des fonctions et deltas gardés

    s = "comparison on functions "
    for i in functions
        s = string(s, i)
    end
    s = string(s, " and deltas ")
    for i in deltas
        s = string(s, i)
    end

    return s
end

function get_piece_number_from_results(filename)
    # retourne un tableau contenant les nombres de morceaux des instances résolues dans filename
    # en les rangeant dans l'ordre des instances de Rebennack (delta, num_func)
    parsed = get_all_parsed_results(filename)
    infos = sum_up_results(parsed,"",true)
    vals = [infos[i][5] for i in 1:length(infos)]
    vals2 = []
    for i in 1:9
        for j in 1:5
            push!(vals2, vals[(j-1)*9+i])
        end
    end
    return vals2
end

function get_vals_from_benchmark(filename)
    # récupère le nombre de polygones de chaque instance du fichier filename du dossier resultats_numeriques
    if !occursin("../resultats_numeriques",filename)
        parsed = get_all_parsed_results("resultats_numeriques/"*filename) # "../" enlevé de la chaîne de caractère
    else
        parsed = get_all_parsed_results(filename)
    end
    infos = sum_up_results(parsed,"",true)
    vals = [infos[i][5] for i in 1:length(infos)]
    vals2 = []
    for i in 1:9
        for j in 1:5
            push!(vals2, vals[(j-1)*9+i])
        end
    end
    return vals2
end

function get_times_from_benchmark(filename)
    # récupère le temps d'exécution de chaque instance du fichier filename du dossier resultats_numeriques
    if !occursin("../resultats_numeriques",filename)
        parsed = get_all_parsed_results("resultats_numeriques/"*filename) # "../" enlevé de la chaîne de caractère
    else
        parsed = get_all_parsed_results(filename)
    end
    infos = sum_up_results(parsed,"",true)
    vals = [infos[i][4] for i in 1:length(infos)]
    vals2 = []
    for i in 1:9
        for j in 1:5
            push!(vals2, vals[(j-1)*9+i])
        end
    end
    return vals2
end

function get_all_vals_from_benchmark(filenames)
    # récupère le nombre de polygones pour chaque instance pour chaque série d'expés de filenames
    res_vals = []
    res_times = []
    for filename in filenames
        push!(res_vals,get_vals_from_benchmark(filename))
        push!(res_times,get_times_from_benchmark(filename))
    end
    return res_vals,res_times
end

function build_table(expr_functions, deltas, results_vals, results_times, choice = "val et temps", best_def = "min", only_time = false)
    # retourne un string permettant de plot proprement un tableau de résultats
    # expr_functions est un vecteur contenant l'expression des fonctions, deltas les delta par fonction (double liste)
    # results_vals les nombres de polygones par heuristiques et results_times les temps correspondant

    n_heur = length(results_vals)

    # round all times to 1 digits
    for i in 1:length(results_times)
        for j in 1:length(results_times[i])
            try
                if results_times[i][j] >= 100
                    results_times[i][j] = Int(round(results_times[i][j], digits=0))
                elseif results_times[i][j] >= 1
                    results_times[i][j] = round(results_times[i][j], digits=1)
                else
                    results_times[i][j] = round(results_times[i][j], digits=3)
                end
            catch e
                println("err catched at the start of build_table: $e")
            end
        end
    end

    # saves the number of best solutions
    card_best = [0 for i in 1:n_heur]
    s = ""
    for i in 1:9
        for j in 1:5
            exprf = expr_functions[i]
            delta = deltas[i][j]
            #s = string(s, "$exprf & $delta & ")
            s = string(s, "$i & $delta & ")

            # compute minimum among heuristics, or maximum, or nothing to put it in bold
            if best_def == "min"
                mini = minimum([results_vals[k][5*(i-1)+j] for k in 1:n_heur])
            elseif best_def == "max"
                mini = maximum([results_vals[k][5*(i-1)+j] for k in 1:n_heur])
            elseif best_def == "nothing"
                mini = Inf
            elseif best_def == "sparse-max"
                if maximum([results_vals[k][5*(i-1)+j] for k in 1:n_heur]) == minimum([results_vals[k][5*(i-1)+j] for k in 1:n_heur])
                    mini = Inf
                else
                    mini = maximum([results_vals[k][5*(i-1)+j] for k in 1:n_heur])
                end
            elseif best_def == "sparse-min"
                if maximum([results_vals[k][5*(i-1)+j] for k in 1:n_heur]) == minimum([results_vals[k][5*(i-1)+j] for k in 1:n_heur])
                    mini = Inf
                else
                    mini = minimum([results_vals[k][5*(i-1)+j] for k in 1:n_heur])
                end
            end

            #choice = "val et temps"
            for k in 1:n_heur
                val = results_vals[k][5*(i-1)+j]
                if val != Inf
                    if !only_time
                        val = Int(val)
                    else
                        # option to just make table of times, by giving them in argument results_vals
                        val = val
                    end
                end
                if val == mini
                    card_best[k] += 1
                    if choice == "val et temps"
                        if typeof(results_times[k][5*(i-1)+j]) <: AbstractFloat
                            digits = 1
                            if results_times[k][5*(i-1)+j] < 0.05
                                digits = 3
                            end
                            s = string(s, "{\\bf $(val)} & {\\color{gray} $(round(results_times[k][5*(i-1)+j],digits=digits))} & ")
                        elseif typeof(results_times[k][5*(i-1)+j]) <: Integer
                            s = string(s, "{\\bf $(val)} & {\\color{gray} $(results_times[k][5*(i-1)+j])} & ")
                        else
                            s = string(s, "{\\bf $(val)} & {\\color{gray} $(results_times[k][5*(i-1)+j])} & ")
                        end
                    elseif choice == "val"
                        s = string(s, "{\\bf $(val)} & ")
                    end
                else
                    if choice == "val et temps"
                        if typeof(results_times[k][5*(i-1)+j]) <: AbstractFloat
                            digits = 1
                            if results_times[k][5*(i-1)+j] < 0.05
                                digits = 3
                            end
                            s = string(s, "$(val) & {\\color{gray} $(round(results_times[k][5*(i-1)+j],digits=digits))} & ")
                        elseif typeof(results_times[k][5*(i-1)+j]) <: Integer
                            s = string(s, "$(val) & {\\color{gray} $(results_times[k][5*(i-1)+j])} & ")
                        else
                            s = string(s, "$(val) & {\\color{gray} $(results_times[k][5*(i-1)+j])} & ")
                        end
                    elseif choice == "val"
                        s = string(s, "$(val) & ")
                    end
                end
            end
            if j == 5
                s = string(s[1:(end-3)]," \\\\ \\hline \n")
            else
                s = string(s[1:(end-3)]," \\\\ \n")
            end
        end
    end
    new_expr_functions = ["\$x^2-y^2\$","\$x^2+y^2\$","\$xy\$","\$x exp(-x^2-y^2)\$","\$x sin(y)\$","\$\\frac{sin(x)}{x}y^2\$","\$x sin(x) sin(y)\$","\$(x^2-y^2)^2\$","\$exp(-10 (x^2-y^2)^2)\$"]
    for i in length(new_expr_functions):-1:1
        s = replace(s,expr_functions[i] => new_expr_functions[i])
    end
    s = replace(s,"Inf \$ {\\color{gray} Inf"=>"- \$ {\\color{gray} TO")
    println(s)

    # print cardinal of best solutions
    s_best = ""
    for k in 1:n_heur
        s_best = string(s_best, "$(Int(card_best[k])) & ")
    end
    s_best = s_best[1:end-3]
    println("-----\nnumber of best solutions:\n$s_best\n")

    return s
end

function get_all_filenames()
    # return a list with all filenames for comparison of my heuristics
    root = "../resultats_numeriques/exps_in_thesis/"
    folders = ["mean_along_edge/PaD/","mean_along_edge/polygon_area/","bisector_direction/PaD/","bisector_direction/polygon_area/","mean_along_edge_near60/PaD/","mean_along_edge_near60/polygon_area/","bisector_direction_near60/PaD/","bisector_direction_near60/polygon_area/"]

    filenames = []
    exps_names = []
    for folder in folders
        list_files = readdir(string(root,folder))
        for file in list_files
            #println(string(root,folder,file))
            split1 = split(folder,"/")[1]
            split2 = split(folder,"/")[2]
            if split1 == "mean_along_edge"
                FDH = "med_"
            elseif split1 == "bisector_direction"
                FDH = "bd_"
            elseif split1 == "mean_along_edge_near60"
                FDH = "near60_med_"
            elseif split1 == "bisector_direction_near60"
                FDH = "near60_bd_"
            else
                FDH = "bd_"
            end
            if split2 == ""
                split2 = "$(split1)"
            end
            # add eff95 to exps without any mention of the value for this parameter
            if !occursin("lip2",file) && !occursin("eff99",file) && !occursin("eff95",file) && !occursin("20", file) && !occursin("50", file) && !occursin("70", file) && !occursin("90", file)
                s_file = string("eff95",file[1:end-4])
            else
                s_file = file[1:end-4]
            end

            # change "20","50","70","90" to there form with eff before
            s_file = replace(s_file,"20"=>"eff20")
            s_file = replace(s_file,"50"=>"eff50")
            s_file = replace(s_file,"70"=>"eff70")
            s_file = replace(s_file,"90"=>"eff90")

            # change Incl_false, Incl_true et hybrid en _false, _true et _hybrid
            s_file = replace(s_file,"Incl_false"=>"_false")
            s_file = replace(s_file,"Incl_true"=>"_true")
            s_file = replace(s_file,"hybrid"=>"_hybrid")

            # change polygon_area en polygonarea
            split2 = replace(split2,"polygon_area"=>"Area")

            # change "lip2" en "lip"
            s_file = replace(s_file,"lip2"=>"lip")

            exp_name = string(FDH,split2,"_",s_file)

            push!(filenames,string(root,folder,file))
            push!(exps_names,exp_name)
        end
    end
    return filenames,exps_names
end

struct instance_results
    obj
    time
end

function get_results_in_dict(filenames, exps_names)
    # return a dictionary with keys the exps_names, and value the results associated to the exp name

    res = Dict()
    results_vals,results_times = get_all_vals_from_benchmark(filenames)
    for i in 1:length(exps_names)
        exp_results = []
        for j in 1:length(results_vals[i])
            push!(exp_results, instance_results(results_vals[i][j],results_times[i][j]))
        end
        res[exps_names[i]] = exp_results
    end
    return res
end

mutable struct score
    win1
    eq
    win2
    log_obj_ratio
    log_time_ratio
    number_of_instances
    obj_scores
    obj_times
end

function add_score(res,sc)
    # add score sc to score res
    res.win1 += sc.win1
    res.eq += sc.eq
    res.win2 += sc.win2
    res.log_obj_ratio += sc.log_obj_ratio
    res.log_time_ratio += sc.log_time_ratio
    res.number_of_instances += sc.number_of_instances
    append!(res.obj_scores[1],sc.obj_scores[1])
    append!(res.obj_scores[2],sc.obj_scores[2])
    append!(res.obj_times[1],sc.obj_times[1])
    append!(res.obj_times[2],sc.obj_times[2])
    return res
end

function print_score(sc)
    # print a structure score
    print("$(sc.win1)/$(sc.eq)/$(sc.win2) - la méthode 1 donne en moyenne $(exp(sc.log_obj_ratio/sc.number_of_instances)) fois plus de morceaux ")
    println("et prend $(exp(sc.log_time_ratio/sc.number_of_instances)) fois plus de temps sur $(sc.number_of_instances) instances")
    println("nombre d'instances non résolues pour chaque méthode : $(sum([i==Inf for i in sc.obj_scores[1]])) / $(sum([i==Inf for i in sc.obj_scores[2]]))")
end

function instance_competition(inst_res1,inst_res2)
    # return a score, that can be multi-objective, for results of instance result inst_res1 and instance result inst_res2
    # [winner (1 or 2),val1/val2,t1/t2]

    win1 = 0
    eq = 0
    win2 = 0
    if inst_res1.obj < inst_res2.obj
        win1 = 1
    elseif inst_res1.obj > inst_res2.obj
        win2 = 1
    else
        eq = 1
    end
    # if one of the instance is not solved, the ratios are put to one so as to not be counted
    if inst_res1.obj == Inf || inst_res2.obj == Inf
        obj_ratio = 1
        time_ratio = 1
    else
        obj_ratio = inst_res1.obj/inst_res2.obj
        time_ratio = inst_res1.time/inst_res2.time
    end
    return score(win1,eq,win2,log(obj_ratio),log(time_ratio),1,[[inst_res1.obj],[inst_res2.obj]],[[inst_res1.time],[inst_res2.time]])
end

function exp_competition(exp1,exp2)
    # return a score structure, with win1 and win2 equal to how many times an instance of exp1 had a better or equal objective for win1, and same for win2

    if length(exp1) != length(exp2)
        error("number of instances different for exp1 ($(length(exp1))) and exp2 ($(length(exp2)))")
    end

    res = score(0,0,0,0,0,0,[[], []],[[], []]) # initialization of the multiobjective score as a score structure
    for i in 1:length(exp1)
        sc = instance_competition(exp1[i],exp2[i])
        res = add_score(res,sc)
    end
    return res
end

function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end

function get_list_of_parameters(exp)
    # return a list of strings which are the parameters of exp described in name of exp
    return split(exp[1],"_")
end

function get_list_of_other_parameters(exps,param,target)
    # return the experience in exps such that its list of parameters is that of target, after dropping the parameter param
    for exp in exps
        l = get_list_of_parameters(exp)
        l = remove!(l,param)
        if l == target
            return exp
        end
    end
    println("no exp found matching the exps $exps with list of parameters $target and with param $param in list $([[i[1]] for i in exps])")
    return -1
end

function create_list_only_with_parameter(param,results)
    # return a list with only results corresponding to exp_name with param
    res = []
    for exp in results
        exp_name = exp[1]
        if occursin(param,exp_name)
            push!(res,exp)
        end
    end
    return res
end

function create_list_without_parameter(param,results)
    # return a list with only results corresponding to exp_name with param
    res = []
    for exp in results
        exp_name = exp[1]
        if !occursin(param,exp_name)
            push!(res,exp)
        end
    end
    return res
end

function parameter_competition(param1,param2,exps_names,results)
    # return a score s1-s2 corresponding to the number of times an exp with param1 beat its counterpart in param2 (s1) and s2 is the reverse
    # a match between two exps is done by exp_competition

    # make two lists with the exps for each parameter
    exps1 = []
    exps2 = []
    for exp in results
        exp_name = exp[1]
        if occursin(param1,exp_name)
            push!(exps1,exp)
        end
        if occursin(param2,exp_name)
            push!(exps2,exp)
        end
    end
    println("nombre d'expés dans exps1 : $(length(exps1))\nexps1 : $([[i[1]] for i in exps1])")
    println("nombre d'expés dans exps2 : $(length(exps2))\nexps2 : $([[i[1]] for i in exps2])")
    # launch competitions
    # WARNING : exps1 and exps2 are lists of pairs with first element the name of the experience and second the results
    # please do not use results here
    res = score(0,0,0,0,0,0,[[], []],[[], []]) # res is the global score
    for exp1 in exps1 # exp1 is a pair name-results
        # find opponent to exp1 with list of parameter
        target = get_list_of_parameters(exp1)
        target = remove!(target,param1)
        exp2 = get_list_of_other_parameters(exps2,param2,target)
        if exp2 != -1
            println("competition between $(exp1[1]) and $(exp2[1])")

            # make the competition
            sc = exp_competition(exp1[2],exp2[2])
            res = add_score(res,sc)
        end
    end

    # best+x plot
    get_best_plus_x_bar_plots(res.obj_scores, [param1,param2], "comparison of parameters $param1 and $param2", "plot", true)
    return res
end

function launch_all_competitions()
    # launch parameter_competition for all pairs of parameters

    filenames,exps_names = get_all_filenames()
    results = get_results_in_dict(filenames, exps_names)

    # remove experiences with parameters undefined in the article
    undefined_parameters = ["SPaD","hybrid","false","near60","random"]
    for param in undefined_parameters
        results = create_list_without_parameter(param,results)
    end
    # exps_names does not correspond to results anymore so we fix it
    exps_names = []
    for exp in results
        push!(exps_names,exp[1])
    end

    param_direction_of_progress = ["med","bisector_direction"]
    param_piece_evaluation = ["Area","PaD"]
    #param_piece_evaluation = ["polygonarea","PaD","SPaD","random","near60"]
    #param_pwl_inner_corridor = ["lip","eff95","eff99"]
    #param_inclination = ["false","true","hybrid"]

    #param_list = [param_direction_of_progress,param_piece_evaluation,param_pwl_inner_corridor,param_inclination]
    param_list = [param_direction_of_progress,param_piece_evaluation]

    for l in param_list
        for i1 in 1:length(l)-1
            for i2 in i1+1:length(l)
                parameter_competition(l[i1],l[i2],exps_names,results)
            end
        end
    end



    return 0
end

function rename_exps(exps_names)
    # use format described in my paper at ISCO 2022 to rename the experiences
    for i in 1:length(exps_names)
        exp_name = exps_names[i]
        splits = split(exp_name,"_")
        exps_names[i] = splits[2]*splits[3]*splits[1]
        exps_names[i] = replace(exps_names[i], "eff20"=>"20")
        exps_names[i] = replace(exps_names[i], "eff50"=>"50")
        exps_names[i] = replace(exps_names[i], "eff70"=>"70")
        exps_names[i] = replace(exps_names[i], "eff90"=>"90")
        exps_names[i] = replace(exps_names[i], "eff95"=>"95")
        exps_names[i] = replace(exps_names[i], "eff99"=>"99")
    end
    return exps_names
end

function launch_all_best_plus_x_plots()
    # launch all best_plus_x plots relevant for selecting my best heuristic with eff95 and with eff99

    filenames,exps_names = get_all_filenames()
    results = get_results_in_dict(filenames, exps_names)
    results = sort(results)

    # remove experiences with parameters undefined in the article
    undefined_parameters = ["SPaD","hybrid","false","near60","random","lip","eff20"]
    for param in undefined_parameters
        results = create_list_without_parameter(param,results)
    end
    # exps_names does not correspond to results anymore so we fix it
    exps_names = []
    for exp in results
        push!(exps_names,exp[1])
    end

    # declare useful arguments for build_table
    expr_functions = []
    for i in 1:9
        push!(expr_functions, get_limits_and_func(i)[6])
    end
    deltas1 = [1.5,1.0,0.5,0.25,0.1]
    deltas2 = deltas1
    deltas3 = [1.0,0.5,0.25,0.1,0.05]
    deltas4 = [0.1,0.05,0.03,0.01,0.001]
    deltas5 = deltas3
    deltas6 = [0.5,0.25,0.1,0.05,0.03]
    deltas7 = deltas3
    deltas8 = deltas3
    deltas9 = deltas3
    deltas = [deltas1,deltas2,deltas3,deltas4,deltas5,deltas6,deltas7,deltas8,deltas9]

    # create lists with specific parameters
    #spec_params = ["eff20","eff50","eff70","eff90","eff95","eff99"]
    #spec_params = ["eff95","eff99"]
    lists_res = []
    lists = []
    lists_times = []
    lists_names = []
    for k in 1:length(spec_params)
        param = spec_params[k]
        res = create_list_only_with_parameter(param,results)
        push!(lists_res,res)
        #res = sort(res,by = x->x[1],f = cmp)
        push!(lists,[])
        push!(lists_times,[])
        push!(lists_names,[])
        for i in 1:length(res)
            exp = res[i]
            push!(lists[k],[])
            push!(lists_times[k],[])
            push!(lists_names[k],exp[1])
            for j in 1:length(exp[2])
                push!(lists[k][i],exp[2][j].obj)
                push!(lists_times[k][i],exp[2][j].time)
            end
        end
    end

    # create best plus x plot for each list of lists
    for i in 1:length(lists)
        lists_names[i] = rename_exps(lists_names[i])
        get_best_plus_x_bar_plots(lists[i],lists_names[i],"comparison between heuristics with $(spec_params[i]) parameter", "plot", true)
    end


    # find best heuristic among the 8 defined
    lists = []
    lists_times = []
    lists_names = []
    for exp in results
        push!(lists,[exp[2][k].obj for k in 1:length(exp[2])])
        push!(lists_times,[exp[2][k].time for k in 1:length(exp[2])])
        push!(lists_names,exp[1])
    end
    lists_names = rename_exps(lists_names)

    functions = [[1,2],[3,4,5,6,7,8,9]]
    titles = ["comparison between heuristics with linearly separable functions","comparison between heuristics with nonlinearly separable functions"]
    for i in 1:length(functions)
        println("lists $i :\n$(lists[i])")
        exps = select_instances(lists, functions[i], [1,2,3,4,5])
        get_best_plus_x_bar_plots(exps,lists_names,titles[i],"plot",true)
    end

    println("table with exps_names $(lists_names)")
    tables = build_table(expr_functions, deltas, lists, lists_times)
    println("\n\n")

    tables = ["",""]
    for i in 1:length(lists)
        println("table $i with exps_names $(lists_names[i])")
        tables[i] = build_table(expr_functions, deltas, lists[i], lists_times[i])
        println("\n\n")
    end

    #=# additional best_plus_x plots without all instances
    functions = [3,4,5,6,7,8,9]
    delta = [3,4,5]
    for i in 1:length(lists)
        println("before selecting instances:\n$(lists[i][1])")
        exps = select_instances(lists[i], functions, delta)
        println("after selecting instances:\n$(exps[1])")
        get_best_plus_x_bar_plots(exps,lists_names[i],"comparison between heuristics with $(spec_params[i]) parameter with nonlinearly separable functions and smaller deltas", "plot", true)
    end=#

    # additional best_plus_x plots with LinA2D
    functions = [1,2]
    delta = [1,2,3,4,5]
    global vals_LinA1D, times_LinA1D
    for i in 1:length(lists)
        println("before selecting instances:\n$(lists[i][1])")
        exps = select_instances(lists[i], functions, delta)
        vals_LinA1D = select_instances([vals_LinA1D], functions, delta)[1]
        times_LinA1D = select_instances([times_LinA1D], functions, delta)[1]
        exp_lina = vals_LinA1D
        #=exp_lina = []
        for j in 1:length(vals_LinA1D)
            push!(exp_lina, instance_results(vals_LinA1D[j],times_LinA1D[j]))
        end=#
        push!(exps,exp_lina)
        push!(lists_names[i],"LinA2D")

        println("after selecting instances:\n$(exps[1])\t$(exps[5])")
        get_best_plus_x_bar_plots(exps,lists_names[i],"comparison between heuristics with $(spec_params[i]) parameter with linearly separable functions", "plot", true)
    end

    # launch pairwise competitions to select best heuristics
    for i in 1:length(lists)
        println("\ncompetition in group $(lists_names[i]):")
        for j in 1:length(lists[i])
            for k in j+1:length(lists[i])
                sc = exp_competition(lists_res[i][j][2],lists_res[i][k][2])
                println("exps used : $(lists_res[i][j][1]) and $(lists_res[i][k][1])")
                println("score for competition $(lists_names[i][j]) vs $(lists_names[i][k]):")
                print_score(sc)
            end
        end
    end

    # print total resolution time for each heuristic
    for i in 1:length(lists_times)
        exps_times = lists_times[i]
        for j in 1:length(exps_times)
            if i == 1
                println("temps d'exécution total de l'expérience $(lists_names[i][j]) : $(sum([lists_times[i][j][k] for k in 1:length(lists_times[i][j])]))")
            elseif i == 2
                println("temps d'exécution total de l'expérience $(lists_names[i][j]) : $(sum([lists_times[i][j][k] for k in 1:length(lists_times[i][j])-1]))")
            end
        end
    end

    return lists,lists_names
end

function launch_comparison_to_state_of_the_art()
    # launch build_table with the state of the art algorithm, PaD99med and PaD95med

    filenames,exps_names = get_all_filenames()
    results = get_results_in_dict(filenames, exps_names)
    results = sort(results)

    # remove experiences with parameters undefined in the article
    undefined_parameters = ["SPaD","hybrid","false","near60","random","lip","bisector_direction","Area"]
    for param in undefined_parameters
        results = create_list_without_parameter(param,results)
    end
    # exps_names does not correspond to results anymore so we fix it
    exps_names = []
    for exp in results
        push!(exps_names,exp[1])
    end

    println(exps_names)
    exps_names = ["PaD95med","PaD99med"]

    results_vals = []
    results_times = []
    true_exps_names = ["med_PaD_eff95_true", "med_PaD_eff99_true"]
    println(typeof(results))
    for i in 1:2
        push!(results_vals,[results[i][2][k].obj for k in 1:length(results[i][2])])
        push!(results_times,[results[i][2][k].time for k in 1:length(results[i][2])])
    end

    # last numbers of ISCO 2022! keep intact!
    sc = exp_competition(results[2][2],results[1][2])
    print_score(sc)

    global vals_reb2D,vals_LinA1D,vals_kazda,times_reb2D,times_LinA1D,times_kazda
    push!(exps_names,"Rebennack 2D")
    push!(exps_names,"LinA 1D")
    push!(exps_names,"Kazda")
    push!(results_vals,vals_reb2D)
    push!(results_vals,vals_LinA1D)
    push!(results_vals,vals_kazda)
    push!(results_times,times_reb2D)
    push!(results_times,times_LinA1D)
    push!(results_times,times_kazda)

    # declare useful arguments for build_table
    expr_functions = []
    for i in 1:9
        push!(expr_functions, get_limits_and_func(i)[6])
    end
    deltas1 = [1.5,1.0,0.5,0.25,0.1]
    deltas2 = deltas1
    deltas3 = [1.0,0.5,0.25,0.1,0.05]
    deltas4 = [0.1,0.05,0.03,0.01,0.001]
    deltas5 = deltas3
    deltas6 = [0.5,0.25,0.1,0.05,0.03]
    deltas7 = deltas3
    deltas8 = deltas3
    deltas9 = deltas3
    deltas = [deltas1,deltas2,deltas3,deltas4,deltas5,deltas6,deltas7,deltas8,deltas9]

    # launch build_table
    build_table(expr_functions, deltas, results_vals, results_times)

    return results,results_vals,results_times
end

function scatter_piece_numbers(results, names, savefile, ylabel = "nombre de pièces")
    # make a scatter plot of the values inside results (more than one heuristic)

    colors = ["blue","red","yellow","grey","orange","green","violet","pink","gold","black","aliceblue","sienna","coral","lawngreen"]

    p = plot(xlabel = "instances", ylabel = ylabel, yscale = :log)
    n = length(results[1])
    for i in 1:length(results)
        p = scatter!(1:n, results[i], color = colors[i], label = names[i], marker = :cross)
    end
    savefig(p, "../resultats_numeriques/exps_in_thesis_best_heuristic_analysis/scatter_plot/"*savefile)
end

function compute_proportion_best_obj(group_results, names)
    # return the proportion of best instances for comparable instances
    # group_results has a specific form: group_results = vector{Vector{Vector{Int64}}}
    # where the first vector contains groups of groups of heuristics
    # the second vector contains groups of comparable heuristics
    # the third vector contains objective values of instances of a common benchmark for a heuristic

    n = length(names)
    best = [0.0 for i in 1:n]
    total_number_of_runs = 0
    for k in 1:length(group_results) # k a group of heuristics
        res = group_results[k]
        total_number_of_runs += length(res[1]) # number of instances in the group k
        for i in 1:length(res) # i a heuristic
            for j in 1:length(res[i]) # res[i][j] a value
                if res[i][j] == minimum([res[ii][j] for ii in 1:length(res)])
                    best[i] += 1
                end
            end
        end
    end
    for i in 1:n
        best[i] /= total_number_of_runs
        println("for parameter $(names[i]): $(best[i]) best number")
    end
    return best
end


filename_exps = ["../resultats_numeriques/exps_in_thesis/mean_along_edge/PaD/50Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/mean_along_edge/PaD/70Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/mean_along_edge/PaD/90Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/mean_along_edge/PaD/Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/mean_along_edge/PaD/eff99Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/mean_along_edge/polygon_area/50Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/mean_along_edge/polygon_area/70Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/mean_along_edge/polygon_area/90Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/mean_along_edge/polygon_area/Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/mean_along_edge/polygon_area/eff99Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/bisector_direction/PaD/50Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/bisector_direction/PaD/70Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/bisector_direction/PaD/90Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/bisector_direction/PaD/Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/bisector_direction/PaD/eff99Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/bisector_direction/polygon_area/50Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/bisector_direction/polygon_area/70Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/bisector_direction/polygon_area/90Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/bisector_direction/polygon_area/Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/bisector_direction/polygon_area/eff99Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/mean_along_edge_near60/PaD/Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/mean_along_edge_near60/PaD/eff99Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/mean_along_edge_near60/polygon_area/Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/mean_along_edge_near60/polygon_area/eff99Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/bisector_direction_near60/PaD/Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/bisector_direction_near60/PaD/eff99Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/bisector_direction_near60/polygon_area/Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/bisector_direction_near60/polygon_area/eff99Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/cutting_corridor/eff99Incl_true.txt",
 "../resultats_numeriques/exps_in_thesis/cutting_corridor/Incl_true.txt"]



name_exps = ["med_PaD_eff50_true", "med_PaD_eff70_true", "med_PaD_eff90_true", "med_PaD_eff95_true", "med_PaD_eff99_true",
"med_Area_eff50_true", "med_Area_eff70_true", "med_Area_eff90_true", "med_Area_eff95_true", "med_Area_eff99_true",
"bd_PaD_eff50_true", "bd_PaD_eff70_true", "bd_PaD_eff90_true", "bd_PaD_eff95_true", "bd_PaD_eff99_true",
"bd_Area_eff50_true", "bd_Area_eff70_true", "bd_Area_eff90_true", "bd_Area_eff95_true", "bd_Area_eff99_true",
"near60_med_PaD_eff95_true", "near60_med_PaD_eff99_true", "near60_med_Area_eff95_true", "near60_med_Area_eff99_true", "near60_bd_PaD_eff95_true",
"near60_bd_PaD_eff99_true", "near60_bd_Area_eff95_true", "near60_bd_Area_eff99_true", "med_PaD_eff99_true_cut", "med_PaD_eff95_true_cut"]

name_exps_with_tsm = copy(name_exps)
for i in 1:length(name_exps_with_tsm)
    if !(occursin("near60", name_exps_with_tsm[i]))
        name_exps_with_tsm[i] = string("tsm_",name_exps_with_tsm[i])
    end
end

results_vals,results_times = get_all_vals_from_benchmark(filename_exps)
#result_vals_cuts = copy(results_vals)
#result_times_cuts = copy(results_times)
#results_vals = results_vals[1:28]
#results_times = results_times[1:28]

push!(name_exps,"Rebennack2D")
push!(name_exps,"Rebennack1D")
push!(name_exps,"LinA1D")
push!(name_exps,"KL21")
push!(name_exps,"KL23-LP")
push!(name_exps,"KL23-LPrelaxed")
push!(results_vals,vals_reb2D)
push!(results_vals,vals_reb1D)
push!(results_vals,vals_LinA1D)
push!(results_vals,vals_kazda)
push!(results_vals,vals_KL23LP)
push!(results_vals,vals_KL23relaxed)
push!(results_times,times_reb2D)
push!(results_times,times_reb1D)
push!(results_times,times_LinA1D)
push!(results_times,times_kazda)
push!(results_times,times_KL23LP)
push!(results_times,times_KL23relaxed)

filename_exps_real_dichotomy = ["../resultats_numeriques/exps_in_thesis/real_dichotomy/eff99Incl_true.txt",
"../resultats_numeriques/exps_in_thesis/real_dichotomy/Incl_true.txt"]
results_vals_real_dichotomy,results_times_real_dichotomy = get_all_vals_from_benchmark(filename_exps_real_dichotomy)
append!(results_vals, results_vals_real_dichotomy) # they should be indices 37 and 38
append!(results_times, results_times_real_dichotomy)
append!(name_exps, ["med_PaD_eff99_true_dich", "med_PaD_eff95_true_dich"])

# near60_bd with bd
index_near60 = [25,26,27,28,14,15,19,20]
results_near60_vals = results_vals[index_near60]
results_near60_times = results_times[index_near60]
name_near60_exps = name_exps[index_near60]
#=res_near60 = []
for i in index_near60
    res = [name_near60_exps[i],[]]
    for i in
parameter_competition("near60","tsm",name_near60_exps,results)=#

col_values = get_best_plus_x_bar_plots(results_vals[[27,28,19,20]], name_exps[[27,28,19,20]], "Comparaison entre pp60 et tsm", "plot", true, true)
bests = compute_proportion_best_obj([results_vals[[27,19]],results_vals[[28,20]]], ["pp60","tsm"])
scatter_piece_numbers(results_vals[[27,28,19,20]], name_exps[[27,28,19,20]], "Comparaison entre pp60 et tsm")
geompp60 = compute_group_geometric_mean(results_times[[27,28]])
geomtsm = compute_group_geometric_mean(results_times[[19,20]])
println("geometric mean time pp60 $geompp60\ngeometric mean time tsm $geomtsm")
geompp60vals = compute_group_geometric_mean(results_vals[[27,28]])
geomtsmvals = compute_group_geometric_mean(results_vals[[19,20]])
println("geometric mean vals pp60 $geompp60vals\ngeometric mean vals tsm $geomtsmvals")
#col_values = get_best_plus_x_bar_plots(results_vals[[27,28]], name_exps[[27,28]], "Comparaison entre pp60 et tsm", "plot", true)

# min eff
index_mineff1 = [5,4,3,2,1]
index_mineff2 = [10,9,8,7,6]
index_mineff3 = [15,14,13,12,11]
index_mineff4 = [20,19,18,17,16]
col_values = get_best_plus_x_bar_plots(results_vals[index_mineff1], name_exps[index_mineff1], "Comparaison sur la valeur η avec ProMA et VaT", "plot", true, false)
col_values = get_best_plus_x_bar_plots(results_vals[index_mineff2], name_exps[index_mineff2], "Comparaison sur la valeur η avec ProMA et Aire", "plot", true, false)
col_values = get_best_plus_x_bar_plots(results_vals[index_mineff3], name_exps[index_mineff3], "Comparaison sur la valeur η avec DB et VaT", "plot", true, false)
col_values = get_best_plus_x_bar_plots(results_vals[index_mineff4], name_exps[index_mineff4], "Comparaison sur la valeur η avec DB et Aire", "plot", true, false)
bests = compute_proportion_best_obj([results_vals[index_mineff1],results_vals[index_mineff2],results_vals[index_mineff3],results_vals[index_mineff4]], ["eta99","eta95","eta90","eta70","eta50"])
scatter_piece_numbers(results_vals[index_mineff1], name_exps[index_mineff1], "Comparaison sur la valeur η avec ProMA et VaT")
scatter_piece_numbers(results_vals[index_mineff2], name_exps[index_mineff2], "Comparaison sur la valeur η avec ProMA et Aire")
scatter_piece_numbers(results_vals[index_mineff3], name_exps[index_mineff3], "Comparaison sur la valeur η avec DB et VaT")
scatter_piece_numbers(results_vals[index_mineff4], name_exps[index_mineff4], "Comparaison sur la valeur η avec DB et Aire")
geom50 = compute_group_geometric_mean(results_times[[1,6,11,16]])
geom70 = compute_group_geometric_mean(results_times[[2,7,12,17]])
geom90 = compute_group_geometric_mean(results_times[[3,8,13,18]])
geom95 = compute_group_geometric_mean(results_times[[4,9,14,19]])
geom99 = compute_group_geometric_mean(results_times[[5,10,15,20]])
println("geometric mean min time eff 50 $geom50")
println("geometric mean min time eff 70 $geom70")
println("geometric mean min time eff 90 $geom90")
println("geometric mean min time eff 95 $geom95")
println("geometric mean min time eff 99 $geom99")
geom50vals = compute_group_geometric_mean(results_vals[[1,6,11,16]])
geom70vals = compute_group_geometric_mean(results_vals[[2,7,12,17]])
geom90vals = compute_group_geometric_mean(results_vals[[3,8,13,18]])
geom95vals = compute_group_geometric_mean(results_vals[[4,9,14,19]])
geom99vals = compute_group_geometric_mean(results_vals[[5,10,15,20]])
println("geometric mean min vals eff 50 $geom50vals")
println("geometric mean min vals eff 70 $geom70vals")
println("geometric mean min vals eff 90 $geom90vals")
println("geometric mean min vals eff 95 $geom95vals")
println("geometric mean min vals eff 99 $geom99vals")


# DB vs ProMA (bd vs med)
index_direction = [14,15,19,20,4,5,9,10]
col_values = get_best_plus_x_bar_plots(results_vals[index_direction], name_exps[index_direction], "Comparaison entre DB et ProMA", "plot", true, true)
bests = compute_proportion_best_obj([results_vals[[14,4]],results_vals[[15,5]],results_vals[[19,9]],results_vals[[20,10]]], ["DB","ProMA"])
scatter_piece_numbers(results_vals[index_direction], name_exps[index_direction], "Comparaison entre DB et ProMA")
geomDB = compute_group_geometric_mean(results_times[[14,15,19,20]])
geomProMA = compute_group_geometric_mean(results_times[[4,5,9,10]])
println("geometric mean time DB $geomDB")
println("geometric mean time ProMA $geomProMA")
geomDBvals = compute_group_geometric_mean(results_vals[[14,15,19,20]])
geomProMAvals = compute_group_geometric_mean(results_vals[[4,5,9,10]])
println("geometric mean vals DB $geomDBvals")
println("geometric mean vals ProMA $geomProMAvals")

# Aire vs VaT (Area vs PaD)
index_score = [9,10,19,20,4,5,14,15]
col_values = get_best_plus_x_bar_plots(results_vals[index_score], name_exps[index_score], "Comparaison entre Aire et VaT", "plot", true, true)
bests = compute_proportion_best_obj([results_vals[[9,4]],results_vals[[10,5]],results_vals[[19,14]],results_vals[[20,15]]], ["Aire","VaT"])
scatter_piece_numbers(results_vals[index_score], name_exps[index_score], "Comparaison entre Aire et VaT")
geomAire = compute_group_geometric_mean(results_times[[9,10,19,20]])
geomVaT = compute_group_geometric_mean(results_times[[4,5,14,15]])
println("geometric mean time Aire $geomAire")
println("geometric mean time VaT $geomVaT")
geomAirevals = compute_group_geometric_mean(results_vals[[9,10,19,20]])
geomVaTvals = compute_group_geometric_mean(results_vals[[4,5,14,15]])
println("geometric mean vals Aire $geomAirevals")
println("geometric mean vals VaT $geomVaTvals")

# meilleure heuristique (pour 95 et pour 99)
index_score1 = [4,9,14,19]
index_score2 = [5,10,15,20]
col_values = get_best_plus_x_bar_plots(results_vals[index_score1], name_exps[index_score1], "Comparaison des heuristiques pour eta=95", "plot", true, false)
col_values = get_best_plus_x_bar_plots(results_times[index_score1], name_exps[index_score1], "Temps des heuristiques pour eta=95", "plot", true, false)
bests = compute_proportion_best_obj([results_vals[index_score1]], ["Proma_VaT","ProMA_Aire","DB_VaT","DB_Aire"])
scatter_piece_numbers(results_vals[index_score1], name_exps[index_score1], "Comparaison des heuristiques pour eta=95")
geomProMA_VaT_95 = compute_group_geometric_mean(results_times[[4]])
geomProMA_Aire_95 = compute_group_geometric_mean(results_times[[9]])
geomDB_VaT_95 = compute_group_geometric_mean(results_times[[14]])
geomDB_Aire_95 = compute_group_geometric_mean(results_times[[19]])
println("geometric mean time proma vat 95 $geomProMA_VaT_95")
println("geometric mean time proma aire 95 $geomProMA_Aire_95")
println("geometric mean time DB vat 95 $geomDB_VaT_95")
println("geometric mean time DB aire 95 $geomDB_Aire_95")
geomProMA_VaT_95vals = compute_group_geometric_mean(results_vals[[4]])
geomProMA_Aire_95vals = compute_group_geometric_mean(results_vals[[9]])
geomDB_VaT_95vals = compute_group_geometric_mean(results_vals[[14]])
geomDB_Aire_95vals = compute_group_geometric_mean(results_vals[[19]])
println("geometric mean vals proma vat 95 $geomProMA_VaT_95vals")
println("geometric mean vals proma aire 95 $geomProMA_Aire_95vals")
println("geometric mean vals DB vat 95 $geomDB_VaT_95vals")
println("geometric mean vals DB aire 95 $geomDB_Aire_95vals")
col_values = get_best_plus_x_bar_plots(results_vals[index_score2], name_exps[index_score2], "Comparaison des heuristiques pour eta=99", "plot", true, false)
col_values = get_best_plus_x_bar_plots(results_times[index_score2], name_exps[index_score2], "Temps des heuristiques pour eta=99", "plot", true, false)
bests = compute_proportion_best_obj([results_vals[index_score2]], ["Proma_VaT","ProMA_Aire","DB_VaT","DB_Aire"])
scatter_piece_numbers(results_vals[index_score2], name_exps[index_score2], "Comparaison des heuristiques pour eta=99")
geomProMA_VaT_99 = compute_group_geometric_mean(results_times[[5]])
geomProMA_Aire_99 = compute_group_geometric_mean(results_times[[10]])
geomDB_VaT_99 = compute_group_geometric_mean(results_times[[15]])
geomDB_Aire_99 = compute_group_geometric_mean(results_times[[20]])
println("geometric mean time proma vat 99 $geomProMA_VaT_99")
println("geometric mean time proma aire 99 $geomProMA_Aire_99")
println("geometric mean time DB vat 99 $geomDB_VaT_99")
println("geometric mean time DB aire 99 $geomDB_Aire_99")
geomProMA_VaT_99vals = compute_group_geometric_mean(results_vals[[5]])
geomProMA_Aire_99vals = compute_group_geometric_mean(results_vals[[10]])
geomDB_VaT_99vals = compute_group_geometric_mean(results_vals[[15]])
geomDB_Aire_99vals = compute_group_geometric_mean(results_vals[[20]])
println("geometric mean vals proma vat 99 $geomProMA_VaT_99vals")
println("geometric mean vals proma aire 99 $geomProMA_Aire_99vals")
println("geometric mean vals DB vat 99 $geomDB_VaT_99vals")
println("geometric mean vals DB aire 99 $geomDB_Aire_99vals")


# print table with comparison to state of the art without times because it is incomparable
index_sota = [4,5,31,33,34,35,36]
println("table with values")
s_vals = build_table(expr_functions, deltas, results_vals[index_sota], results_times[index_sota],"val")
println("table with times")
s_times = build_table(expr_functions, deltas, results_times[index_sota], [[]], "val", "nothing", true)
name_sota = ["DN95","DN99","RK2D","LinA1D","KL21","KL23-LP","KL23-LPrelaxed"]

col_values = get_best_plus_x_bar_plots(results_vals[index_sota], name_exps[index_sota], "Comparaison sota", "plot", true, false)
#col_values = get_best_plus_x_bar_plots(results_times[index_sota], name_exps[index_sota], "Comparaison sota times", "plot", true, false, true)
bests = compute_proportion_best_obj([results_vals[index_sota]], name_sota)
#scatter_piece_numbers(results_vals[index_score], name_exps[index_score], "Comparaison entre Aire et VaT")
for i in 1:length(index_sota)
    geom = compute_group_geometric_mean(results_vals[index_sota[i]])
    println("geometric mean vals $(name_sota[i]) $geom")
end

# display DN95 and DN99 times
println("temps DN95")
for i in 1:9
    for j in 1:5
        print("$(results_times[4][(i-1)*5+j]) ")
    end
    println()
end
println("temps DN99")
for i in 1:9
    for j in 1:5
        print("$(results_times[5][(i-1)*5+j]) ")
    end
    println()
end


# print table with comparison to sota without DN95 because it decreases the results of DN99
index_sota = [5,31,33,34,35,36]
println("table with values without DN95")
s_vals = build_table(expr_functions, deltas, results_vals[index_sota], results_times[index_sota],"val")
println("table with times without DN95")
s_times = build_table(expr_functions, deltas, results_times[index_sota], [[]], "val", "nothing", true)


# print table with cutting versions and comparison to sota without DN95 because it decreases the results of DN99
index_sota = [37,31,33,34,35,36]
println("table cutting version with values without DN95")
s_vals = build_table(expr_functions, deltas, results_vals[index_sota], results_times[index_sota],"val")
println("table cutting version with times without DN95")
s_times = build_table(expr_functions, deltas, results_times[index_sota], [[]], "val", "nothing", true)


# print table with comparison between two best heuristics without cutting pieces in the PWL inner corridor (results in the manuscript thesis), and two with cutting pieces (option cutting_corridor)
index_sota = [4,5,29,30]
println("table with cutting tests")
s_vals = build_table(expr_functions, deltas, results_vals[index_sota], results_times[index_sota],"val")
println("table with times of cutting tests")
s_times = build_table(expr_functions, deltas, results_times[index_sota], [[]], "val", "nothing", true)
col_values = get_best_plus_x_bar_plots(results_vals[index_sota], name_exps[index_sota], "Comparaison sota", "plot", true, false)

# print table with comparison between two best heuristics without cutting pieces in the PWL inner corridor, two with cutting pieces (option cutting_corridor), and two with real dichotomy (option real_dichotomy)
index_sota = [4,30,38,5,29,37]
println("table with real dichotomy tests")
s_vals = build_table(expr_functions, deltas, results_vals[index_sota], results_times[index_sota],"val")
println("table with times of real dichotomy tests")
s_times = build_table(expr_functions, deltas, results_times[index_sota], [[]], "val", "nothing", true)
col_values = get_best_plus_x_bar_plots(results_vals[index_sota], name_exps[index_sota], "Comparison improvements of heuristic", "plot", true, false)

#=
col_values = get_best_plus_x_bar_plots(results_vals[[27,19]], name_exps[[27,19]], "Comparaison entre pp60 et tsm mineff95", "plot", true, true)
col_values = get_best_plus_x_bar_plots(results_vals[[28,20]], name_exps[[28,20]], "Comparaison entre pp60 et tsm mineff99", "plot", true, true)
geompp601 = compute_group_geometric_mean(results_times[27])
geomtsm1 = compute_group_geometric_mean(results_times[[19]])
geompp602 = compute_group_geometric_mean(results_times[28])
geomtsm2 = compute_group_geometric_mean(results_times[[20]])
println("geometric mean pp601 $geompp601\ngeometric mean tsm1 $geomtsm1")
println("geometric mean pp602 $geompp602\ngeometric mean tsm2 $geomtsm2")
#col_values = get_best_plus_x_bar_plots(results_vals[[27,28]], name_exps[[27,28]], "Comparaison entre pp60 et tsm", "plot", true)

# min eff
index_mineff1 = [5,4,3,2,1]
index_mineff2 = [10,9,8,7,6]
index_mineff3 = [15,14,13,12,11]
index_mineff4 = [20,19,18,17,16]
col_values = get_best_plus_x_bar_plots(results_vals[index_mineff1], name_exps[index_mineff1], "Comparaison sur la valeur η avec ProMA et VaT", "plot", true, false)
col_values = get_best_plus_x_bar_plots(results_vals[index_mineff2], name_exps[index_mineff2], "Comparaison sur la valeur η avec ProMA et Aire", "plot", true, false)
col_values = get_best_plus_x_bar_plots(results_vals[index_mineff3], name_exps[index_mineff3], "Comparaison sur la valeur η avec DB et VaT", "plot", true, false)
col_values = get_best_plus_x_bar_plots(results_vals[index_mineff4], name_exps[index_mineff4], "Comparaison sur la valeur η avec DB et Aire", "plot", true, false)
geom50 = compute_group_geometric_mean(results_times[[1,6,11,16]])
geom70 = compute_group_geometric_mean(results_times[[2,7,12,17]])
geom90 = compute_group_geometric_mean(results_times[[3,8,13,18]])
geom95 = compute_group_geometric_mean(results_times[[4,9,14,19]])
geom99 = compute_group_geometric_mean(results_times[[5,10,15,20]])
println("geometric mean min eff 50 $geom50")
println("geometric mean min eff 70 $geom70")
println("geometric mean min eff 90 $geom90")
println("geometric mean min eff 95 $geom95")
println("geometric mean min eff 99 $geom99")


# DB vs ProMA (bd vs med)
index_direction = [14,15,19,20,4,5,9,10]
col_values = get_best_plus_x_bar_plots(results_vals[index_direction], name_exps[index_direction], "Comparaison entre DB et ProMA", "plot", true, true)
geomDB = compute_group_geometric_mean(results_times[[14,15,19,20]])
geomProMA = compute_group_geometric_mean(results_times[[4,5,9,10]])
println("geometric mean DB $geomDB")
println("geometric mean ProMA $geomProMA")

# Aire vs VaT (Area vs PaD)
index_score = [9,10,19,20,4,5,14,15]
col_values = get_best_plus_x_bar_plots(results_vals[index_score], name_exps[index_score], "Comparaison entre Aire et VaT", "plot", true, true)
geomAire = compute_group_geometric_mean(results_times[[9,10,19,20]])
geomVaT = compute_group_geometric_mean(results_times[[4,5,14,15]])
println("geometric mean Aire $geomAire")
println("geometric mean VaT $geomVaT")
=#

# SOTA
#index_sota = find values within my heuristics. One with 95 and one with 99

#=
#=filenames_surf = ["polygon_surface_new/lip2Incl_false.txt","polygon_surface_new/lip2Incl_true.txt","polygon_surface_new/lip2hybrid.txt","polygon_surface_new/Incl_false.txt","polygon_surface_new/Incl_true.txt","polygon_surface_new/hybrid.txt","polygon_surface_new/eff99Incl_false.txt","polygon_surface_new/eff99Incl_true.txt","polygon_surface_new/eff99hybrid.txt"]
filenames_PDTV = ["PDTV/lip2Incl_false.txt","PDTV/lip2Incl_true.txt","PDTV/lip2hybrid.txt","PDTV/Incl_false.txt","PDTV/Incl_true.txt","PDTV/hybrid.txt","PDTV/eff99Incl_false.txt","PDTV/eff99Incl_true.txt","PDTV/eff99hybrid.txt"]
filenames_SPDTV = ["SPDTV/lip2Incl_false.txt","SPDTV/lip2Incl_true.txt","SPDTV/lip2hybrid.txt","SPDTV/Incl_false.txt","SPDTV/Incl_true.txt","SPDTV/hybrid.txt","SPDTV/eff99Incl_false.txt","SPDTV/eff99Incl_true.txt","SPDTV/eff99hybrid.txt"]
extra_filenames = ["near60/test_near60.txt"]
filenames = filenames_surf
append!(filenames,filenames_PDTV)
append!(filenames,filenames_SPDTV)
#filenames = filenames_surf[4:end]
#append!(filenames,filenames_PDTV[4:end])
#append!(filenames,filenames_SPDTV[4:end])
append!(filenames,extra_filenames)

name_exps_surf = ["SurfLip2","SurfLip2Incl","SurfLip2Hybrid","SurfEff95","SurfEff95Incl","SurfEff95Hybrid","SurfEff99","SurfEff99Incl","SurfEff99Hybrid"]
name_exps_PDTV = ["PDTVLip2","PDTVLip2Incl","PDTVLip2Hybrid","PDTVEff95","PDTVEff95Incl","PDTVEff95Hybrid","PDTVEff99","PDTVEff99Incl","PDTVEff99Hybrid"]
name_exps_SPDTV = ["SPDTVLip2","SPDTVLip2Incl","SPDTVLip2Hybrid","SPDTVEff95","SPDTVEff95Incl","SPDTVEff95Hybrid","SPDTVEff99","SPDTVEff99Incl","SPDTVEff99Hybrid"]
extra_names = ["near60PDTV95Incl"]
name_exps = name_exps_surf
append!(name_exps,name_exps_PDTV)
append!(name_exps,name_exps_SPDTV)
#name_exps = name_exps_surf[4:end]
#append!(name_exps,name_exps_PDTV[4:end])
#append!(name_exps,name_exps_SPDTV[4:end])
append!(name_exps,extra_names)=#

filenames = ["final_exps_old/near60PDTVEff95.txt","final_exps_old/PDTVEff95.txt","final_exps_old/PDTVEff99.txt","final_exps_old/SurfEff95.txt","final_exps_old/SurfEff99.txt","final_exps_old/PDTVEff95med.txt","final_exps_old/PDTVEff99med.txt"]
name_exps = ["near60Eff95","PDTVEff95","PDTVEff99","SurfEff95","SurfEff99","PDTVEff95med","PDTVEff99med"]
bonus_filenames = ["PDTV/eff99Incl_true.txt","PDTV/eff99hybrid.txt","SDI/eff99hybrid.txt"]
bonus_name_exps = ["PDTVE99Incl","PDTVEff99Hybrid","SDIEff99Hybrid"]

results_vals,results_times = get_all_vals_from_benchmark(filenames)


push!(name_exps,"Rebennack 2D")
push!(name_exps,"Rebennack 1D")
push!(name_exps,"LinA 1D")
push!(name_exps,"Kazda")
push!(results_vals,vals_reb2D)
push!(results_vals,vals_reb1D)
push!(results_vals,vals_LinA1D)
push!(results_vals,vals_kazda)
push!(results_times,times_reb2D)
push!(results_times,times_reb1D)
push!(results_times,times_LinA1D)
push!(results_times,times_kazda)

functions = [4,5,6,7,8,9]
deltas = [1,2]
res_vals = select_instances(results_vals,functions,deltas)
col_values = get_best_plus_x_bar_plots(res_vals[[2,3,4,5]], name_exps[[2,3,4,5]], "Best heuristics and state of the art")
deltas = [4,5]
res_vals = select_instances(results_vals,functions,deltas)
col_values = get_best_plus_x_bar_plots(res_vals[[2,3,4,5]], name_exps[[2,3,4,5]], "Best heuristics and state of the art")
deltas = [1,2]
res_vals = select_instances(results_vals,functions,deltas)
col_values = get_best_plus_x_bar_plots(res_vals[[2,4]], name_exps[[2,4]], "Best heuristics and state of the art")
deltas = [4,5]
res_vals = select_instances(results_vals,functions,deltas)
col_values = get_best_plus_x_bar_plots(res_vals[[2,4]], name_exps[[2,4]], "Best heuristics and state of the art")
deltas = [1,2]
res_vals = select_instances(results_vals,functions,deltas)
col_values = get_best_plus_x_bar_plots(res_vals[[3,5]], name_exps[[3,5]], "Best heuristics and state of the art")
deltas = [4,5]
res_vals = select_instances(results_vals,functions,deltas)
col_values = get_best_plus_x_bar_plots(res_vals[[3,5]], name_exps[[3,5]], "Best heuristics and state of the art")
deltas = [1]
res_vals = select_instances(results_vals,functions,deltas)
col_values = get_best_plus_x_bar_plots(res_vals[[2,4]], name_exps[[2,4]], "Best heuristics and state of the art")
deltas = [5]
res_vals = select_instances(results_vals,functions,deltas)
col_values = get_best_plus_x_bar_plots(res_vals[[2,4]], name_exps[[2,4]], "Best heuristics and state of the art")
deltas = [1]
res_vals = select_instances(results_vals,functions,deltas)
col_values = get_best_plus_x_bar_plots(res_vals[[3,5]], name_exps[[3,5]], "Best heuristics and state of the art")
deltas = [5]
res_vals = select_instances(results_vals,functions,deltas)
col_values = get_best_plus_x_bar_plots(res_vals[[3,5]], name_exps[[3,5]], "Best heuristics and state of the art")


functions = [1,2,3,4,5,6,7,8,9]
deltas = [1,2,3]
res_vals = select_instances(results_vals,functions,deltas)
col_values = get_best_plus_x_bar_plots(res_vals[[6,7,8,9,11]], name_exps[[6,7,8,9,11]], "Best heuristics and state of the art")

deltas = [4,5]
res_vals = select_instances(results_vals,functions,deltas)
col_values = get_best_plus_x_bar_plots(res_vals[[6,7,8,9,11]], name_exps[[6,7,8,9,11]], "Best heuristics and state of the art")

#results = select_instances(results,functions,deltas)
#title = create_title_comparison(functions, deltas)
#col_values = get_best_plus_x_bar_plots(results, name_exps, title)
#col_times = get_best_plus_x_bar_plots(results_times[1:(end-2)], name_exps[1:(end-2)])
col_values = get_best_plus_x_bar_plots(results_vals[[2,3,4,5,6,7]], name_exps[[2,3,4,5,6,7]], "Best heuristics")
col_values = get_best_plus_x_bar_plots(results_vals[[2,3,6,7]], name_exps[[2,3,6,7]], "Best PDTV heuristics")
col_values = get_best_plus_x_bar_plots(results_vals[[3,7]], name_exps[[3,7]], "Best PDTVEff99 heuristics")
col_values = get_best_plus_x_bar_plots(results_vals[[2,6]], name_exps[[2,6]], "Best PDTVEff95 heuristics")
col_values = get_best_plus_x_bar_plots(results_vals[[9,10]], name_exps[[9,10]], "Best PDTVEff95 heuristics")
col_values = get_best_plus_x_bar_plots(results_vals[[6,7,8,9,11]], name_exps[[6,7,8,9,11]], "Best heuristics and state of the art")
#col_values = get_best_plus_x_bar_plots(results_vals[[7,8]], name_exps[[7,8]], "Best heuristics and state of the art")
#col_values = get_best_plus_x_bar_plots(results_vals[[3,6,7,8]], name_exps[[3,6,7,8]], "Eff99 and co")
#col_values = get_best_plus_x_bar_plots(results_vals[[1,2,3,10,11,12,19,20,21]], name_exps[[1,2,3,10,11,12,19,20,21]], "Lip2 heuristics")
#col_values = get_best_plus_x_bar_plots(results_vals[[4,5,6,13,14,15,22,23,24]], name_exps[[4,5,6,13,14,15,22,23,24]], "Eff95 heuristics")
#col_values = get_best_plus_x_bar_plots(results_vals[[7,8,9,16,17,18,25,26,27]], name_exps[[7,8,9,16,17,18,25,26,27]], "Eff99 heuristics")
#col_values = get_best_plus_x_bar_plots(results_vals[[12,15,18,29,30]], ["Regular refinement", "95% efficiency refinement", "99% efficiency refinement", "Rebennack 1D heuristic", "Rebennack 2D heuristic"], "best heuristics and state of the art")
#col_values = get_best_plus_x_bar_plots(results_vals[[15,28]], ["95% efficiency refinement", "Nearest 60 degree angle first vertex"], "Interest of choosing all vertices as first vertex")
#results = select_instances(results_vals,[1,2,3,4,5,6,7,8],[1,2,3,4,5])
#col_values = get_best_plus_x_bar_plots(results[[13,14,15]], name_exps[[13,14,15]], "inclination option comparison for 95% efficiency refinement")
col_values = get_best_plus_x_bar_plots(results_vals[[1,2,3,4,5,6,7]], name_exps[[1,2,3,4,5,6,7]], "Some of my heuristics")
#col_values = get_best_plus_x_bar_plots(results_vals[28:30], name_exps[28:30])

#=p = plot()
for i in 1:length(results_times)
    p = plot!(results_times[i], label="$(name_exps[i])", yaxis=:log)
    println("expérience $(name_exps[i])\t de temps total\t $(sum(results_times[i]))")
end
display(p)

p = plot(results_vals[28],label="$(name_exps[28])")
p = plot!(results_vals[29],label="$(name_exps[29])")
#p = plot!(results_vals[30],label="$(name_exps[30])")
display(p)
#p = plot(results_times[28],label="$(name_exps[28])")
#p = plot!(results_times[30],label="$(name_exps[30])")
=#

expr_functions = []
for i in 1:9
    push!(expr_functions, get_limits_and_func(i)[6])
end

deltas1 = [1.5,1.0,0.5,0.25,0.1]
deltas2 = deltas1
deltas3 = [1.0,0.5,0.25,0.1,0.05]
deltas4 = [0.1,0.05,0.03,0.01,0.001]
deltas5 = deltas3
deltas6 = [0.5,0.25,0.1,0.05,0.03]
deltas7 = deltas3
deltas8 = deltas3
deltas9 = deltas3
deltas = [deltas1,deltas2,deltas3,deltas4,deltas5,deltas6,deltas7,deltas8,deltas9]

s = build_table(expr_functions, deltas, results_vals[[6,7,8,9,10,11]], results_times[[6,7,8,9,10,11]])
s = build_table(expr_functions, deltas, results_vals[[6,7,8,10,11]], results_times[[6,7,8,10,11]])
s = build_table(expr_functions, deltas, results_vals[[1,2,3,4,5,6,7]], results_times[[1,2,3,4,5,6,7]])

sum_95 = sum(results_vals[6][1:end])
sum_99 = sum(results_vals[7][1:end])
time_95 = sum(results_times[6][1:end])
time_99 = sum(results_times[7][1:end])
println("somme 95 $sum_95 temps $time_95")
println("somme 99 $sum_99 temps $time_99")
println("ratio pièces total $(sum_99/sum_95) temps total $(time_99/time_95)")
ratio_pieces = sum([results_vals[7][i]/results_vals[6][i] for i in 1:45])/45
ratio_times = sum([results_times[7][i]/results_times[6][i] for i in 1:45])/45
println("ratio pièces moyen : $ratio_pieces temps moyen : $ratio_times")

=#
