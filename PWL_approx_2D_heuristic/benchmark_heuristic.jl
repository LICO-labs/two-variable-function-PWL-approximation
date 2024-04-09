using Plots

include("main.jl")

function benchmark_flex_heuristic3(INCLINED_BOUNDING = true, EVAL_F_STR = "PDTV", BOUNDING_METHOD = "eff", MIN_EFF = 0.95, firstv = "all", filename = "../resultats_numeriques/resultats.txt", FDH = "mean_along_edge", compute_piece_version = "real_dichotomy"; delta_values = 1:5, function_number_values = 1:9)
    # fait tourner les 45 instances de Rebennack avec flex_heuristic_procedure_for_benchmark
    # en changeant en priorité les fonctions et ensuite les deltas (instances courtes résolues d'abord)

    # ajout de deux saut de ligne dans le fichier résultat
    file = open("../resultats_numeriques/resultats.txt", "a")
    println(file,"\n")
    close(file)

    # compilation des fonctions
    println("blank launch before the benchmark to load functions")
    launch_flex_heuristic(1,Absolute(1.5))
    println("\n\n\nfunctions compiled. Start of real experiments.\n\n\n")

    # définition des deltas
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

    # répertorie les instances résolues correctement
    well_solved = zeros(5,9)

    # répertorie les temps
    times = zeros(5,9)

    # répertorie le nombre de polygones
    number_of_polygons = zeros(5,9)


    for i in delta_values
        for j in function_number_values
            num_func = j
            A,B,C,D,f,str_exprf = get_limits_and_func(num_func)

            delta = deltas[num_func][i]
            err = Absolute(delta)

            domain = A,B,C,D
            if BOUNDING_METHOD == "lipschitz2"
                n_corridor = [20,20]
            elseif BOUNDING_METHOD == "eff"
                n_corridor = [1,1]
            end

            DSR = 0.0625
            inclined_bounding = INCLINED_BOUNDING # false ou true ou "hybrid"
            bounding_method = BOUNDING_METHOD
            eval_f_str = EVAL_F_STR # "polygon_surface" or"remaining_error" or "SDI" or "SSDI"
            min_eff = MIN_EFF
            n_eval_sample = 200

            arguments = Dict("num_func"=>num_func,"err"=>err,"str_exprf"=>str_exprf,"f"=>f,"domain"=>domain,"LP_SOLVER"=>"Gurobi","n_corridor"=>n_corridor,"DSR"=>DSR,"time_limit"=>3600.0,"firstv"=>firstv)
            arguments["bonus_plot_infos"] = []
            arguments["BM"] = bounding_method
            arguments["inclined_bounding"] = inclined_bounding
            arguments["LP_SOLVER"] = "Gurobi"
            arguments["FDH"] = FDH # anciennement new_piece_heuristic
            arguments["eval_f_str"] = eval_f_str
            if eval_f_str == "SDI" || eval_f_str == "SSDI" || eval_f_str == "APCI" || eval_f_str == "PDTV" || eval_f_str == "PaD" || eval_f_str == "SPDTV"
                arguments["n_eval_sample"] = n_eval_sample
            end

            # useful in function efficiency_refinement_bounding
            if bounding_method == "eff"
                arguments["min_eff"] = min_eff
            end

            if compute_piece_version == "real_dichotomy"
                compute_piece = compute_piece_by_forward_heuristic_real_dichotomy
            elseif compute_piece_version == "cutting_corridor"
                compute_piece = compute_piece_by_forward_heuristic_cutting_corridor
            else
                error("unrecognised compute_piece_version argument; real_dichotomy or cutting_corridor are the only options")
            end
            arguments["compute_piece"] = compute_piece
            arguments["v3"] = "opt"

            name_func_new_piece = "forward_heuristic+test_all_vertex_for_new_piece_v5"
            display_and_save = true
            disp = false
            plot_partial_pwl = false


            bool,t,polygons = flex_heuristic_procedure_for_benchmark(f, arguments, filename, name_func_new_piece, display_and_save = display_and_save, plot_partial_pwl = plot_partial_pwl)

            if bool
                well_solved[i,j] = 1
                times[i,j] = t
                number_of_polygons[i,j] = length(polygons)
            else
                well_solved[i,j] = 0
            end
        end
    end
    return well_solved,times,number_of_polygons
end

# experiments for real_dichotomy version:
#ws,ts,ps = benchmark_flex_heuristic3(true, "PDTV", "eff", 0.95, "all", "../resultats_numeriques/exps_in_thesis/real_dichotomy/Incl_true.txt","mean_along_edge");
#ws,ts,ps = benchmark_flex_heuristic3(true, "PDTV", "eff", 0.99, "all", "../resultats_numeriques/exps_in_thesis/real_dichotomy/eff99Incl_true.txt","mean_along_edge");
