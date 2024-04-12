This folder contains files with the results of different experiments for the algorithm coded in this git:
- results mentioned in [1] (for the edited version: https://link.springer.com/chapter/10.1007/978-3-031-18530-4_9; for the unedited version: https://hal.laas.fr/hal-03629850v2/document) are in folder "ISCO2022_AD_SUN_solutions",
- results mentioned in my manuscript thesis [2] (in french, dating from november 2023, see https://theses.hal.science/tel-04474172) are in folder "thesis results",
- results more recent than the manuscript thesis are in folder "improvement from thesis results".

The .txt files contain the result of one instance per line, with the string "::" as separator between informations.

Most of the result files contain the full results, you can convert them to simplified result files with function "simplify_ISCO2022_results" as explained below. Result files with "_simplified" in the name are simplified result files. They contain the essential informations for comparison purposes:
        - the name of the method
        - the function name (names: L1,L2,N1,...,N7)
        - the function expression (X*X-Y*Y...)
        - the error allowed (Absolute(delta) for absolute error of delta)
        - the domain (in format [A,B,C,D] meaning [A,B]x[C,D])
        - the number of pieces of the solution
        - the execution time in seconds
        - the solution is a Vector{Vector{Vector{Float64}}} in julia language (described below)

solution is a Vector{Vector{Vector{Float64}}}: a vertex of a piece is a vector of $\mathbb{R}^3, the first two coordinates give the position in the domain, the third one the value of the PWl function in that position. A piece is a set of vertices, and a PWl function is a set of pieces. Thus, a solution is represented as a vector of pieces, which are vectors of vertices, which are vectors of Float64.

The julia language file "parse_and_plot.jl" contains code in julia with a function called "parse_and_plot_instance", that plots and saves the result of an instance.

    Syntax of the "parse_and_plot_instance" function call (see line 257 in 'parse_and_plot.jl' file):
        T_set, domain, delta, n_pieces, cpu_time, str_exprf = parse_and_plot_instance(method_name, function_name, precision)

    Accepted input parameters values:
        - method_name (string): Any string you want to describe the method.
        - function_name (string): "L1, "L2", "N1", "N2", "N3", "N4", "N5", "N6", "N7"
        - precision (float64): should be equal to the requested delta value (please refer to the table of results in appendix 'https://homepages.laas.fr/sungueve/Docs/IJ/Appendix_ISCO_2022.pdf')

    Examples:
        T_set, domain, delta, n_pieces, cpu_time, str_exprf = parse_and_plot_instance("DN99", "N1", 0.05)
        T_set, domain, delta, n_pieces, cpu_time, str_exprf = parse_and_plot_instance("DN99", "N1", 0.1)

    (Figures of the solutions domains are displayed but also saved in the folder 'images/' in a 'png' format)


[1]: Duguet, A., Ngueveu, S.U. (2022). Piecewise Linearization of Bivariate Nonlinear Functions: Minimizing the Number of Pieces Under a Bounded Approximation Error. In: Ljubić, I., Barahona, F., Dey, S.S., Mahjoub, A.R. (eds) Combinatorial Optimization. ISCO 2022. Lecture Notes in Computer Science, vol 13526. Springer, Cham. https://doi.org/10.1007/978-3-031-18530-4_9
[2]: Aloïs Duguet. Approximation linéaire par morceaux de fonctions de deux variables avec erreur bornée pour la résolution de problèmes d'optimisation non linéaire mixte en nombres entiers. Autre [cs.OH]. Institut National Polytechnique de Toulouse - INPT, 2023. Français.
