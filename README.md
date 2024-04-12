# two-variable-function-PWL-approximation

## License
License CC BY-NC-SA
Developed by Aloïs Duguet

## Goal
This document describes how to use the function PWL2D_heuristic to build a piecewise linear function approximating a function within a given approximating error.
It is the implementation of the algorithm described in [1] and [2], with some further changes decreasing in average the number of pieces but increasing the computation time.

[1]: Duguet, A., Ngueveu, S.U. (2022). Piecewise Linearization of Bivariate Nonlinear Functions: Minimizing the Number of Pieces Under a Bounded Approximation Error. In: Ljubić, I., Barahona, F., Dey, S.S., Mahjoub, A.R. (eds) Combinatorial Optimization. ISCO 2022. Lecture Notes in Computer Science, vol 13526. Springer, Cham. https://doi.org/10.1007/978-3-031-18530-4_9
[2]: Aloïs Duguet. Approximation linéaire par morceaux de fonctions de deux variables avec erreur bornée pour la résolution de problèmes d'optimisation non linéaire mixte en nombres entiers. Autre [cs.OH]. Institut National Polytechnique de Toulouse - INPT, 2023. Français. (in french, access [here](https://theses.hal.science/tel-04474172))

## Software
Software needed:
+ julia version 1.6.2 (more recent version of julia possible with minor changes in libraries)
+ Gurobi as LP solver (others may be used with a slight change in the code if a julia library link them to the JuMP library; line "model = Model(Gurobi.optimizer)" of functions feasibility_model2 and feasibility_model_1D need to be adapted to your LP solver), a free student version is available on the optimizer website

Julia library needed: LinA, JuMP v0.21.10, Gurobi v0.9.14, LinearAlgebra, Plots v1.19.3, Clipper v0.6.1, ForwardDiff v0.10.21, IntervalArithmetic v0.18.2, AutoGrad v1.2.4, PolygonOps v0.1.1, TimerOutputs v0.5.12, StringEncodings v0.3.5

Install those libraries with:

```
]
add https://github.com/LICO-labs/LinA.jl.git
add JuMP Gurobi LinearAlgebra Plots Clipper ForwardDiff IntervalArithmetic@0.18.2 AutoGrad PolygonOps GLPK Cbc AmplNLWriter Couenne_jll AngleBetweenVectors Polyhedra TimerOutputs StringEncodings
```


## Description of the main function

```
function PWL2D_heuristic(f, str_exprf, err, domain; LP_SOLVER = "GLPK", not_rectangular_domain = false, num_func = "UNK", n_corridor = [1,1], DSR = 0.0625, LP_time_limit = 3600.0, firstv = "all", bounding_method = "eff", min_eff = 0.95, inclined_bounding = true, FDH = "mean_along_edge", eval_f_str = "PDTV", n_eval_sample = 200, save_solution = false, plot_pwl = false, plot_partial_pwl = false)
```

Easy to use function that builds a piecewise linear function approximating function f of two variables with approximation error err on domain domain. See [1] for more details.

### Inputs:
+ **f** (generic function): function to piecewise linearize. Must take a single Vector{Float64} as argument.
  
Example:    f(x) = x[1]^2-x[2]^2
+ **str_exprf** (String): describes the formula of function f with "X" as first variable and "Y" as second variable.

Example:    "X*X-Y*Y".
+ **err** (LinA.ErrorType): describes the type of error wanted (absolute or relative) and its characteristics.
  + Absolute(delta), delta>0 means piecewise linear function g constructed respect g(x) in [f(x)-delta,f(x)+delta]
  + Relative(epsilon), 0<epsilon<100 means piecewise linear function g constructed respect g(x) in [f(x)*(1-epsilon/100),f(x)*(1+epsilon/100)]

WARNING: the Relative error type is only supported if function f is strictly positive.

Examples:   Absolute(1.5), Relative(10.0) (for 10%).
+ **domain** (4-element Vector{Float64}): describes the variables domains. [A,B,C,D] means that X ∈ [A,B] and Y ∈ [C,D].
        If the domain is not a cartesian product, put as domain the smallest cartesian product containing the domain wanted and use the not_rectangular_domain option to describe the real domain.

Example:    [0,7.5,1.5,3] means X ∈ [0,7.5] and Y ∈ [1.5,3].
+ **LP_SOLVER** (String): describes the name of the LP solver to use.
        Options: "GLPK", "Gurobi" and "Cbc".
        Gurobi and GLPK will not output anything while Cbc will output information on each LP.
        Gurobi is likely the fastest, followed by Cbc and far away by GLPK.

Example:    "GLPK"
+ **not_rectangular_domain** (Bool or Vector{Vector{Float64}}): Complement the option domain. If the domain of the linearization is a cartesian product, put false.
        If the domain of the linearization is not a cartesian product, put the list of vertices of the polygonal convex domain in counter clockwise fashion.

WARNING: the algorithm needs counter clockwise fashion to work as intended.

Examples:   triangle with vertices (0,0), (0,1) and (2,2) => [[0,0],[2,2],[0,1]] or [[2,2],[0,1],[0,0]] or [[0,1],[0,0],[2,2]]; the first vertex may change the solution.
+ **num_func** (String or Int64, or anything that can be printed with function print): function name printed if option save_solution == true.

Examples:    "inverse_function" or 10.
+ **n_corridor** (2-element Vector{Int64}): starting grid size when computing a piecewise linear inner corridor (PWL inner corridor).

Example:    [20,10].
+ **DSR** (Float64): number ∈ ]0,1] that is used to find the maximum parameter σ of problem (7) in a dichotomy search.

Example:    0.3 means that the first σ tested in the dichotomy search is equal to max(1,Int(Ceil(0.3*|I|))).
+ **LP_time_limit** (Float64): sets the time limit in seconds for each LP subproblem. If a time limit is reached, the execution of this function stops with the error "time_limit atteint dans le modèle LP".

Example:    100.0 sets the time limit of LP solve to 100 seconds.
+ **firstv** (String): sets the vertices of the current domain to be used to produce pieces. Options are "all" and "near60".
  + option "all" produces a piece for each vertex of the current domain, and an evaluation function (parameter eval_f_str) selects the "best piece" to keep among them.
  + option "near60" produces one piece with first vertex the vertex of the current domain with associated angle nearest to 60°.

Example:    "all".
+ **bounding_method** (String): sets the method to build the pieces of a PWL inner corridor. Options are "eff" and "lipschitz2".
  + option "eff" is an iterative subdivision scheme enforcing each piece of the PWL inner corridor to respect a given efficiency (parameter min_eff)
  + option "lipschitz2" ensure a valid PWL inner corridor with a regular grid. WARNING: failures might arise with this option due to the less controlled approximation of the corridor. In such a case, try increasing n_corridor or change this option to "eff".

Example:    "lipschitz2"
+ **min_eff** (Float64): number ∈ ]0,1] giving the efficiency to satisfy if option bounding_method is set to "eff".

Example:    O.95
+ **inclined_bounding** (Bool): sets the orientation of rectangles making up the pieces of the PWL inner corridor. Options are true and false.
  + option true directs the rectangles according to the forward direction
  + option false directs the rectangles according to the x-axis and the y-axis

Example:    true
+ **FDH** (String): stands for Forward Direction Heuristic. It selects the method to chose the forward direction. Options are "bisector_direction" and "mean_along_edge"
  + option "bisector_direction" selects the forward direction as the one following the bisector and pointing at the interior of the current domain
  + option "mean_along_edge", corresponding to "med" in the article, selects the forward direction by taking into account how much the piece can extend along the two edges of the current domain around first vertex
  + option ""

Example:    "mean_along_edge"
+ **eval_f_str** (String): selects the method to choose the "best piece" among a set of piece produced with the option firstv = "all". Options are "Area", "PaD", "SPDTV", "SDI", "SSDI" and "APCI".
  + option "Area" scores the pieces with the area. The bigger the area, the better it is.
  + option "PaD" for "Partial Derivatives total variation" approximates the total variation of the partial derivatives of the function to piecewise linearize. The bigger the area, the better it is.
  + option "SPDTV" uses a convex combination of "Area" and "PaD"
  + option "SDI" stands for "Second Derivative Integral" and approximates the integral of the hessian norm
  + option "SSDI" uses a convex combination of "Area" and "SDI"
  + option "APCI" stands for "Absolute Principal Curvatures Integral" and uses the principal curvatures of the hessian of the function to piecewise linearize

Example: "Area"
+ **compute_piece_version** (String): selects the method to compute a candidate piece, corresponding to the technical way to compute 'key idea 3' from article [1]. Options are "cutting_corridor" and "real_dichotomy". 
  + option "cutting_corridor" computes value σ of optimization problem (7) in [1] among a discrete set of values {σ_k}, where each k corresponds to an element of I∪J
  + option "real_dichotomy" computes value σ of optimization problem (7) in [1] by dichotomy on its bounds. Those bounds are computed from the domain C^d_σ.

Example: "real_dichotomy"
+ **n_eval_sample** (Int64): number of samples used in the scoring function when necessary.

Example:    50
+ **save_solution** (Bool): controls the save of the input and output in a file. Options are true and false.
  + option true adds the input and output of the run at the end of the file "resultats.txt"
  + option false does not save anything

Example:    true
+ **plot_pwl** (Bool): controls if the computed PWL is plotted. Options are true and false.

Example:    false
+ **plot_partial_pwl** (Bool): controls if the PWL built is plotted after each new piece. Options are true and false.

Example:    true


### Output:
+ **polygons** (Vectors{Vector{Vector{Float64}}}): list of pieces of the piecewise linear function g approximating function f on domain domain with approximation error err.
        Each piece is a list of vertices of \\mathbb{R}^3: the first two coordinates are variables X and Y, while the third is the value of g(X,Y).

Example:    [Any[[2.0, 3.61768018018, 8.21726025077892], [2.0, 2.0, 3.014892513924398], [6.853040540541, 2.0, 14.688938483459221]], 
        Any[[8.0, 4.0, 31.011094492736703], [6.304054054054, 4.0, 25.824988820772575], [3.578547297297, 3.091497747748, 11.118987553419728], [6.853040540541, 2.0, 13.47724901735819], [8.0, 2.0, 16.984585422997753]], 
        Any[[6.304054054054, 4.0, 24.262039374051042], [2.0, 4.0, 7.045823157835027], [2.0, 3.61768018018, 6.281183518195025], [3.578547297297, 3.091497747748, 11.543007842519028]]]
    


## Pipeline
To use function PWL2D_heuristic:
+ open a julia REPL in folder PWL_approx_2D_heuristic
+ enter "include("main.jl")" to charge the functions
+ launch function with "PWL2D_heuristic({INSERT ARGUMENTS})"

## Example
```f(x) = x[1]^2-x[2]^2
str_exprf = "X*X-Y*Y"
err = Absolute(1.5)
domain = [0.5,7.5,0.5,3.5]
PWL2D_heuristic(f,str_exprf,err,domain=domain)
```

## Benchmark
function benchmark_flex_heuristic3 from file benchmark_heuristic.jl launches flex_heuristic_procedure on all instances of the benchmark of Rebennack and Kallrath (2014).
Two examples of how to call this function are available on this file.

## Analysis
Most tables analysing the numerical results from this heuristic and the state of the art solvers are computed with functions from file managing_benchmark_results.jl.
See examples lines 1407-1421 from the same file.

## Result files and produce figures of the solutions
Result files of the experiments in [1] and [2] are available in folder "results".
The solutions of the algorithm written in the result files can be displayed using function parse_and_plot_instance from file results/parse_and_plot_solutions.jl. More details in results/README.txt.

Examples:

![Solution with method ProMA_Aire_95_dich of function L1 with delta equal to 1.5](https://github.com/LICO-labs/two-variable-function-PWL-approximation/blob/main/results/improvement%20from%20thesis%20results/real_dichotomy/images/solution_ProMA_VaT_95_dich_L1_delta1.5.png)
