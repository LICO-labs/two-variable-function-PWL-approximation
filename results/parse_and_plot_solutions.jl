using Pkg
Pkg.add("Plots")
Pkg.add(url="https://github.com/LICO-labs/LinA.jl.git")
using LinA, Plots

SEP1 = "::"
SEP2 = " "
SEP3 = "_"
SEP_MAT = ["___","__","_"]

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

function parse_one_line_result(s)
    """
    Parse the information on one instance.
    INPUT:
    s: string describing one instance (one line, so there is no '\n')
    OUTPUT:
    str_exprf: expression of the function
    domain: domain of the function in the format [A,B,C,D] meaning [A,B]x[C,D]
    err: error allowed for this instance
    n_pieces: number of pieces (not all are triangles) of the PWL function
    t: computing time in seconds
    choices: list of arguments of the function in the format "name_arg1 val_arg1 name_arg2 val_arg2 ..."
    T_set: list of pieces of the PWL function; each piece is described as a set of extreme points (x,y,g(x,y))
    """

    # separator recovery
    sep1 = SEP1
    sep2 = SEP2
    sep3 = SEP3
    sep_mat = SEP_MAT

    # parsing
    elements = split(s,sep1)

    method = elements[1]

    str_exprf = elements[2]

    s_domain = split(elements[3],sep2)
    domain = [parse(Float64,s_domain[1]),parse(Float64,s_domain[2]),parse(Float64,s_domain[3]),parse(Float64,s_domain[4])]

    s_err = split(elements[4],sep2)
    delta,epsilon = [parse(Float64,s_err[1][6:end]),parse(Float64,s_err[2][8:end])]
    if epsilon == 0.0
        err = Absolute(delta)
    elseif delta == 0.0
        err = Relative(epsilon/100)
    else
        error("parsing of the error not handled : $(elements[4])")
    end

    s_n = split(elements[5],sep2)
    n_pieces = parse(Int64, s_n[1])

    s_t = split(elements[6],sep2)
    t = 0
    try
        t = parse(Float64, s_t[1])
    catch e
        t = s_t[1]
    end

    elements[7] = replace(elements[7], ", " => ",")
    s_choices = split(elements[7],sep2)
    if s_choices[end] == ""
        s_choices = s_choices[1:end-1]
    end
    choices = []
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

    return str_exprf, domain, err, n_pieces, t, choices, T_set
end

function simplify_ISCO22_results(file_out, file_to_simplify, method_name)
    """Write in file 'file_out' a simplified version of the results in file 'file_to_simplify' which are of the model of those obtained for the ISCO 2022 article
    the line at the start of file_out stores file_to_simplify for tracking purpose because information is lost
    the only other informations are, in this order:
    - the name of the method name (input)
    - the function name (names: L1,L2,N1,...,N7)
    - the function expression (X*X-Y*Y...)
    - the error allowed
    - the domain
    - the number of pieces
    - the time
    - the solution"""

    lines = readlines(file_to_simplify)
    file = open(file_out,"w")

    # s will be the simplified file
    s = ""

    conv_func = ["L1","L2","N1","N2","N3","N4","N5","N6","N7"]
    cpt = 0
    for line in lines
        cpt += 1
        str_exprf, domain, err, n_pieces, t, choices, T_set = parse_one_line_result(line)
        num_func = mod(cpt,9)
        if num_func == 0
            num_func = 9
        end
        s = string(s,method_name,SEP1,conv_func[num_func],SEP1,str_exprf,SEP1,err,SEP1,domain,SEP1,n_pieces,SEP1,t,SEP1,T_set,"\n")
    end
    print(file,s)

    close(file)
end

function list_to_array(l)
    # return a matrix of size (length(l),length(l[1])) from a list l of elements of \mathbb{R}^2
    # this function is used to display polygons

    M = zeros(length(l),length(l[1]))

    for i in 1:length(l)
        M[i,:] = [l[i][j] for j in 1:length(l[1])]
    end
    return M
end

function plot_polygon(p, P, polygon_color = :blue, draw_inside = false, transparency = 1,label_poly = "")
    # return plot p with a polygon P of color polygon_color
    PP = copy(P)
    push!(PP,P[1])
    global mat_P
    try
        global mat_P = list_to_array(PP)
    catch e
        println(e)
        println(PP)
    end
    if draw_inside
        p = plot!(mat_P[:,1], mat_P[:,2], color = polygon_color, seriestype = :shape, seriesalpha = transparency, label = label_poly)
    else
        p = plot!(mat_P[:,1], mat_P[:,2], color = polygon_color, label = label_poly)
    end
    return p
end

function plot_pieces(P, p = plot(), pieces_color = :blue, draw_inside = false)
    """return plot p with all the pieces in P in the color pieces_color

    EXAMPLE OF CALL:
    P = [[[3.0, 1.75, 0.003410107400922513], [3.0, 2.0, 0.4900990588365719], [1.0, 2.0, 2.9196613880744238], [1.0, 1.0, 0.9729055823318264], [1.5, 1.0, 0.3655150000223635]], [[3.0, 1.75, 0.47862098270737996], [1.5, 1.0, 0.198614808630589], [3.0, 1.0, -0.447522220363739]]]
    p = plot(title="Subdomains constituting a PWL function") # choose the title of the plot
    pieces_color = :green
    draw_inside = true
    p = plot_pieces(P, p, pieces_color, draw_inside)
    savefig("images/file_name_for_the_plot.png")
    display(p)"""

    p = plot!(legend = false)
    try
        for i in 1:size(P)[1]
            piece = P[i,:,:]
            piece = [piece[1,:],piece[2,:],piece[3,:],piece[4,:]]
            p = plot_polygon(p, piece, pieces_color, draw_inside)
        end
    catch e
        for i in 1:length(P)
            piece = P[i]
            p = plot_polygon(p, piece, pieces_color, draw_inside)
        end
    end
    return p
end

function parse_and_plot_instance(method_name, function_name, delta)
    """
    Plot the piecewise linear function solution of instance (function_name, delta) with method method_name
    and return the informations of the instance
    INPUT:
    method_name (string): name of the method used for the resolution of the CFP instance. Example "DN95".
    function_name (string): name of the function as in the article (Ref in the article). "L2", "N5"...
    delta (Float64): absolute approximation error allowed in the CFP instance (delta in the article)
    name_figure: string with the name under which a plot of the solution will be saved
    OUTPUT:
    T_set: piecewise linear function solution of the instance in a format handled by plot_pieces (used to plot the solution)
    domain: domain of the piecewise linear function in the format [A,B,C,D] for [A,B]x[C,D] \\in \\mathbb{R}^2
    delta: absolute approximation error allowed for this instance
    n_pieces: number of pieces in the piecewise linear function solution
    cpu_time: time in seconds needed by the algorithm to solve the instance
    str_exprf: string containing the expression of the function function_name
    EXAMPLE OF CALL:
    T_set, domain, delta, n_pieces, cpu_time, str_exprf = parse_and_plot_instance("DN95", "N1", 0.05)
    """

    # read lines of the .txt file with the solutions
    filename = method_name * ".txt"
    lines = readlines(filename)

    name_figure = "solution_"*method_name*"_"*function_name*"_delta"*string(delta)*".png"

    # create dictionary with the function names employed in the .txt files as values and function_name as keys
    real_func_names = Dict("L1"=>"num_func 1","L2"=>"num_func 2","N1"=>"num_func 3","N2"=>"num_func 4","N3"=>"num_func 5","N4"=>"num_func 6","N5"=>"num_func 7","N6"=>"num_func 8","N7"=>"num_func 9",)

    # find the line of instance (function_name, delta)
    global str_exprf, domain, err, n_pieces, cpu_time, choices, T_set
    line_found = false
    cpt = 0
    for line in lines
        cpt += 1
        real_func_name = real_func_names[function_name]
        res1 = findfirst(real_func_name, line)
        res2 = findfirst("delta$delta", line)
        if res1 != nothing && res2 != nothing
            str_exprf, domain, err, n_pieces, cpu_time, choices, T_set = parse_one_line_result(line)
            delta = err.delta
            line_found = true
            break
        end
    end

    if line_found
        # plot the piecewise linear function solution of the instance
        p = plot(title = "solution with $method_name of instance ($function_name, $delta), $n_pieces pieces")
        p = plot_pieces(T_set, p)
        print("\n save solution figure ''", name_figure, "'' in folder images/ \n\n" )
        if !("images" in readdir())
            run(`mkdir images`)
        end
        savefig("images/"*name_figure) # save figure with name name_figure
        display(p)
    else
        error("instance ($function_name, $delta) with method $method_name was not found")
    end

    return T_set, domain, delta, n_pieces, cpu_time, str_exprf
end

#T_set, domain, delta, n_pieces, cpu_time, str_exprf = parse_and_plot_instance("DN99", "N1", 0.05)
