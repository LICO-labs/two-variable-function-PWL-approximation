using LinA

function launch_lin(x1 = -1, str_f = "x*x*x", err = Absolute(0.1), x2max = 1)
    expr_f = Meta.parse(str_f)
    PWL1 = LinA.ExactLin(expr_f,x1,x2max,err)
    for i in 1:length(PWL1)
        println("piece $i: ", PWL1[i])
    end
    return PWL1
end
