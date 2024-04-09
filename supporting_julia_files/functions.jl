include("test_functions.jl")

function get_limits_and_func(num_func)
    if num_func == 1 || num_func == "L1"
        expression = "X*X-Y*Y"
        A,B,C,D,func = 0.5,7.5,0.5,3.5,f1
    elseif num_func == 2 || num_func == "L2"
        expression = "X*X+Y*Y"
        A,B,C,D,func = 0.5,7.5,0.5,3.5,f2
    elseif num_func == 201
        expression = "X*X+Y*Y"
        A,B,C,D,func = 0,1,0,1,f2
    elseif num_func == 3 || num_func == "N1"
        expression = "X*Y"
        A,B,C,D,func = 2.0,8.0,2.0,4.0,f3
    elseif num_func == 301
        expression = "X*Y"
        A,B,C,D,func = 0.0,1.0,0.0,1.0,f3
    elseif num_func == 302
        expression = "X*Y"
        A,B,C,D,func = 0.0,6.0,0.0,2.0,f3
    elseif num_func == 4 || num_func == "N2"
        expression = "X*exp(-X*X-Y*Y)"
        A,B,C,D,func = 0.5,2.0,0.5,2.0,f4
    elseif num_func == 5 || num_func == "N3"
        expression = "X*sin(Y)"
        A,B,C,D,func = 1.0,4.0,0.05,3.1,f5
    elseif num_func == 6 || num_func == "N4"
        expression = "sin(X)/X*Y*Y"
        A,B,C,D,func = 1.0,3.0,1.0,2.0,f6
    elseif num_func == 7 || num_func == "N5"
        expression = "X*sin(X)*sin(Y)"
        A,B,C,D,func = 0.05,3.1,0.05,3.1,f7
    elseif num_func == 8 || num_func == "N6"
        expression = "(X^2-Y^2)^2"
        A,B,C,D,func = 1.0,2.0,1.0,2.0,f8
    elseif num_func == 9 || num_func == "N7"
        expression = "exp(-10*(X^2-Y^2)^2)"
        A,B,C,D,func = 1.0,2.0,1.0,2.0,f9
    elseif num_func == 21
        expression = "X*X-Y*Y"
        A,B,C,D,func = 0.0,1.0,0.0,1.0,f1
    elseif num_func == 22
        expression = "X*X+Y*Y"
        A,B,C,D,func = 0.0,1.0,0.0,1.0,f2
    elseif num_func == 23
        expression = "X*Y"
        A,B,C,D,func = 0.0,1.0,0.0,1.0,f3
    elseif num_func == 0
        expression = "sqrt(X^2+Y^2)"
        A,B,C,D,func = -1,1,-1,1,f0
    elseif num_func == 10
        expression = "sqrt(X^2+Y^2)"
        A,B,C,D,func = 0.01,1,0.01,1,f0
    elseif num_func == "inverse"
        expression = "X/Y"
        A,B,C,D,func = 100,300,100,300,f_inverse
    elseif num_func == 1011
        expression = "X*X-Y*Y"
        A,B,C,D,func = 0.5,4,0.5,2,f1
    elseif num_func == 1012
        expression = "X*X-Y*Y"
        A,B,C,D,func = 0.5,4,2,3.5,f1
    elseif num_func == 1013
        expression = "X*X-Y*Y"
        A,B,C,D,func = 4,7.5,2,3.5,f1
    elseif num_func == 1014
        expression = "X*X-Y*Y"
        A,B,C,D,func = 4,7.5,0.5,2,f1
    elseif num_func == 2011
        expression = "X*X-Y*Y"
        A,B,C,D,func = 0.5,4,0.5,2,f1
    elseif num_func == 2012
        expression = "X*X-Y*Y"
        A,B,C,D,func = 0.5,4,2,3.5,f1
    elseif num_func == 2013
        expression = "X*X-Y*Y"
        A,B,C,D,func = 4,7.5,2,3.5,f1
    elseif num_func == 2014
        expression = "X*X-Y*Y"
        A,B,C,D,func = 4,7.5,0.5,2,f1
    elseif num_func == "brekelmans" # Brekelmans12 : x*exp(-y)
        expression = "X*exp(-Y)"
        # (c_l y_lk + b_l u_lk = 2000 because c_l may be bounded by 500, b_l by 3 and u_lk by 500
        # lambda_l * sum of u_li = 5 because sum of u_li is bounded by around 500 centimeters (fig 2 Brekelmans12 is 200 centimeters), and lambda may be bounded by 0.01, using the lambda value of Eijgenraam16
        A,B,C,D,func = 1,2000,0,5,f_brekelmans
    elseif num_func == "brekelmans_absolute" # Brekelmans12 : x*exp(-y)
        # defined on [0,1]x[0,5]
        expression = "X*exp(-Y)"
        A,B,C,D,func = 0,1,0,5,f_brekelmans
    end
    return A,B,C,D,func,expression
end
