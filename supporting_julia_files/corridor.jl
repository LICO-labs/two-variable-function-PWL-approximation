using LinA

function get_corridor_functions(f,err,corridor_type = "double")
    # retourne fm et fp les fonctions évaluant respectivement le bas et le
    # haut du corridor de f avec erreur err (qui peut être absolue ou relative mais pas encore mix)
    if typeof(err) == Absolute
        delta = err.delta
        fm = x -> f(x) - delta
        fp = x -> f(x) + delta
	    if corridor_type == "low"
	        fp = x -> f(x)
	    end
	    if corridor_type == "high"
	        fm = x -> f(x)
	    end
    elseif typeof(err) == Relative
        epsilon = err.percent/100
        fm = x -> f(x)*(1-epsilon)
        fp = x -> f(x)*(1+epsilon)
		if corridor_type == "low"
	        fp = x -> f(x)
	    end
	    if corridor_type == "high"
	        fm = x -> f(x)
	    end
    else
        error("unknown ErrorType :\n$err")
    end
    return fp,fm
end

function eval_2DPWL(polygons, x, domain)
	# evaluate at point x the PWL 2D function represented by polygons defined in domain
	n = length(polygons)
	A,B,C,D = domain # domain [A,B]x[C,D]

	for i in 1:n
		P = polygons[i]
		# check if x in domain of P
		bool_inside = true
		m = length(P)
		for j in 1:m
			x1 = P[j][1:2]
			x2 = P[mod(j,m)+1][1:2]
			n,b = line_characteristics(x1,x2)
			if !(scalar_product(n,x) >= b)
				bool_inside = false
				break
			end
		end
		if bool_inside
			println("inside piece $P with point $x and plane")
			a,b,c = get_planar_equation(P)
			println("[$a,$b,$c]")
			return a*x[1]+b*x[2]+c
		end
	end
	return "not found"
end

function check_error_by_grid(f, polygons, err, domain, nx = 100, ny = 100, eps = 1e-9)
	# check if PWL 2D function represented by polygons respect error err in domain domain

	A,B,C,D = domain
	u,l = get_corridor_functions(f,err)
    incr_x = (B-A)/(nx-1)
    incr_y = (D-C)/(ny-1)
	cpt = 0
	L_errors = []
	for i in 0:nx-1
		for j in 0:ny-1
			x = [A+incr_x*i, C+incr_y*j]
			val = eval_2DPWL(polygons, x, domain)
			valu = u(x)
			vall = l(x)
			if vall <= val <= valu
				# error respected at x
				if nx*ny <= 100
					println("$vall <= $val <= $valu is respected")
				end
			else
				println()
				println("error not respected at point $x")
				println("$vall <= $val <= $valu is not respected")
				println("for function\n$polygons\nwith error $err")
				cpt += 1
				push!(L_errors, x)
			end
		end
	end
	println("number of errors: $cpt\nat positions:\n$L_errors")
	return cpt, L_errors
end
