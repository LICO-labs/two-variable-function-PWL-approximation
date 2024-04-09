function approximate_1D_corridor(f, delta, X1,X2, n, corridor_type = "double")
	# approxime le corridor 1D de f avec erreur absolue delta de t1 à t2 avec le type de corridor corridor_type
    # donne un encadrement linéaire entre f-delta et f (donc création de la pwl inférieur l) si corridor_type == "low"
    # f et f+delta si corridor_type == "high"
    # f-delta et f+delta si corridor_type == "double"
    # définition de fm et fp les corridors bas et haut respectivement, puis de leur dérivées partielles dfxm, dfym, dfxp, dfyp
    fm = x -> f(x) - delta
    fp = x -> f(x) + delta
    if corridor_type == "low"
        fp = x -> f(x)
    end
    if corridor_type == "high"
        fm = x -> f(x)
    end
    samples = zeros(n,4)
    for i in 1:n
		x = X1 .+ (i-1)/(n-1) .* X2
		samples[i,:] = [x[1],x[2],fm(x),fp(x)]
    end
    return samples
end
