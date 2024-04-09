function mat_to_txt(mat, filename)
    # convertit une matrice mat de taille quelconque (dimension <= 5) en un string dans lequel les éléments de mat
    # sont séparés par un espace et sont passés dans l'ordre lexicographique
    N = size(mat)
    n = length(N)
    if n == 1
        s = ""
        for i in 1:N[1]
            s = string(s,"$mat[i] ")
        end
        s = s[1:(length(s)-1)]
    elseif n == 2
        s = ""
        for i1 in 1:N[1]
            for i2 in 1:N[2]
                s = string(s,"$mat[i1,i2] ")
            end
        end
        s = s[1:(length(s)-1)]
    elseif n == 3
        s = ""
        for i in 1:N[1]
            s = string(s,"$mat[i] ")
        end
        s = s[1:(length(s)-1)]
    elseif n == 4
        s = ""
        for i in 1:N[1]
            s = string(s,"$mat[i] ")
        end
        s = s[1:(length(s)-1)]
    elseif n == 5
        s = ""
        for i in 1:N[1]
            s = string(s,"$mat[i] ")
        end
        s = s[1:(length(s)-1)]
    end
    file = open(filename,"w")
    file.write(s)
    file.close()
end

function print_to_file()
    # test pour output dans un fichier
    file = open("heuristique triangle+rectangle/log.txt", "w")
    println(file,"salut1")
    close(file)
    for i in 1:100
        file = open("heuristique triangle+rectangle/log.txt", "a")
        println(file,i)
        close(file)
        sleep(0.2)
    end
end
