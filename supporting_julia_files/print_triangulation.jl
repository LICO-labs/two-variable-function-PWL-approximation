function print_polygon(P, digits=3)
    # retourne un string représentant le polygone P
    s = ""
    for i in 1:length(P)
        S = P[i]
        v1 = round(S[1],digits=digits)
        v2 = round(S[2],digits=digits)
        v3 = round(S[3],digits=digits)
        s = string(s, "[$v1, $v2, $v3] ")
    end
    return s
end

function print_piece(pieces, x, y)
    # affiche la piece qui contient [x,y]
    n = size(pieces)[1]
    for i in 1:n
        x1 = pieces[i,1,1]
        x2 = pieces[i,2,1]
        y1 = pieces[i,1,2]
        y2 = pieces[i,3,2]
        if (x1 <= x <= x2) && (y1 <= y <= y2)
            return pieces[i,:,:]
        end
    end
    return -1
end
