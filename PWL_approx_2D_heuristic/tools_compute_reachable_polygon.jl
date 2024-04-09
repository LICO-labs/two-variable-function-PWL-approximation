function get_polygon_diameter(polygon)
    # mesure du diamètre du polygone polygon avec le max de la distance entre deux de ses sommets
    n = length(polygon)
    max_dist = 0
    for i in 1:(n-1)
        for j in (i+1):n
            distance = dist2(polygon[i],polygon[j])
            if max_dist < distance
                max_dist = distance
            end
        end
    end
    return max_dist
end
