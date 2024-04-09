using Polyhedra

struct CutoffPointIndex
    cutoff::Int
    index::Int
end
Base.show(io::IO, p::CutoffPointIndex) = print(io, "p[(p.index)]")

struct CutoffRayIndex
    cutoff::Int
    index::Int
end
Base.show(io::IO, r::CutoffRayIndex) = print(io, "r[(r.index)]")

struct DoubleDescriptionData{PointT, RayT, LineT, HST}
    fulldim::Int
    halfspaces::Vector{HST}
    # Elements ordered by first halfspace cutting it off
    points::Vector{PointT}
    pz::Vector{BitSet}
    cutpoints::Vector{Vector{PointT}}
    cutpz::Vector{Vector{BitSet}}
    pin::Vector{Vector{CutoffPointIndex}}
    rays::Vector{RayT}
    rz::Vector{BitSet}
    cutrays::Vector{Vector{RayT}}
    cutrz::Vector{Vector{BitSet}}
    rin::Vector{Vector{CutoffRayIndex}}
    lines::Vector{LineT}
    cutline::Vector{Union{Nothing, LineT}}
    lineray::Vector{Union{Nothing, CutoffRayIndex}}
    nlines::Vector{Int}
end

function Base.show(io::IO, data::DoubleDescriptionData)
    println(io, "DoubleDescriptionData in $(data.fulldim) dimension:")
    println(io, data.points)
    println(io, data.rays)
    println(io, data.lines)
    for i in reverse(eachindex(data.cutpoints))
        println(io, " Halfspace (data.halfspaces[i]):")
        if !isempty(data.cutpoints[i])
            println(io, "  Cut points:")
            for j in eachindex(data.cutpoints[i])
                println(io, "  $j: ", data.cutpoints[i][j], " zero at: ", data.cutpz[i][j])
            end
        end
        if !isempty(data.pin[i])
            println(io, "  In: ", data.pin[i])
        end
        if !isempty(data.cutrays[i])
            println(io, "  Cut rays:")
            for j in eachindex(data.cutrays[i])
                println(io, "  $j: ", data.cutrays[i][j], " zero at: ", data.cutrz[i][j])
            end
        end
        if !isempty(data.rin[i])
            println(io, "  In: ", data.rin[i])
        end
        if data.cutline[i] !== nothing
            println(io, "  Cut line: ", data.cutline[i])
            if data.lineray[i] !== nothing
                println(io, "  Line ray: ", data.lineray[i])
            end
        end
        if !iszero(data.nlines[i])
            println(io, "  $(data.nlines[i]) uncut lines left")
        end
    end
end

function DoubleDescriptionData{PointT, RayT, LineT}(fulldim::Integer, halfspaces) where {PointT, RayT, LineT}
    n = length(halfspaces)
    return DoubleDescriptionData{PointT, RayT, LineT, eltype(halfspaces)}(
        fulldim,
        halfspaces,
        PointT[],
        BitSet[],
        [PointT[] for i in 1:n],
        [BitSet[] for i in 1:n],
        [CutoffPointIndex[] for i in 1:n],
        RayT[],
        BitSet[],
        [RayT[] for i in 1:n],
        [BitSet[] for i in 1:n],
        [CutoffRayIndex[] for i in 1:n],
        LineT[],
        Union{Nothing, LineT}[nothing for i in 1:n],
        Union{Nothing, CutoffRayIndex}[nothing for i in 1:n],
        zeros(Int, n)
    )
end

function Base.getindex(data::DoubleDescriptionData, p::CutoffPointIndex)
    if p.cutoff == 0
        return data.points[p.index]
    else
        return data.cutpoints[p.cutoff][p.index]
    end
end

function Base.getindex(data::DoubleDescriptionData, r::CutoffRayIndex)
    if r.cutoff == 0
        return data.rays[r.index]
    else
        return data.cutrays[r.cutoff][r.index]
    end
end

function _bitdot_range(b1::BitSet, b2::BitSet, i, n)
    count = 1 # They share the hyperplance `i`
    for j in (i+1):n
        if j in b1 && j in b2
            count += 1
        end
    end
    return count
end

function isadjacent(data, i::Integer, p1::CutoffPointIndex, p2::CutoffPointIndex)
    pz1 = data.cutpz[p1.cutoff][p1.index]
    pz2 = data.cutpz[p2.cutoff][p2.index]
    n = length(data.halfspaces)
    return _bitdot_range(pz1, pz2, i, n) + data.lines[i] + 1 == data.fulldim
end

function isadjacent(data, i::Integer, p::CutoffPointIndex, r::CutoffRayIndex)
    pz = data.cutpz[p.cutoff][p.index]
    rz = data.cutrz[r.cutoff][r.index]
    n = length(data.halfspaces)
    return _bitdot_range(pz, rz, i, n) + data.lines[i] + 1 == data.fulldim
end

function isadjacent(data, i::Integer, r::CutoffRayIndex, p::CutoffPointIndex)
    return isadjacent(data, i, p, r)
end

function isadjacent(data, i::Integer, r1::CutoffPointIndex, r2::CutoffPointIndex)
    rz1 = data.cutrz[r1.cutoff][r1.index]
    rz2 = data.cutrz[r2.cutoff][r2.index]
    n = length(data.halfspaces)
    return _bitdot_range(rz1, rz2, i, n) + data.lines[i] + 2 == data.fulldim
end

function isin(data, i, p::CutoffPointIndex)
    if p.cutoff == 0
        return i in data.pz[p.index]
    else
        return i in data.cutpz[p.cutoff][p.index]
    end
end

function isin(data, i, r::CutoffRayIndex)
    if r.cutoff == 0
        return i in data.rz[r.index]
    else
        return i in data.cutrz[r.cutoff][r.index]
    end
end

resized_bitset(data) = sizehint!(BitSet(), length(data.halfspaces))

function add_index!(data, cutoff::Nothing, p::AbstractVector)
    push!(data.points, p)
    push!(data.pz, resized_bitset(data))
    return CutoffPointIndex(0, length(data.points))
end

function add_index!(data, cutoff::Integer, p::AbstractVector)
    push!(data.cutpoints[cutoff], p)
    push!(data.cutpz[cutoff], resized_bitset(data))
    return CutoffPointIndex(cutoff, length(data.cutpoints[cutoff]))
end

function add_index!(data, cutoff::Nothing, r::Polyhedra.Ray)
    push!(data.rays, r)
    push!(data.rz, resized_bitset(data))
    return CutoffRayIndex(0, length(data.rays))
end

function add_index!(data, cutoff::Integer, r::Polyhedra.Ray)
    push!(data.cutrays[cutoff], r)
    push!(data.cutrz[cutoff], resized_bitset(data))
    return CutoffRayIndex(cutoff, length(data.cutrays[cutoff]))
end

function add_in!(data, i, index::CutoffPointIndex)
    push!(data.pin[i], index)
    if index.cutoff == 0
        push!(data.pz[index.index], i)
    else
        push!(data.cutpz[index.cutoff][index.index], i)
    end
end

function add_in!(data, i, index::CutoffRayIndex)
    push!(data.rin[i], index)
    if index.cutoff == 0
        push!(data.rz[index.index], i)
    else
        push!(data.cutrz[index.cutoff][index.index], i)
    end
end

function set_in!(data, I, el, index)
    for i in I
        if el in Polyhedra.hyperplane(data.halfspaces[i])
            add_in!(data, i, index)
        end
    end
end

function add_element!(data, k, el)
    cutoff = nothing
    for i in reverse(1:k)
        if data.cutline[i] !== nothing
            el = line_project(el, data.cutline[i], data.halfspaces[i])
            index = add_adjacent_element!(data, i - 1, el, data.lineray[i])
            set_in!(data, i:k, el, index)
            return index
        end
        if !(el in data.halfspaces[i])
            cutoff = i
            break
        end
    end
    index = add_index!(data, cutoff, el)
    set_in!(data, (index.cutoff+1):k, el, index)
    return index
end

function add_adjacent_element!(data, k, el, parent)
    index = add_element!(data, k, el)
    # Condition (c_k) in [FP96]
    if index.cutoff != parent.cutoff
        addintersection!(data, index, parent)
    end
    return index
end

using LinearAlgebra

function combine(β, p1::AbstractVector, value1, p2::AbstractVector, value2)
    λ = (value2 - β) / (value2 - value1)
    return λ * p1 + (1 - λ) * p2
end

function combine(β, p::AbstractVector, pvalue, r::Polyhedra.Ray, rvalue)
    λ = (β - pvalue) / rvalue
    return p + λ * r
end

combine(β, r::Polyhedra.Ray, rvalue, p::AbstractVector, pvalue) = combine(β, p, pvalue, r, rvalue)

function combine(r1::Polyhedra.Ray, value1, r2::Polyhedra.Ray, value2)
    # should take
    # λ = value2 / (value2 - value1)
    @assert 0 <= value2 / (value2 - value1) <= 1
    # By homogeneity we can avoid the division and do
    #newr = value2 * r1 - value1 * r2
    # but this can generate very large numbers (see JuliaPolyhedra/Polyhedra.jl#48)
    # so we still divide
    newr = (value2 * r1 - value1 * r2) / (r2 - r1)
    # In CDD, it does value2 * r1 - value1 * r2 but then it normalize the ray
    # by dividing it by its smallest nonzero entry (see dd_CreateNewRay)
    return Polyhedra.simplify(newr)
end

function addintersection!(data, idx1, idx2)
    if idx1.cutoff > idx2.cutoff
        return addintersection!(data, idx2, idx1)
    end
    i = idx2.cutoff
    if isin(data, i, idx1)
        return
    end
    el1 = data[idx1]
    el2 = data[idx2]
    h = data.halfspaces[i]
    newel = combine(h.β, el1, h.a ⋅ el1, el2, h.a ⋅ el2)
    add_adjacent_element!(data, i - 1, newel, idx1)
end

function add_if_adjacent!(data, i::Integer, el1, el2)
    # Condition (c_k) in [FP96]
    if el1.cutoff != el2.cutoff
        if isadjacent(data, i, el1, el2)
            addintersection!(data, el1, el2)
        end
    end
end

_shift(el::AbstractVector, line::Polyhedra.Line) = el + Polyhedra.coord(line)
_shift(el::Polyhedra.Line, line::Polyhedra.Line) = el + line
_shift(el::Ray, line::Polyhedra.Line) = el + Polyhedra.Ray(Polyhedra.coord(line))

function line_project(el, line, h)
    # (line + λ * cutline) ⋅ h.a == h.β
    # λ = (h.β - line ⋅ h.a) / (cutline ⋅ h.a)
    λ = (h.β - el ⋅ h.a) / (line ⋅ h.a)
    return Polyhedra.simplify(_shift(el, λ * line))
end

function hline(data, line::Polyhedra.Line, i, h)
    value = h.a ⋅ line
    if !Polyhedra.isapproxzero(value)
        if data.cutline[i] === nothing
            if value > 0
                line = -line # Make `lineray` point inward
            end
            data.cutline[i] = line
            cut = true
            return true
        else
            line = line_project(line, data.cutline[i], hs)
        end
    end
    data.nlines[i] += 1
    return false
end

function double_description(hr::HRepresentation)
    v = Polyhedra.dualfullspace(hr)
    hps = Polyhedra.lazy_collect(hyperplanes(hr))
    hss = Polyhedra.lazy_collect(halfspaces(hr))
    data = DoubleDescriptionData{pointtype(v), raytype(v), linetype(v)}(fulldim(hr), hps, hss)
    for line in lines(v)
        cut = false
        for i in reverse(eachindex(hps))
            cut = hline(data, line, nhalfspaces(hr) + i, hss[i])
            if cut
                break
            end
        end
        if !cut
            for i in reverse(eachindex(hss))
                cut = hline(data, line, i, hss[i])
            end
            if cut
                break
            end
        end
        if !cut
            push!(data.lines, line)
        end
    end
    # Add line rays after all lines are added so that the rays can be `line_project`ed.
    # We only do that for halfspaces, hyperplanes do not create rays from cutoff lines.
    # We use increasing index order since higher index may need the `lineray` of lower index.
    for i in eachindex(hss)
        line = data.cutline[i]
        if line !== nothing
            ray = Polyhedra.Ray(Polyhedra.coord(line))
            data.lineray[i] = add_element!(data, i - 1, ray)
        end
    end
    @assert isone(npoints(v))
    add_element!(data, nhalfspaces(hr), first(points(v))) # Add the origin
    for i in reverse(eachindex(hss))
        if isempty(data.cutpoints[i]) && isempty(data.cutrays[i])
            # Redundant, remove its contribution to avoid incorrect `isadjacent`
            for p in data.pin
                if p.cutoff == 0
                    delete!(data.pz, i)
                else
                    delete!(data.cutpz[p.cutoff], i)
                end
            end
            for r in data.rin
                if r.cutoff == 0
                    delete!(data.rz, i)
                else
                    delete!(data.cutrz[pr.cutoff], i)
                end
            end
            continue
        end
        if i > 1
            for p1 in data.pin[i], p2 in data.pin[i]
                add_if_adjacent!(data, i, p1, p2)
            end
            for p in data.pin[i], r in data.rin[i]
                add_if_adjacent!(data, i, p, r)
            end
        end
        deleteat!(data.cutpoints, i)
        deleteat!(data.cutpz, i)
        if i > 1
            for r1 in data.rin[i], r2 in data.rin[i]
                add_if_adjacent!(data, i, r1, r2)
            end
        end
        deleteat!(data.cutrays, i)
        deleteat!(data.cutrz, i)
        deleteat!(data.pin, i)
        deleteat!(data.rin, i)
    end
    if isempty(data.points)
        # Empty polyhedron, there may be rays left,
        # Example 1: for 0x_1 + x_2 = -1 ∩ 0x_1 + x_2 = 1, the line (0, 1) is detected as correct
        # Example 2: for 0x_1 + 0x_2 = 1, the lines (1, 0) and (0, 1) are detected as correct
        # but since there is no point, the polyhedron is empty and we should drop all rays/lines
        empty!(data.lines)
        empty!(data.rays)
    end
    similar(v, data.points, data.lines, data.rays)
end
