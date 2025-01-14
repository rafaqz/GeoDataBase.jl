# Tracks the burning status for each column
struct BurnStatus
    ic::Int
    burn::Bool
    hasburned::Bool
end
BurnStatus() = BurnStatus(1, false, false)

# _burn_polygon!
# Burn `true` values into a raster
# `boundary` determines how edges are handled 
function _burn_polygon!(B::AbstractDimArray, geom; kw...)::Bool
    B1 = _prepare_for_burning(B)
    _burn_polygon!(B1::AbstractDimArray, GI.geomtrait(geom), geom; kw...)
end
function _burn_polygon!(B::AbstractDimArray, trait, geom;
    fill=true, boundary=:center, geomextent, verbose=false, allocs=Allocs(B), kw...
)::Bool
    allocs = _get_alloc(allocs)
    edges = Edges(geom, dims(B); allocs)
    
    hasburned::Bool = _burn_polygon!(B, edges, allocs.crossings; fill)

    # Lines
    n_on_line = 0
    if boundary !== :center
        _check_intervals(B, boundary, verbose)
        if boundary === :touches 
            if _check_intervals(B, boundary, verbose)
                # Add line pixels
                n_on_line = _burn_lines!(B, geom; fill)::Int
            end
        elseif boundary === :inside 
            if _check_intervals(B, boundary, verbose)
                # Remove line pixels
                n_on_line = _burn_lines!(B, geom; fill=!fill)::Int
            end
        else
            throw(ArgumentError("`boundary` can be :touches, :inside, or :center, got :$boundary"))
        end
        if verbose
            (n_on_line > 0) || @info "$n_on_line pixels were on lines"
        end
    end

    hasburned |= (n_on_line > 0)

    return hasburned
end
function _burn_polygon!(A::AbstractDimArray, edges::Edges, crossings::Vector{Float64};
    offset=nothing, verbose=true, fill=true
)::Bool
    local prev_ypos = 0
    hasburned = false
    # Loop over each index of the y axis
    for iy in axes(A, YDim)
        # Calculate where on the x axis iy is crossed
        ncrossings, prev_ypos = _set_crossings!(crossings, edges, iy, prev_ypos)
        # Burn between alternate crossings
        status = _burn_crossings!(A, crossings, ncrossings, iy; fill)
        hasburned |= status.hasburned
    end
    return hasburned
end

function _set_crossings!(crossings::Vector{Float64}, edges::Edges, iy::Int, prev_ypos::Int)
    # max_ylen tells us how big the largest y edge is.
    # We can use this to jump back from the last y position
    # rather than iterating from the start of the edges
    ypos = max(edges.min_y, prev_ypos - edges.max_ylen - 1)
    start_ypos = searchsortedfirst(edges, ypos)
    # We know the maximum size on y, so we can start from ypos 
    prev_ypos = start_ypos
    ncrossings = 0
    for i in start_ypos:lastindex(edges)
        e = @inbounds edges[i]
        # Edges are sorted on y, so we can skip
        # some at the end once they are larger than iy
        if iy < e.iystart 
            prev_ypos = iy
            break
        end
        if iy <= e.iystop 
            ncrossings += 1
            if ncrossings <= length(crossings)
                @inbounds crossings[ncrossings] = _x_at_y(e, iy)
            else
                push!(crossings, _x_at_y(e, iy))
            end
        end
    end
    # For some reason this is much faster than `partialsort!`
    sort!(view(crossings, 1:ncrossings))
    return ncrossings, prev_ypos
end

function _burn_crossings!(A, crossings, ncrossings, iy; 
    status::BurnStatus=BurnStatus(), fill::Bool=true
) 
    stop = false
    # Start burning loop from outside any rings
    (; ic, burn) = status
    ix = firstindex(A, X())
    hasburned = false
    while ic <= ncrossings
        crossing = crossings[ic]
        # Burn/skip until we hit the next edge crossing
        while ix < crossing
            if ix > lastindex(A, X()) 
                stop = true
                break
            end
            if burn
                @inbounds A[X(ix), Y(iy)] = fill
                hasburned = true
            end
            ix += 1
        end
        if stop
            break
        else
            # Alternate burning/skipping with each edge crossing
            burn = !burn
            ic += 1
        end
    end
    # Maybe fill in the end of the row
    if burn
        for x in ix:lastindex(A, X())
            @inbounds A[X(ix), Y(iy)] = fill
        end
    end
    return BurnStatus(ic, burn, hasburned)
end

const INTERVALS_INFO = "makes more sense on `Intervals` than `Points` and will have more correct results. You can construct dimensions with a `X(values; sampling=Intervals(Center()))` to achieve this"

@noinline _check_intervals(B, verbose::Bool) = 
    _chki(B) ? true : (verbose && @info "burning lines $INTERVALS_INFO"; false)
@noinline _check_intervals(B, boundary, verbose::Bool) =
    _chki(B) ? true : (verbose && @info "`boundary=:$boundary` $INTERVALS_INFO"; false)

_chki(B) = all(map(s -> s isa Intervals, sampling(dims(B)))) 
