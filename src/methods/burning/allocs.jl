struct Allocs{B}
    buffer::B
    edges::Vector{Edge}
    scratch::Vector{Edge}
    crossings::Vector{Float64}
end
function Allocs(buffer)
    edges = Vector{Edge}(undef, 0)
    scratch = Vector{Edge}(undef, 0)
    crossings = Vector{Float64}(undef, 0)
    return Allocs(buffer, edges, scratch, crossings)
end

function _burning_allocs(x; 
    nthreads=_nthreads(), 
    threaded=true, 
    burncheck_metadata=Metadata(),
    kw...
) 
    if threaded
        if isnothing(x)
            [Allocs(nothing) for _ in 1:nthreads]
        else
            dims = commondims(x, DEFAULT_POINT_ORDER)
            [Allocs(_init_bools(dims; metadata=deepcopy(burncheck_metadata))) for _ in 1:nthreads]
        end
    else
        if isnothing(x)
            Allocs()
        else
            dims = commondims(x, DEFAULT_POINT_ORDER)
            Allocs(_init_bools(dims; metadata=burncheck_metadata))
        end
    end
end

_get_alloc(allocs::Vector{<:Allocs}) = _get_alloc(allocs[Threads.threadid()])
_get_alloc(allocs::Allocs) = allocs

# TODO include these in Allocs
_alloc_burnchecks(n::Int) = fill(false, n)
_alloc_burnchecks(x::AbstractArray) = _alloc_burnchecks(length(x))

function _set_burnchecks(burnchecks, metadata::Metadata{<:Any,<:Dict}, verbose)
    metadata["missed_geometries"] = .!burnchecks
    verbose && _burncheck_info(burnchecks)
end
_set_burnchecks(burnchecks, metadata, verbose) = verbose && _burncheck_info(burnchecks)
function _burncheck_info(burnchecks)
    nburned = sum(burnchecks)
    nmissed = length(burnchecks) - nburned
    nmissed > 0 && @info "$nmissed geometries did not affect any pixels. See `metadata(raster)[\"missed_geometries\"]` for a vector of misses"
end
