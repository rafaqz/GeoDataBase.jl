"""
    GeometryLookup(data, dims = (X(), Y()); geometrycolumn = nothing)


A lookup type for geometry dimensions in vector data cubes.

`GeometryLookup` provides efficient spatial indexing and lookup for geometries using an STRtree (Sort-Tile-Recursive tree).
It is used as the lookup type for geometry dimensions in vector data cubes, enabling fast spatial queries and operations.

It spans the dimensions given to it in `dims`, as well as the dimension it's wrapped in - you would construct a DimArray with a GeometryLookup
like `DimArray(data, Geometry(GeometryLookup(data, dims)))`.  Here, `Geometry` is a dimension - but selectors in X and Y will also eventually work!


# Examples

```julia
using Rasters

using NaturalEarth
import GeometryOps as GO

# construct the polygon lookup
polygons = NaturalEarth.naturalearth("admin_0_countries", 110).geometry
polygon_lookup = GeometryLookup(polygons, (X(), Y()))

# create a DimArray with the polygon lookup
dv = rand(Geometry(polygon_lookup))

# select the polygon with the centroid of the 88th polygon
dv[Geometry(Contains(GO.centroid(polygons[88])))] == dv[Geometry(88)] # true
```

"""
struct GeometryLookup{T, D, CRS} <: Lookups.Lookup{T, 1}
    data::Vector{T}
    tree::SortTileRecursiveTree.STRtree
    dims::D
    crs::CRS
end

function GeometryLookup(data, dims = (X(), Y()); geometrycolumn = nothing, crs = nothing)
    # First, retrieve the geometries - from a table, vector of geometries, etc.
    geometries = _get_geometries(data, geometrycolumn)
    # Check that the geometries are actually geometries
    if any(!GI.isgeometry, geometries)
        throw(ArgumentError("""
        The collection passed in to `GeometryLookup` has some elements that are not geometries 
        (`GI.isgeometry(x) == false` for some `x` in `data`).
        """))
    end
    # Make sure there are only two dimensions
    if length(dims) != 2
        throw(ArgumentError("""
        The `dims` argument to `GeometryLookup` must have two dimensions, but it has $(length(dims)) dimensions (`$(dims)`).
        Please make sure that it has only two dimensions, like `(X(), Y())`.
        """))
    end
    # Build the lookup accelerator tree
    tree = SortTileRecursiveTree.STRtree(geometries)

    true_crs = isnothing(crs) ? GI.crs(first(geometries)) : crs
        

    GeometryLookup(geometries, tree, dims, crs)
end

crs(l::GeometryLookup) = l.crs

#=

## DD methods for the lookup

Here we define DimensionalData's methods for the lookup.
This is broadly standard except for the `rebuild` method, which is used to update the tree accelerator when the data changes.
=#

DD.dims(l::GeometryLookup) = l.dims

DD.dims(d::DD.Dimension{<: GeometryLookup}) = val(d).dims
 
DD.order(::GeometryLookup) = Lookups.Unordered()

DD.parent(lookup::GeometryLookup) = lookup.data

DD.Dimensions.format(l::GeometryLookup, D::Type, values, axis::AbstractRange) = l

# Make sure that the tree is rebuilt if the data changes
function DD.rebuild(lookup::GeometryLookup; data = lookup.data, tree = nokw, dims = lookup.dims, crs = nokw)
    new_tree = if data == lookup.data
        lookup.tree
    else
        SortTileRecursiveTree.STRtree(data)
    end
    new_crs = if isnokw(crs)
        GI.crs(first(data))
    else
        crs
    end

    GeometryLookup(data, new_tree, dims, new_crs)
end

#=
## Lookups methods

Here we define the methods for the Lookups API.
The main entry point is `selectindices`, which is used to select the indices of the geometries that match the selector.

We need to define methods that take selectors and convert them to extents, then GeometryOps needs 
=#

# Return an `Int` or  Vector{Bool}
Lookups.selectindices(lookup::GeometryLookup, sel::DD.DimTuple) =
    selectindices(lookup, map(_val_or_nothing, sortdims(sel, dims(lookup))))
function Lookups.selectindices(lookup::GeometryLookup, sel::NamedTuple{K}) where K
    dimsel = map(rebuild, map(name2dim, K), values(sel))
    selectindices(lookup, dimsel) 
end
Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.StandardIndices) = sel
function Lookups.selectindices(lookup::GeometryLookup, sel::Tuple)
    if (length(sel) == length(dims(lookup))) && all(map(s -> s isa At, sel))
        i = findfirst(x -> all(map(_matches, sel, x)), lookup)
        isnothing(i) && _coord_not_found_error(sel)
        return i
    else
        return [_matches(sel, x) for x in lookup]
    end
end

"""
    _maybe_get_candidates(tree, selector_extent)

Get the candidates for the selector extent.  If the selector extent is disjoint from the tree rootnode extent, return an error.
"""
function _maybe_get_candidates(tree, selector_extent)
    if Extents.disjoint(tree.rootnode.extent, selector_extent)
        return error("""
            The geometry with extent $(GI.extent(val(sel))) is outside of the extent of the geometry lookup.

            The geometry lookup has extent $(GI.extent(tree.rootnode.extent))
        """)
    end
    potential_candidates = SortTileRecursiveTree.query(tree, selector_extent)

    isempty(potential_candidates) && return error("""
        The geometry with extent $(GI.extent(val(sel))) does not interact with any of the geometries in the lookup.
    """)

    return potential_candidates
end

function Lookups.selectindices(lookup::GeometryLookup, sel::Contains)

    potential_candidates = _maybe_get_candidates(lookup.tree, sel_ext)

    for candidate in potential_candidates
        if GO.contains(lookup.data[candidate], val(sel))
            return candidate
        end
    end
    return error("""
        The geometry with extent $(GI.extent(val(sel))) is not contained by any of the geometries in the lookup.
    """)
end


function Lookups.selectindices(lookup::GeometryLookup, sel::Touches)
    lookup_ext = lookup.tree.rootnode.extent
    sel_ext = GI.extent(val(sel))
    potential_candidates = _maybe_get_candidates(lookup.tree, sel_ext)

    for candidate in potential_candidates
        if GO.intersects(lookup.data[candidate], val(sel))
            return candidate
        end
    end
    return error("""
        The geometry with extent $(GI.extent(val(sel))) does not touch any of the geometries in the lookup.
    """)
end

function Lookups.selectindices(lookup::GeometryLookup, (xs, ys)::Tuple{Union{<: At, <: Contains}, Union{<: At, <: Contains}})
    xval, yval = val(xs), val(ys)

    lookup_ext = lookup.tree.rootnode.extent

    if lookup_ext.X[1] <= xval <= lookup_ext.X[2] && lookup_ext.Y[1] <= yval <= lookup_ext.Y[2]
        potential_candidates = SortTileRecursiveTree.query(lookup.tree, (xval, yval))
        if isempty(potential_candidates)
            return error("""
            The point ($xval, $yval) is not within any of the geometries in the lookup.
        """) # no geometry intersects with it
        else
            for candidate in potential_candidates
                if GO.contains(lookup.data[candidate], (xval, yval))
                    return candidate
                end
            end
            return error("""
                The point ($xval, $yval) is not within any of the geometries in the lookup.
            """) # no geometry intersects with it
        end
    else
        return  error("""
        The point ($xval, $yval) is outside of the extent of the geometry lookup. 
        """) # outside of extent / geometrylookup bbox
    end
end

@inline Lookups.reducelookup(l::GeometryLookup) = NoLookup(OneTo(1))

function Lookups.show_properties(io::IO, mime, lookup::GeometryLookup)
    print(io, " ")
    show(IOContext(io, :inset => "", :dimcolor => 244), mime, DD.basedims(lookup))
end

# Dimension methods

@inline _reducedims(lookup::GeometryLookup, dim::DD.Dimension) =
    rebuild(dim, [map(x -> zero(x), dim.val[1])])

function DD._format(dim::DD.Dimension{<:GeometryLookup}, axis::AbstractRange)
    checkaxis(dim, axis)
    return dim
end

# Local functions
_val_or_nothing(::Nothing) = nothing
_val_or_nothing(d::DD.Dimension) = val(d)


#=
## Reproject

Reproject just forwards to `GO.reproject`.
=#

function reproject(target::GeoFormat, l::GeometryLookup)
    source = crs(l)
    return rebuild(l; data = GO.reproject(l.data; source_crs = source, target_crs = target, always_xy = true))
end

# I/O utils
function _geometry_cf_encode(::Union{GI.PolygonTrait, GI.MultiPolygonTrait}, geoms)
    n_points_per_geom_vec = GI.npoint.(geoms)
    total_n_points = sum(n_points_per_geom_vec)

    # Create a vector of the total number of points
    xs = fill(0.0, total_n_points)
    ys = fill(0.0, total_n_points)

    ngeoms = length(geoms)
    nrings = GO.applyreduce(GI.nring, +, GI.PolygonTrait(), geoms; init = 0)

    node_count_vec = fill(0, ngeoms)
    part_node_count_vec = fill(0, nrings)
    interior_ring_vec = fill(0, nrings)

    current_xy_index = 1
    current_ring_index = 1

    for (i, geom) in enumerate(geoms)

        this_geom_npoints = GI.npoint(geom)
        node_count_vec[i] = this_geom_npoints

        # push individual components of the ring
        for poly in GO.flatten(GI.PolygonTrait, geom)
            exterior_ring = GI.getexterior(poly)
            for point in GI.getpoint(exterior_ring)
                xs[current_xy_index] = GI.x(point)
                ys[current_xy_index] = GI.y(point)
                current_xy_index += 1
            end
            part_node_count_vec[current_ring_index] = GI.npoint(exterior_ring)
            interior_ring_vec[current_ring_index] = 0
            current_ring_index += 1

            if GI.nring(poly) == 1
                continue
            else
                for hole in GI.gethole(poly)
                    for point in GI.getpoint(hole)
                        xs[current_xy_index] = GI.x(point)
                        ys[current_xy_index] = GI.y(point)
                        current_xy_index += 1
                    end
                    part_node_count_vec[current_ring_index] = GI.npoint(hole)
                    interior_ring_vec[current_ring_index] = 1
                    current_ring_index += 1
                end
            end
        end
    end

    return xs, ys, node_count_vec, part_node_count_vec, interior_ring_vec

end

function _geometry_cf_encode(::Union{GI.LineStringTrait, GI.MultiLineStringTrait}, geoms)
    error("Not implemented yet")
end

function _geometry_cf_encode(::Union{GI.PointTrait, GI.MultiPointTrait}, geoms)
    error("Not implemented yet")
end


function _geometry_cf_decode(::Union{GI.PolygonTrait, GI.MultiPolygonTrait}, xs, ys, node_count_vec, part_node_count_vec, interior_ring_vec)
    current_xy_index = 1
    current_ring_index = 1
    current_geom_index = 1

    geoms = Vector{GO.GeometryBasics.MultiPolygon{2, Float64}}(undef, length(node_count_vec))

    for (i, npoints) in enumerate(node_count_vec)

        this_geom_npoints = npoints
        # this_geom_nholes = 
    end
end


# total_area_of_intersection = 0.0
# current_area_of_intersection = 0.0
# last_point = nothing
# apply_with_signal(trait, geom) do subgeom, state
#     if state == :start
#         total_area_of_intersection += current_area_of_intersection
#         current_area_of_intersection = 0.0
#         last_point = nothing
#     elseif state == :continue
#         # shoelace formula for this point
#     elseif state == :end
#         # finish off the shoelace formula
#     end
# end