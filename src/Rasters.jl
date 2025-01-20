module Rasters

# Use the README as the module docs
@doc let
    # path = joinpath(dirname(@__DIR__), "README.md")
    # include_dependency(path)
    # read(path, String)
end Rasters

using Dates

# Load first to fix StaticArrays invalidations
import DimensionalData

import Adapt,
       ColorTypes,
       CommonDataModel,
       ConstructionBase,
       DiskArrays,
       Extents,
       FillArrays,
       Flatten,
       GeoInterface,
       GeometryOps,
       GeometryOpsCore,
       OffsetArrays,
       ProgressMeter,
       Missings,
       Mmap,
       RecipesBase,
       Reexport,
       Setfield,
       SortTileRecursiveTree,
       Statistics

Reexport.@reexport using DimensionalData, GeoFormatTypes

using DimensionalData.Tables,
      DimensionalData.Lookups,
      DimensionalData.Dimensions,
      DimensionalData.Lookups.IntervalSets

using DimensionalData: Name, NoName
using .Dimensions: StandardIndices, DimTuple
using .Lookups: LookupTuple

using Statistics: mean
using RecipesBase: @recipe, @series
using Base: tail, @propagate_inbounds

import GeoInterface: crs
import Extents: Extent, extent

using Setfield: @set, @set!
using ColorTypes: RGB

using CommonDataModel: AbstractDataset, AbstractVariable

using DiskArrays: @implement_diskarray

using GeometryOpsCore: Planar, Spherical
export Planar, Spherical

export AbstractRaster, Raster
export AbstractRasterStack, RasterStack
export AbstractRasterSeries, RasterSeries
export Projected, Mapped, GeometryLookup
export Band
export missingval, boolmask, missingmask, replace_missing, replace_missing!,
       aggregate, aggregate!, disaggregate, disaggregate!, mask, mask!,
       resample, warp, zonal, crop, extend, trim, slice, combine, points,
       classify, classify!, mosaic, mosaic!, extract, rasterize, rasterize!,
       coverage, coverage!, setcrs, setmappedcrs, smapseries, cellsize, cellarea
export crs, mappedcrs, mappedindex, mappedbounds, projectedindex, projectedbounds
export reproject, convertlookup
export Extent, extent

const DD = DimensionalData
const DA = DiskArrays
const GI = GeoInterface
const GO = GeometryOps
const LA = Lookups

# DimensionalData documentation urls
const DDdocs = "https://rafaqz.github.io/DimensionalData.jl/stable/api"
const DDdimdocs = joinpath(DDdocs, "#DimensionalData.Dimension")
const DDarraydocs = joinpath(DDdocs, "#DimensionalData.AbstractDimensionalArray")
const DDabssampleddocs = joinpath(DDdocs, "#DimensionalData.AbstractSampled")
const DDsampleddocs = joinpath(DDdocs, "#DimensionalData.Sampled")
const DDlocusdocs = joinpath(DDdocs, "#DimensionalData.Locus")
const DDselectordocs = joinpath(DDdocs, "#DimensionalData.Selector")
const DDtidocs = joinpath(DDdocs, "#DimensionalData.Ti")

const EXPERIMENTAL = """
    WARNING: This feature is experimental. It may change in future versions, and may
    not be 100% reliable in all cases. Please file github issues if problems occur.
    """

const DEFAULT_POINT_ORDER = (X(), Y())
const DEFAULT_TABLE_DIM_KEYS = (:X, :Y)

include("nokw.jl")
include("methods/shared_docstrings.jl")
include("lookup.jl")
include("dimensions.jl")
include("sources/sources.jl")
include("filearray.jl")
include("filestack.jl")
include("openstack.jl")
include("array.jl")
include("stack.jl")
include("series.jl")
include("crs.jl")
include("geometry_lookup.jl")

const RasterStackOrArray = Union{AbstractRasterStack,AbstractRaster}
const RasterSeriesOrStack = Union{AbstractRasterSeries,AbstractRasterStack}

include("utils.jl")
include("skipmissing.jl")

include("table_ops.jl")
include("create.jl")
include("read.jl")
include("write.jl")
include("show.jl")
include("plotrecipes.jl")

include("methods/burning/edges.jl")
include("methods/burning/allocs.jl")
include("methods/burning/array_init.jl")
include("methods/burning/geometry.jl")
include("methods/burning/point.jl")
include("methods/burning/line.jl")
include("methods/burning/polygon.jl")
include("methods/burning/extents.jl")
include("methods/burning/utils.jl")

include("methods/mask.jl")
include("methods/rasterize.jl")
include("methods/aggregate.jl")
include("methods/classify.jl")
include("methods/crop_extend.jl")
include("methods/coverage.jl")
include("methods/extract.jl")
include("methods/mosaic.jl")
include("methods/points.jl")
include("methods/replace_missing.jl")
include("methods/reproject.jl")
include("methods/slice_combine.jl")
include("methods/trim.jl")
include("methods/zonal.jl")

include("sources/grd.jl")
include("sources/commondatamodel.jl")
include("extensions.jl")

end
