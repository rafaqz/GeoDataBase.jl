# Share common docstrings here to keep things consistent

const NAME_KEYWORD = """
- `name`: a `Symbol` name for a Raster, which will also retrieve the 
    a named layer if `Raster` is used on a multi-layered file like a NetCDF. 
"""
const METADATA_KEYWORD = """
- `metadata`: `Dict` or `Metadata` object for the array, or `NoMetadata()`.
"""
const REFDIMS_KEYWORD = """
- `refdims`: `Tuple of` position `Dimension`s the array was sliced from, defaulting to `()`.
    Usually not needed.
"""

const TO_KEYWORD = """
- `to`: a `Raster`, `RasterStack`, `Tuple` of `Dimension` or `Extents.Extent`.
    If no `to` object is provided the extent will be calculated from the geometries,
    Additionally, when no `to` object or an `Extent` is passed for `to`, the `size`
    or `res` keyword must also be used.
"""
const SIZE_KEYWORD = """
- `size`: the size of the output array, as a `Tuple{Int,Int}` or single `Int` for a square.
    Only required when `to` is not used or is an `Extents.Extent`, and `res` is not used.
"""
const RES_KEYWORD = """
- `res`: the resolution of the dimensions (often in meters or degrees), a `Real` or `Tuple{<:Real,<:Real}`.
    Only required when `to` is not used or is an `Extents.Extent`, and `size` is not used.
"""
const CRS_KEYWORD = """
- `crs`: a `crs` which will be attached to the resulting raster when `to` not passed
   or is an `Extent`. Otherwise the crs from `to` is used.
"""
const BOUNDARY_KEYWORD = """
- `boundary`: for polygons, include pixels where the `:center` is inside the polygon,
    where the polygon `:touches` the pixel, or that are completely `:inside` the polygon.
    The default is `:center`.
"""
const SHAPE_KEYWORD = """
- `shape`: Force `data` to be treated as `:polygon`, `:line` or `:point` geometries.
    using points or lines as polygons may have unexpected results.
"""
const THREADED_KEYWORD = """
- `threaded`: run operations in parallel, `false` by default. In some circumstances `threaded` 
    can give large speedups over single-threaded operation. This can be true for complicated 
    geometries written into low-resolution rasters, but may not be for simple geometries with 
    high-resolution rasters. With very large rasters threading may be counter productive due 
    to excessive memory use. Caution should also be used: `threaded` should not be used in in-place 
    functions writing to `BitArray` or other arrays where race conditions can occur. 
"""

const GEOM_KEYWORDS = """
$TO_KEYWORD
$RES_KEYWORD
$SIZE_KEYWORD
$CRS_KEYWORD
$BOUNDARY_KEYWORD
$SHAPE_KEYWORD
"""

const DATA_ARGUMENT = """
- `data`: a GeoInterface.jl `AbstractGeometry`, a nested `Vector` of `AbstractGeometry`,
    or a Tables.jl compatible object containing a `:geometry` column or points and values columns,
    in which case `geometrycolumn` must be specified.
"""

const OBJ_ARGUMENT = """
a [`Raster`](@ref) or one or multiple geometries. Geometries can be
    a GeoInterface.jl `AbstractGeometry`, a nested `Vector` of `AbstractGeometry`,
    or a Tables.jl compatible object containing a `:geometry` column or points and values columns,
    in which case `geometrycolumn` must be specified."""

const GEOMETRYCOLUMN_KEYWORD = """
- `geometrycolumn`: `Symbol` to manually select the column the geometries are in
    when `data` is a Tables.jl compatible table, or a tuple of `Symbol` for columns of
    point coordinates.
"""


const LAZY_KEYWORD = """
- `lazy`: A `Bool` specifying if to load data lazily from disk. `false` by default.
"""
const FILENAME_KEYWORD = """
- `filename`: a filename to write to directly, useful for large files.
"""
const SUFFIX_KEYWORD = """
- `suffix`: a string or value to append to the filename.
    A tuple of `suffix` will be applied to stack layers. `keys(stack)` are the default.
"""
const PROGRESS_KEYWORD = """
- `progress`: show a progress bar, `true` by default, `false` to hide.
"""
const VERBOSE_KEYWORD = """
- `vebose`: whether to print messages about potential problems. `true` by default.
"""
const SOURCE_KEYWORD = """
- `source`: Usually automatically detected from filepath extension. 
    To manually force, a `Symbol` can be passed `:gdal`, `:netcdf`, `:grd`, `:grib`.
    The internal [`Rasters.Source`](@ref) objects, such as `Rasters.GDALsource()`, 
    `Rasters.GRIBsource()` or `Rasters.NCDsource()` can also be used.
"""
const EXT_KEYWORD = """
- `ext`: filename extension such as ".tiff" or ".nc". 
    Used to specify specific files if only a directory path is used.
"""
const FORCE_KEYWORD = """
- `force`: `false` by default. If `true` it force writing to a file destructively, even if it already exists.
"""
const DROPBAND_KEYWORD = """
- `dropband`: drop single band dimensions when creating stacks from filenames. `true` by default.
"""

const CONSTRUCTOR_CRS_KEYWORD = """
- `crs`: the coordinate reference system of  the objects `XDim`/`YDim` dimensions.
    Only set this if you know the detected crs is incorrect, or it is not present in
    the file. The `crs` is expected to be a GeoFormatTypes.jl `CRS` or `Mixed` mode `GeoFormat` object,
    like `EPSG(4326)`.
"""
const CONSTRUCTOR_MAPPEDCRS_KEYWORD = """
- `mappedcrs`: the mapped coordinate reference system of the objects `XDim`/`YDim` dimensions.
    for `Mapped` lookups these are the actual values of the index. For `Projected` lookups
    this can be used to index in eg. `EPSG(4326)` lat/lon values, having it converted automatically.
    Only set this if the detected `mappedcrs` in incorrect, or the file does not have a `mappedcrs`,
    e.g. a tiff. The `mappedcrs` is expected to be a GeoFormatTypes.jl `CRS` or `Mixed` mode `GeoFormat` type.
"""
const GROUP_KEYWORD = """
- `group`: the group in the dataset where `name` can be found. Only needed for nested datasets.
    A `String` or `Symbol` will select a single group. Pairs can also used to access groups
    at any nested depth, i.e `group=:group1 => :group2 => :group3`.
"""

const CHECKMEMORY_KEYWORD = """
- `checkmemory`: if `true` (the default), check if there is enough memory for the operation. 
    `false` will ignore memory needs.
"""

const SCALE_KEYWORD = """
- `scale`: set `scale` for `x * scale + offset` transformations. 
"""

const OFFSET_KEYWORD = """
- `offset`: set `offset` for `x * scale + offset` transformations. 
"""

const RAW_KEYWORD = """
- `raw`: turn of all scaling and masking and load the raw values from disk.
    `false` by default. If `true`, `scaled` will be set to `false` and `missingval`
    will to the existing missing value in the file. A warning will be printed if 
    `scaled` or `missingval` are manually set to another value.
"""

const SCALED_KEYWORD = """
- `scaled`: apply scale and offset as `x * scale + offset` where 
    `scale` and/or `offset` are found in file metadata. `true` by default.
    This is common where data has been convert to e.g. UInt8 to save disk space.
    To ignore `scale` and `offset` metadata, use `scaled=false`. 
    Note 1: If `scale` and `offset` are `1.0` and `0.0` they will be ignored and the 
    original type will be used even when `scaled=true`. This is because these values 
    may be fallback defaults and we do not want to convert every `Real` array to larger
    `Float64` values. 
    Note 2: `raw=true` will ignore `scaled` and `missingval` and return
    the raw values.
"""

const COERCE_KEYWORD = """
- `coerce`: where `scale` and/or `offset` are present during `setindex!` to disk, 
    coerce values to the element type used on dist. `convert` is the default, 
    but `round`, `trunc` or or `ceil` or other functions with `f(::Type{T}, x)`
    signature may be needed where the values are not exact.
"""

const MISSINGVAL_KEYWORD = """
- `missingval`: value representing missing data, normally detected from the file and 
    automatically converted to `missing`. Setting to an alternate value, such as `0` 
    or `NaN` may be desirable for improved perfomance. `nothing` specifies no missing value. 
    Using the same `missingval` the file already has removes the overhead of replacing it,
    this can be done by passing the `missingval` function as `missingval`. 
    If the file has an incorrect value, we can manually define the transformation
    as a pair like `correct_value => missing` or `correct_value => NaN`.
    `correct_value => correct_value` will keep remove the overhead of changing it. 
    Note: When `raw=true` is set, `missingval` is not changed from the value specified
    in the file.
"""

const WRITE_MISSINGVAL_KEYWORD = """
- `missingval`: set the missing value (i.e. FillValue / nodataval) of the written raster,
    as Julia's `missing` cannot be stored. If not passed in, an appropriate `missingval` 
    will be detected from the objects `missingval`, its `metadata`, or a default will be 
    chosen base on the array element type(s).
"""

const CHUNKS_KEYWORD = """
- `chunks`: a `NTuple{N,Int}` specifying the chunk size for each dimension. 
    To specify only specific dimensions, a Tuple of `Dimension` wrapping `Int` 
    or a `NamedTuple` of `Int` can be used. Other dimensions will have a chunk
    size of `1`. `true` can be used to mean: use the original 
    chunk size of the lazy `Raster` being written or X and Y of 256 by 256.
    `false` means don't use chunks at all.
"""