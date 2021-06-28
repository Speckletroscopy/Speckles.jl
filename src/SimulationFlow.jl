#-------------------------------------------------------------------------------

abstract type Correlation end
abstract type CorrelationVector <: Correlation end
abstract type CorrelationMatrix <: Correlation end

struct DenseCorrelationVector <: CorrelationVector
    n::Int # column of correlation matrix
    data::Vector
end

struct SparseCorrelationVector <: CorrelationVector
    n::Int # column of correlation matrix
    data::SparseVector
end

struct DenseCorrelationMatrix <: CorrelationMatrix
    data::Matrix
end

struct SparseCorrelationMatrix <: CorrelationMatrix
    data::SparseMatrixCSC
end

#-------------------------------------------------------------------------------

struct SpeckleSim{U<:SpeckleReadout,V<:Correlation}
    dt::DateTime
    id::UUID
    params::SpeckleParams
    bs::Beamsplitter
    readout::Vector{U}
    corr::Vector{V}
end

#-------------------------------------------------------------------------------
# this is what enters the simulation for each "run"
struct SpeckleInstance
    source::LightSource
    detect::Detector
    bs::Beamsplitter
    γint::γIntensity
    params::SpeckleParams
end

function SpeckleInstance(params::SpeckleParams)
    source = LightSource(
                         params.n,
                         params.Em,
                         params.νm,
                         params.σ,
                         params.fγ
                        )
    detect = Detector(
                      params.deadtime,
                      params.resolution,
                      params.jitter,
                      params.efficiency,
                      params.darkcounts
                     )
    bs = Beamsplitter(1,1)
    nb = nbar(params.duration,source)*bs.t
    γint = γIntensity(nb,params.duration,detect.resolution,source)

    return SpeckleInstance(source,  detect, bs, γint, params)
end

export SpeckleInstance
#-------------------------------------------------------------------------------
