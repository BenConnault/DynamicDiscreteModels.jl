module DynamicDiscreteModels

importall ParametricModels

import Optim
import Distributions: Dirichlet, wsample 



export calibrate, simulate, loglikelihood, mle, dim,
	DynamicDiscreteModel, estep, em, viterbi, hmm2ddm!



include("utils.jl")

include("dynamicdiscretemodel.jl")
include("estep.jl")
include("emalgorithm.jl")
include("viterbi.jl")
include("loglikelihood.jl")

include("hmm.jl")


end
