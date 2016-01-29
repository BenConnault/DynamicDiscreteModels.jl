module DynamicDiscreteModels

importall ParametricModels

import Optim
import Distributions: Dirichlet, wsample 



export calibrate!, simulate, loglikelihood, mle, dim,
	DynamicDiscreteModel, estep, em, viterbi, hmm2ddm!



include("dynamicdiscretemodel.jl")
include("simulate.jl")
include("estep.jl")
include("emalgorithm.jl")
include("viterbi.jl")
include("loglikelihood.jl")

#non-core convenience function for front-end packages
include("hmm.jl")


end
