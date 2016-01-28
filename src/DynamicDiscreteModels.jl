module DynamicDiscreteModels

importall ParametricModels

import Optim
import Distributions: Dirichlet, wsample 



export calibrate, simulate, loglikelihood, mle, dim,
	DynamicDiscreteModel, viterbi, EMalgorithm, filtersmoother



#source files
include("utils.jl")

include("dynamicdiscretemodel.jl")
include("emalgorithm.jl")
include("viterbi.jl")
include("loglikelihood.jl")



# include("dev.jl")
# export loglikelihood2, mle2, calibrate2

end
