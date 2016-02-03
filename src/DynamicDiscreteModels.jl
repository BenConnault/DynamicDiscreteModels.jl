module DynamicDiscreteModels

include("statisticalmodels-stopgap.jl")
include("statisticalmodels-sugar.jl")

#need to import StatsBase.StatisticalModel once I make the transition
import Optim
import Distributions: Dirichlet, wsample 



export coef!, rand, loglikelihood, mle, dim,
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
