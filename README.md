# DynamicDiscreteModels

Back-end support for doing statistical inference with "partially observed Markov chain": z_t=(x_t,y_t) is Markov, y_t is observed but x_t is not observed. Hidden Markov Models are a popular example of such models.


## Installation

~~~julia
julia> Pkg.clone("git://github.com/BenConnault/ParametricModels.jl.git")
julia> Pkg.clone("git://github.com/BenConnault/DynamicDiscreteModels.jl.git")
~~~

## Goals and Usage 

The objectives of the packages include:

- filter/smoother (aka. forward/backward) algorithms
- parameter estimation: maximum likelihood and EM algorithm (known as the Baum-Welch algorithm in this context)
- hidden state inference: Viterbi filtering
- support for deep modelling of the transition matrices (the transition matrices can be arbitrary functions of a deeper statistical parameter)
- support for sparse transition matrices
- support for Jacobians in all of the above

The `DynamicDiscreteModels` module mostly defines an abstract type `DynamicDiscreteModel` which inherits from `ParametricModel`.  Any concrete instance of a `DynamicDiscreteModel` must implement:

- a field `m` which is the transition matrix for (x,y) indexed as m[x,y,x',y']
- a field `mu` which is the initial distribution of (x_1,y_1) indexed as mu[x_1,y_1]
- several other fields used for efficient implementations of filtering/smoothing.




