include("../examples/toymodel.jl")

theta0=(.65, .5)
model=toymodel()
calibrate(model,theta0)
data=simulate(model,100,100)
thetahat=eta2theta(mle(model,data,etai))
tetahat2=eta2theta(em(model,data,etai))

@test norm(collect(thetahat)-collect(thetahat2))<1e-3

#make tests deterministics (using a seed)
#add tests for viterbi()
#add tests for jacobian methods