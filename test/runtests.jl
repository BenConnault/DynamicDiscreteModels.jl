using Base.Test
include("../examples/toymodel.jl")


theta0=(.65, .5)
model=toymodel()
coef!(model,theta0)
data=rand(model,100,100)

thetahat=eta2theta(mle(model,data))
thetahat2=eta2theta(em(model,data))
@test norm(collect(thetahat)-collect(thetahat2))<1e-3


### JACOBIANS

using Calculus

data=rand(model,10)
eta0=theta2eta(theta0)
ff(eta)=loglikelihood(model,data,eta)
@test norm(vec(Calculus.gradient(ff,eta0))-vec(loglikelihood_jac(model,data,eta0)[2])) < 1e-5

data=rand(model,10,100)
eta0=theta2eta(theta0)
ff(eta)=loglikelihood(model,data,eta)
@test norm(vec(Calculus.gradient(ff,eta0))-vec(loglikelihood_jac(model,data,eta0)[2])) < 1e-5




#make tests deterministics (using a seed)
#add tests for jacobian methods

#add tests for viterbi()
# a=fill(.5,(2,2))
# b=[ .99 .01; .01 .99]	# very little noise

# model=hmm([.5 0; 0 .5],a,b)
# data=[2,2,2,2,2,1,1,1,1,1]
# @test viterbi(model,data) == data