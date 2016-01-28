dx=3
dy=2
mu0=rand(dx,dy)
m0=reshape(rsm(dx*dy,dx*dy),(dx,dy,dx,dy))
model=ddm(mu0,m0)
model::DynamicDiscreteModels.DDM
model::DynamicDiscreteModels.DynamicDiscreteModel
model::DynamicDiscreteModels.ParametricModel

@test (try dim(model,1) catch; 57 end) == 57

data=simulate(model,1000)
data2=simulate(model,[1000,50])

@test sort(unique(data)) == collect(1:dy)
@test sort(unique(data2[1])) == collect(1:dy)
@test loglikelihood(model,data2,reshape(rsm(dx*dy,dx*dy),(dx,dy,dx,dy))) < loglikelihood(model,data2,m0)
@test loglikelihood(model,data2,reshape(rsm(dx*dy,dx*dy),(dx,dy,dx,dy))) < loglikelihood(model,data2,m0)
@test loglikelihood(model,data2,reshape(rsm(dx*dy,dx*dy),(dx,dy,dx,dy))) < loglikelihood(model,data2,m0)


# (2,2) HMM
# a=[.99 .01 ; .01 .99] 	# very persistent
a=fill(.5,(2,2))
b=[ .99 .01; .01 .99]	# very little noise

model=hmm([.5 0; 0 .5],a,b)
data=[2,2,2,2,2,1,1,1,1,1]
@test viterbi(model,data) == data