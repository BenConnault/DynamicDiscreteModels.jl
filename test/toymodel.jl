noptimtest=2
mu0=rand(3,5)
b0=rsm(3,5)
model=toymodel(mu0,b0)
model2=toymodel(6,5)
calibrate(model,(.6,.1))
calibrate(model2,(.6,.1))


data=simulate(model,1000)
data2=simulate(model2,(1000,50))

@test sort(unique(data)) == collect(1:5)
@test sort(unique(data2[1])) == collect(1:5)
@test loglikelihood(model,data) != loglikelihood(model,data2)

viterbipath=viterbi(model,data)


@test sort(unique(viterbipath)) == collect(1:3)




md=toymodel(2,6)
theta0=[.7]
data=simulate(md,(50,1000),theta0)
println(" ** computing mle 1/$noptimtest...")
ml=mle(md,data)

@test norm(ml-theta0) <.03


md=toymodel(4,6)
#[diagonal term,sum of antidiagonal terms]: sum must be <1
theta0=[.7,.1]
data=simulate(md,(50,1000),theta0)
println(" ** computing mle 2/$noptimtest...")
ml=mle(md,data)

@test norm(ml-theta0) <.03



