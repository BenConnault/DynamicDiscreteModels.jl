import StatsBase.StatisticalModel
import Base.rand
using Optim

export coef!,rand,loglikelihood,dim,mle

coef!(model::StatisticalModel,parameter)=error("Method coef!(model::$(typeof(model)),parameter) must be implemented for $(typeof(model)).")
rand(model::StatisticalModel,T::Int)=error("Method rand(model::$(typeof(model)),T::Int) must be implemented for $(typeof(model)).")
loglikelihood(model::StatisticalModel,data)=error("Method loglikelihood(model::$(typeof(model)),data,parameter) must be implemented for $(typeof(model)).")


#return the dimension of the parameter: used to provide a random starting value when calling mle()
dim(model::StatisticalModel)=error("If you want to use the default MLE method on a $(typeof(model)), you must define dim(model::$(typeof(model))).")


function mle_nojac(model::StatisticalModel,data,thetai::Array{Float64,1}=rand(dim(model)),L=15000)
	ff(theta)=-loglikelihood(model,data,theta)
	df=Optim.DifferentiableFunction(ff)	
	ret = Optim.optimize(df, thetai,method=:cg,iterations=L)
	println("  no jacobian, $(ret.iterations) iterations, $(ret.f_calls) evaluations, final log-likelihood: $(round(-ret.f_minimum,4))")
	ret.minimum
end

function mle_jac(model::StatisticalModel,data,thetai::Array{Float64,1}=rand(dim(model)),L=15000)
	ff(theta)=-loglikelihood(model,data,theta)
	function fj!(theta,jac)
		jac[:]=-loglikelihood_jac(model,data,theta)[2]
	end
	function ffj!(theta,jac)
		res=-loglikelihood_jac(model,data,theta)
		jac[:]=res[2]
		res[1]
	end
	df=Optim.DifferentiableFunction(ff,fj!,ffj!)
	ret = Optim.optimize(df, thetai,method=:cg,iterations=L)
	println("  with jacobian, $(ret.iterations) iterations, $(ret.f_calls) evaluations, final log-likelihood: $(round(-ret.f_minimum,4))")
	ret.minimum
end

function mle(model::StatisticalModel,data,thetai::Array{Float64,1}=rand(dim(model)),L=15000)
	println()
	println(" generic optimization of the likelihood...")
	try
		mle_jac(model,data,thetai,L)
	catch
		mle_nojac(model,data,thetai,L)
	end
end

function fit!(model::StatisticalModel,data)
	thetahat=mle(model,data)
	coef!(model,thetahat)
end
