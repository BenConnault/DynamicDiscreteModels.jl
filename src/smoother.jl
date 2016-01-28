function filter!(model::DynamicDiscreteModel,data::Array{Int,1})
	rho=sum(model.mu[:,data[1]])
	model.philambda[1]=log(rho)				#log-normalization factor for numerical stability
	model.phi[:,1]=model.mu[:,data[1]]/rho		#actual filter
	dx=length(pi)
	for t=2:length(data)
		mm=reshape(model.m[:,data[t-1],:,data[t]],(dx,dx))'
		model.phi[:,t]=mm*model.phi[:,t-1]
		rho=sum(model.phi[:,t])
		model.philambda[t]=model.philambda[t-1]+log(rho)
		model.phi[:,t]=model.phi[:,t]/rho
	end
end

function smoother!(model::DynamicDiscreteModel,data::Array{Int,1})
	T=length(data)
	model.sigma[:,T]=vec(sum(model.m[:,data[T-1],:,data[T]],3))
	rho=sum(model.sigma[:,T])
	model.sigma[:,T]=model.sigma[:,T]/rho
	model.sigmalambda[T]=log(rho)				#log-normalization factor for numerical stability
	dx=length(pi)
	for t=T-1:-1:1
		mm=reshape(model.m[:,data[t],:,data[t+1]],(dx,dx))
		model.sigma[:,t]=mm*model.sigma[:,t+1]
		rho=sum(model.sigma[:,t])
		model.sigmalambda[t]=model.sigmalambda[t+1]+log(rho)
		model.sigma[:,t]=model.sigma[:,t]/rho
	end
end



