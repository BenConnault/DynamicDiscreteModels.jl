function loglikelihood(model::StatisticalModel,data,parameter)
	coef!(model,parameter)
	loglikelihood(model,data)
end

#data: several iid individuals
function loglikelihood(model::StatisticalModel,data::Array{Array,1})
	n=length(data)
	llk=0.0
	for i=1:n
		llk+=loglikelihood(model,data[i])
	end
	llk/n
end



function loglikelihood_jac(model::StatisticalModel,data,parameter)
	coef_jac!(model,parameter)
	loglikelihood_jac(model,data)
end

#data: several iid individuals
function loglikelihood_jac(model::StatisticalModel,data::Array{Array,1})
	n=length(data)
	llk,llkjac=loglikelihood_jac(model,data[1])
	for i=2:n
		temp=loglikelihood_jac(model,data[i])
		llk+=temp[1]
		llkjac+=temp[2]
	end
	llk,llkjac
end



#simulate iid individuals with heterogeneous number of periods of observation
function rand(model::StatisticalModel,T::Array{Int,1})
	n=length(T) 
	data=Array(Array,n)
	for i=1:n
		data[i]=rand(model,T[i])
	end
	data
end

# (T,n) -> n iid individuals with T periods of observation each
rand(model::StatisticalModel,T,n)=rand(model,fill(T,n))
# rand(model::StatisticalModel,Tn::Tuple{Int,Int})=rand(model,Tn[1],Tn[2])

