# filtersmoother() builds the weights ("E step")
# EMalgorithm() does numerical optimization using the weights to build the objective function ("M-steps").
# In examples with specific structure (eg Hidden Markov Models), the M-step might be solved in closed-form. When this is the case, one should not use this generic EMalgorithm() implementation. filtersmoother() remains useful.

# EMalgorithm() currently does not use Jacobians. If profiling shows that most of the time is spent in the M-step rather than teh E-step, this might something useful to do. 


function filtersmoother(model::DynamicDiscreteModel,data::Array{Int,1},w::Array{Float64,4})
	T=length(data)
	dx,dy=size(model.mu)

	filter=Array(Float64,dx,T)
	filter[:,1]=model.mu[:,data[1]]/sum(model.mu[:,data[1]])
	for t=2:T
		rho=0.0
		for jx=1:dx
			filter[jx,t]=0
			for ix=1:dx
				filter[jx,t]+=model.m[ix,data[t-1],jx,data[t]]*filter[ix,t-1]
			end
			rho+=filter[jx,t]
		end
		for jx=1:dx
			filter[jx,t]/=rho
		end
	end

	smoother=Array(Float64,dx,T)
	smoother[:,T]=1
	for t=T-1:-1:1
		rho=0.0
		for ix=1:dx
			smoother[ix,t]=0
			for jx=1:dx
				smoother[ix,t]+=model.m[ix,data[t],jx,data[t+1]]*smoother[jx,t+1]
			end
			rho+=smoother[ix,t]
		end
		for ix=1:dx
			smoother[ix,t]/=rho
		end
	end

	conditional=Array(Float64,dx,dx)
	for t=1:T-1
		tempsum=0.0
		for jx=1:dx
			for ix=1:dx
				conditional[ix,jx]=model.m[ix,data[t],jx,data[t+1]]
				conditional[ix,jx]*=filter[ix,t]
				conditional[ix,jx]*=smoother[jx,t+1]
				tempsum+=conditional[ix,jx]
			end
		end
		for jx=1:dx
			for ix=1:dx
				#note that there is no risk of division by zero 
				#reason: this is called only on _observed_ data which must thus have non-zero probability
				w[ix,data[t],jx,data[t+1]]+=conditional[ix,jx]/tempsum
			end
		end
	end

end

function filtersmoother(model::DynamicDiscreteModel,data::Array{Array,1},w::Array{Float64,4})
	for i=1:length(data)
		filtersmoother(model,data[i],w)
	end
end

function EMalgorithm(model::DynamicDiscreteModel,data,thetai,L=1000)
	tol=0.01
	dx,dy=size(model.mu)
	w=zeros(dx,dy,dx,dy)
	llks=zeros(L)
	theta=copy(thetai)
	calibrate(model,theta)
	
	nonzeroindicies=find(model.m.!=0.0)

	llks[1]=loglikelihood(model,data)
	l=1
	go=true
	while go
		l+=1
		w[:]=0
		filtersmoother(model,data,w)
		function ff(xtheta)
			calibrate(model,xtheta)
			#w[] has zero at lest wherever model.m has zeros, and maybe more
			#we still must be careful of not calling log(0)
			#this is what nonzeroindicies is for
			res=0.0
			for i in nonzeroindicies
				res+=w[i]*log(model.m[i])
			end
			-res
		end
		ret = Optim.optimize(ff, theta,method=:cg,iterations=L)
		theta[:]=ret.minimum
		calibrate(model,theta)
		llks[l]=loglikelihood(model,data)
		go=(abs(llks[l]-llks[l-1])>tol && l<L)
	end
	println()
	println(" estimating dynamic discrete model via EM algorithm...")
	println("  $l iterations, final log-likelihood: $(round(llks[l],2))")
	# plot(llks[1:l])
	theta
end	