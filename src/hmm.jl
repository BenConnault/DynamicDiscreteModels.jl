type HiddenMarkovModel <: DynamicDiscreteModel
	#HMM specific field

	#DynamicDiscreteModel fields
	m::Array{Float64,4}			  	#the transition matrix given as m[x,y,x',y'] 
	mu::Array{Float64,2}  			#initial distribution (dx,dy)
	rho::Array{Float64,1}
	phi1::Array{Float64,1}	
	phi2::Array{Float64,1}
end


function hmm(mu::Array{Float64,2},a::Array{Float64,2},b::Array{Float64,2})
	dx,dy=size(mu)
	model=HiddenMarkovModel(Array(Float64,dx,dy,dx,dy),mu,Array(Float64,1),Array(Float64,dx),Array(Float64,dx))
	hmm2ddm!(model,a,b)
	model
end

#ab=(a,b)
calibrate(model::HMM,ab)=hmm2ddm!(model,ab[1],ab[2])

# m[x,y,x',y']=a[x,x']* b[x',y']
function hmm2ddm!(model::DynamicDiscreteModel,a,b)
	dx,dy=size(b)
	for jy=1:dy
		for jx=1:dx
			for iy=1:dy
				for ix=1:dx
					model.m[ix,iy,jx,jy]=a[ix,jx]*b[jx,jy]
				end
			end
		end
	end
end

#multiple dispatch will detect jacobian or no jacobian
function hmm2ddm!(model::DynamicDiscreteModel,a,b,ajac,bjac)
	dx,dy,dtheta=size(bjac)
	for jy=1:dy
		for jx=1:dx
			for iy=1:dy
				for ix=1:dx
					model.m[ix,iy,jx,jy]=a[ix,jx]*b[jx,jy]
					for itheta=1:dtheta
						model.mjac[ix,iy,jx,jy,itheta]=ajac[ix,jx,itheta]*b[jx,jy]+a[ix,jx]*bjac[jx,jy,itheta]
					end
				end
			end
		end
	end
end



