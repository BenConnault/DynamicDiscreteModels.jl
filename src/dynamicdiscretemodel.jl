abstract DynamicDiscreteModel <: ParametricModel

	# m::Array{Float64,4}			  	#the transition matrix given as m[x,y,x',y'] 
	# mu::Array{Float64,2}  			#initial distribution (dx,dy)
	# mjac::Array{Float64,5}			#jacobian
	
	# #discrete filter variables
	# rho::Array{Float64,1}				#container for the constant rho used in the discrete filter
	# phi::Array{Float64,1}			#filter value today
	# psi::Array{Float64,1}			#filter value tomorrow


	#if you want to implement jacobian too
	# rhojac::Array{Float64,1}				
	# phijac::Array{Float64,2}			
	# psijac::Array{Float64,2}
		

function simulate(model::DynamicDiscreteModel,T::Int)
	#throw error if not calibrated
	data=Array(Int,T)
	dx,dy=size(model.mu)
	x,data[1]=wsample2(model.mu)
	for t=2:T
		x,data[t]=wsample2(reshape(model.m[x,data[t-1],:,:],(dx,dy)))
	end
	data
end


