importall DynamicDiscreteModels 
#need to `importall` rather than `using` because we will extend `calibrate()` and `dim()`
#in real applications this would typically be handled in a front-end package ToyModel.jl

type ToyModel <: DynamicDiscreteModel

	#DynamicDiscreteModel fields	
	m::Array{Float64,4}			  	#the transition matrix given as m[x,y,x',y'] 
	mu::Array{Float64,2}  			#initial distribution (dx,dy)

	#DynamicDiscreteModel's technical fields	
	rho::Array{Float64,1}
	phi::Array{Float64,1}	
	psi::Array{Float64,1}
end

dx=2
dy=3

toymodel()=ToyModel(Array(Float64,dx,dy,dx,dy),fill(1/6,2,3),Array(Float64,1),Array(Float64,dx),Array(Float64,dx))

function calibrate(model::ToyModel,theta::Tuple)
	p1,p2=theta[1],theta[2]
	a=[p1 1-p1;1-p1 p1]
	p3=(1-p2)/2
	b=[p2 p3 p3; p3 p3 p2]
	hmm2ddm!(model,a,b)
end

theta2eta(theta::Tuple)=[log(theta[1]/(1-theta[1])),log(theta[2]/(1-theta[2]))]
eta2theta(eta::Array)=(exp(eta[1])/(1+exp(eta[1])),exp(eta[2])/(1+exp(eta[2])))
calibrate(model::ToyModel,eta::Array)=calibrate(model,eta2theta(eta))

dim(model::ToyModel)=2




