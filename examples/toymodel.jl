importall DynamicDiscreteModels 
#need to `importall` rather than `using` because we will extend `coef!()` and `dim()`
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

function coef!(model::ToyModel,theta::Tuple)
	p1,p2=theta[1],theta[2]
	a=[p1 1-p1;1-p1 p1]
	p3=(1-p2)/2
	b=[p2 p3 p3; p3 p3 p2]
	hmm2ddm!(model,a,b)
end

theta2eta(theta::Tuple)=[log(theta[1]/(1-theta[1])),log(theta[2]/(1-theta[2]))]
eta2theta(eta::Array)=(exp(eta[1])/(1+exp(eta[1])),exp(eta[2])/(1+exp(eta[2])))
coef!(model::ToyModel,eta::Array)=coef!(model,eta2theta(eta))

dim(model::ToyModel)=2


function lkcontour()
	theta0=(.65, .5)
	model=toymodel()
	coef!(model,theta0)
	data=rand(model,100,100)
	thetahat=eta2theta(mle(model,data))
	thetahat2=eta2theta(em(model,data))
	thetahat3=eta2theta(mle(model,data[1:10]))


	xx=linspace(.1,.9,50)
	yy=linspace(.1,.9,50)
	function ff(x,y,n)
		coef!(model,(x,y))
		loglikelihood(model,data[1:n])
	end
	gg1=[ff(x,y,10)::Float64 for x=xx,y=yy]
	gg2=[ff(x,y,100)::Float64 for x=xx,y=yy]

	qlevels=[0,.25,.5,.75,0.8,0.9,0.95,0.975,0.99,1]
	levels=quantile(vec(gg),qlevels)

	contour(xx,yy,gg,fill=true,levels=levels,color=ColorGradient(:heat,[0,0.98,1]))
	scatter!([theta0[1]],[theta0[2]],c=:green,leg=false)
	scatter!([thetahat[1]],[thetahat[2]],c=:blue,leg=false)
	contour(xx,yy,gg,fill=true,levels=levels,color=ColorGradient(:heat,[0,0.98,1]))
	scatter!([theta0[1]],[theta0[2]],c=:green,leg=false)
	scatter!([thetahat3[1]],[thetahat3[2]],c=:blue,leg=false)
end


