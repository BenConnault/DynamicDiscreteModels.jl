#sample (x,y) from a matrix of joint probabilities.
function wsample2(mu)
	dx,dy=size(mu)
	ind2sub((dx,dy),wsample(1:dx*dy,vec(mu)))
end

# draw random stochastic matrix
rsm(dx::Int,dy::Int)=mapslices(x->rand(Dirichlet(x)),ones(dx,dy),2)
rsm(k::Int)=rsm(k,k)

# draw a sparse random stochastic matrix
function rssm(dim::Int,density=.2)
	m=Array(Float64,(dim,dim))
	for i=1:dim
		m[i,:]=Base.sprand(1,dim,density)
		m[i,rand(1:dim)]+=rand()
	end
	sparse(m./sum(m,2))	
end

#stochastic matrix normalization, ie each row is normalized to sum to 1
function nsm!(m::Array{Float64,2})
	for i in 1:size(m)[1]
		m[i,:]/=sum(m[i,:])
	end
end

#theta_1 on the diagonal and theta_2/2 on the subdiagonals
function fatdiagonal(theta,dx::Int)
	if dx==1
		eye(1)
	elseif dx==2
		[theta[1] 1-theta[1];1-theta[1] theta[1]]
	elseif dx==3
		[theta[1] (1-theta[1])/2 (1-theta[1])/2;(1-theta[1])/2 theta[1] (1-theta[1])/2;(1-theta[1])/2 (1-theta[1])/2 theta[1]]
	else
		small=(1-theta[1]-theta[2])/(dx-3)
		q=fill(small,dx,dx)
		for ix=2:dx-1
			q[ix,ix-1]=theta[2]/2
			q[ix,ix]=theta[1]
			q[ix,ix+1]=theta[2]/2
		end
		q[1,1]=theta[1]
		q[1,2]=theta[2]/2
		q[1,3]=theta[2]/2
		q[dx,dx]=theta[1]
		q[dx,dx-1]=theta[2]/2
		q[dx,dx-2]=theta[2]/2
		q
	end
end

function fatdiagonaljac(theta,dx::Int)
	if dx==1
		fill(0,1,1,2)
	elseif dx==2
		cat(3,2*eye(2)-ones(2,2),zeros(2,2))
	elseif dx==3
		cat(3,1.5*eye(3)-.5*ones(3,3),zeros(3,3))
	else
		q=fill(-1/(dx-3),dx,dx,2)
		for ix=2:dx-1
			q[ix,ix-1,1]=0
			q[ix,ix-1,2]=1/2
			q[ix,ix,1]=1
			q[ix,ix,2]=0
			q[ix,ix+1,1]=0
			q[ix,ix+1,2]=1/2
		end
		q[1,1,1]=1
		q[1,1,2]=0
		q[1,2,1]=0
		q[1,2,2]=1/2
		q[1,3,1]=0
		q[1,3,2]=1/2
		q[dx,dx,1]=1
		q[dx,dx,2]=0
		q[dx,dx-1,1]=0
		q[dx,dx-1,2]=1/2
		q[dx,dx-2,1]=0
		q[dx,dx-2,2]=1/2
		q
	end
end




fatdiagonal(dx::Int)=fatdiagonal([.6,.3],dx)



# not obvious how to do this but I don't need it for MLE:
# function ddm2hmm(ddm)
# end



function z2q!(z::Array{Float64,1},q::Array{Float64,1})
	q[:]=exp(z)/(1+sum(exp(z)))
end

#WATCH OUT: there are inconsistencies in the above and below. Should refactor. The above is used for ToyModel.
function z2q!(z::Array{Float64,1},q::Array{Float64,2})
	dx=size(q)[1]
	for ix=1:dx
		iz=(ix-1)*(dx-1)+1
		jz=ix*(dx-1)
		q[ix,1:dx-1]=exp(z[iz:jz])/(1+sum(exp(z[iz:jz])))
		q[ix,dx]=1/(1+sum(exp(z[iz:jz])))
	end
end

function q2z!(q::Array{Float64,2},z::Array{Float64,1},dtheta::Int=0)
	dx=size(q)[1]
	for ix=1:dx
		iz=dtheta+(ix-1)*(dx-1)+1
		jz=dtheta+ix*(dx-1)
		z[iz:jz]=log(q[ix,1:dx-1]/(1-sum(q[ix,1:dx-1])))
	end
end



