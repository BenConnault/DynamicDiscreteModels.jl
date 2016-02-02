#sample (x,y) from a matrix of joint probabilities.
function wsample2(mu)
	dx,dy=size(mu)
	ind2sub((dx,dy),wsample(1:dx*dy,vec(mu)))
end



function rand(model::DynamicDiscreteModel,T::Int)
	#throw error if not calibrated
	data=Array(Int,T)
	dx,dy=size(model.mu)
	x,data[1]=wsample2(model.mu)
	for t=2:T
		x,data[t]=wsample2(reshape(model.m[x,data[t-1],:,:],(dx,dy)))
	end
	data
end


