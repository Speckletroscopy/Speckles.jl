"""
	function poissonCount(nbar::Real)

Returns Poisson distributed counts for average count rate nbar
"""
function poissonCount(nbar::Real)
    p = exp(-nbar)
	s = p
	r = rand()
	count = 0
	while r > s
		count += 1
		p *= nbar/count
		s += p
	end
	return count
end

export poissonCount

"""
	function beCount(nbar::Real)

Returns Bose-Einstein distributed counts for average count rate nbar
"""
function beCount(nbar::Real)
	p = 1/(nbar+1)
	fnbar = p*nbar
	f = p*nbar
	s = p
	r = rand()
	count = 0
	while r>s
		count += 1
		p *= f
		s += p
	end
	return count
end

export poissonCount

"""
	function γIntensity(intensity::Vector,nbar::Real)

Calculates the photon count rate in each bin of an intensity histogram
"""
function γIntensity(intensity::Vector,nbar::Real)
	nintensity = intensity/sum(intensity)
	return nbar*nintensity
end

export γIntensity

