"""
	function beCount(nbar::Real)

Returns Bose-Einstein distributed counts for average count rate nbar
"""
function beCount(nbar::Real)
	p = 1/(nbar+1)
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
