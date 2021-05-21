"""
	function σTemp(ω0::Real,temp::Number)

Returns the standard deviation due to Doppler broadening of frequency ω0 at temp
"""
function σTemp(ω0::Real,temp::Number)
	kbOverMhC2 = 9.178e-14;
	return sqrt(kbOverMhC2*temp)*ω0
end

export σTemp
