# Run 1

| Parameter                | Value |
| ------------------------ | ----- |
| time resolution (ns)     | 0.01  |
| maximum time (ns)        | 20    |
| number of atoms          | 10    |
| Doppler broadening (GHz) | 10    |
| Number of lines          | 2     |                    |       |

## Classical Plots
### Instances: 1
![](202105013_tres=0.01_tmax=10.0_len-%CE%BDM=2_ntot=100000_nbar=10_bigN=10_reset=1.0_temp=5000_len-mag=2_seed=-1_classical-single.svg)

### Instances: $10^4$


![](202105013_tres=0.01_tmax=10.0_len-%CE%BDM=2_ntot=100000_nbar=10_bigN=10_reset=1.0_temp=5000_len-mag=2_seed=-1_classical-sum.svg)


## Photon Counting Plots

### Approach \#1
This approach treats the photon counting case similarly to the classical field case. The plot below was acquired by the following:
- Taking the intensity from one instance of frequencies and generating an average of 10 Poisson-distributed photons per instance in two beams
- Calculating the correlation between each beam 
$$
g_i = \sum_{j=0}^{N_t/2}\gamma_j\gamma_{i+j}	
$$

- where above, I define the following
	- $\gamma_j$ is the number of counts in bin $j$
	- $N_t$ is the number of time bins
- If $g_i$ is non-zero, I register one correlation count at a time of $t_{res}i$
- Repeat the above for each instance until the desired number of photons is reached
- Make a histogram of all time correlations
- Take the Fourier transform of the resulting histogram
	
![](202105013_tres=0.01_tmax=10.0_len-%CE%BDM=2_ntot=100000_nbar=10_bigN=10_reset=1.0_temp=5000_len-mag=2_seed=-1_photon-correlation.svg)

We can see that there is no visible correlation at the line-difference frequency. This is probably because (a) I use a boolean correlation function which wipes out the magnitude of each correlation instance, (b) the Poisson variance is large since the number of counts is small, and (c) the re-instancing the frequencies and phases does not allow the Poisson-distributed photon counts to approach an average.

### Approach \#2


![](202105013_tres=0.01_tmax=10.0_len-%CE%BDM=2_ntot=100000_nbar=10_bigN=10_reset=1.0_temp=5000_len-mag=2_seed=-1_coincident-counts-vs-time.svg)