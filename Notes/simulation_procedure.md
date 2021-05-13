# Simulation Procedure

We wish to measure $g^{(2)}(\tau)$ using photon counting, not field measurements. Field values are translated into photon counts as follows:

1. Generate frequencies, phases, and magnitudes for a time dependent electric field $E(t)$, given by equation (1) in [(Merlin 2021)](https://arxiv.org/abs/2104.10200).
2. Generate a time array $\tilde{t} = \{t_0,t_1,...,t_K\}$ where values are evenly spaced by the timing resolution $\delta t = t_{k+1}-t_k$ up to twice the measurement period $2T = t_K$. Doubling the measurement period allows us to calculate correlations out to $T$ without any zero-padding.
3. Generate field values $E(\tilde{t}) = \tilde{E} = \{E_0,E_1,...,E_K\}$
4. Apply a beamsplitter to the field to get $\tilde{E}\underset{b.s.}{\rightarrow}\left\{\tilde{E}_A,\tilde{E}_B\right\}$ 
5. For now, $\tilde{E}_A = \tilde{E}_B = \frac{\sqrt{2}}{2}\tilde{E}$.
6. Calculate the intensity of each beam $\tilde{I}_{\text{beam}} = \tilde{E}_{\text{beam}}\cdot\tilde{E}_{\text{beam}}^*$
7. Renormalize the beam intensity to represent the average photon counts over time $T$ with $\tilde{I}_{\bar{n}} = T\frac{d\bar{n}}{dt}\frac{\tilde{I}_{\text{beam}}}{\sum_{k=1}^KI_{\text{beam},k}}$ where $\frac{d\bar{n}}{dt}$ is the average photon count rate.
8. Generate counts for each beam using the inverted Poisson CDF to get $\tilde{I}_\gamma = \{n_0,n_1,...,n_K\}$ where $n_k\in \mathbb{Z}\ge 0$. 
9. Calculate the correlation time $\tau$ between beams $A$ and $B$. 
	- When generating all photon counts in each measurement period, this is calculated by offsetting beams $A$ and $B$ by $\tau/\delta t = k$ such that $G(\tau_i) = \sum_{k=0}^{K/2}I_{A,k}I_{B,k+i}$
	- When generating only photon pairs, we use sparse arrays so we can directly calculate the correlation time by accessing the indices where photon counts are non-zero. 
10. Repeat steps 1-9 for the re-instancing case. Only repeat steps 8 and 9 if using the same field instances. Continue until the desired number of photon counts is reached.
11. Generate a histogram of all correlation times
12. Apply any time-domain noise-reduction algorithms
13. Take the Fourier transform
