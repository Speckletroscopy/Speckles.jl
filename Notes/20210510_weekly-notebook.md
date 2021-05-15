# Week of May 10, 2021

## Notes

- Wrote down the [procedure](simulation_procedure.md) followed by the simulation to generate photon correlations

- Is persistent homology a good tool to find the strongest spectral lines in correlation measurements?

  

### Speckle Noise Analysis
The results of our paper show that the average intensity correlations are given by 
$$
\label{eq:g2}
g^{(2)}(\tau) = 1-\frac{1}{N}\frac{\sum_{m=1}^M\mathcal{E}_m^4}{\left(\sum_{m=1}^M\mathcal{E}_m^2\right)^2}+\left|\frac{\sum_{m=1}^M\mathcal{E}_m^2e^{-i\Delta_m\tau}}{\sum_{m=1}^M\mathcal{E}_m^2}\right|^2\frac{\langle S(\tau)\rangle_\omega}{N^2}
$$
where
$$
\begin{align}
\label{eq:stau}
S(\tau) &= \left|\sum_{n=1}^Ne^{-i\omega_n\tau}\right|^2\\
\label{eq:stauAvg}
\langle S(\tau)\rangle_\omega &= N+N(N-1)e^{-\sigma^2\tau^2}\\
\label{eq:stauVar}
\text{Var}[S(\tau)]_\omega &= 8N(N-1)e^{-2\sigma^2\tau^2}\left[N-1+\cosh(\sigma^2\tau^2)\right]\sinh^2(\sigma^2\tau^2/2)~.
\end{align}
$$
The quantity in equation $(\ref{eq:g2})$ is what we desire to calculate in the simulation, but the noise from equation $(\ref{eq:stauVar})$ is just the contribution from Doppler broadening. We must include additional sources of noise for a realistic simulation. 

#### Shot Noise

Each bin in the intensity time series is treated as the average photon counts for a Poisson distribution over the time interval of the bin. Hence, for average counts $\bar{n}_k$ in bin $k$, we expect fluctuations of $\sqrt{\bar{n}_k}$. 



## Test program structure

- Set up parameter dictionary
- Instancing or not?
- Set up objects that persist over all instances
  - Photon counts
  - Field dataframes
  - Correlation times
- Loop over instances

 

