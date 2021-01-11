# InfDimMCMC.jl
Julia package for MCMC algorithms well-defined for infinite-dimensional unknowns

To run, 

1. Setup the MCMC problem via the `mcmcProb` structure:
  - The number of samples: `nsamp`
  - The number of burnin samples: `nburn`
  - The prior measure (from the Distributions package): `prior`
  - The MCMC method: `step` (typically set with `mcmcSetSampler()`)

2. Define the following functions:
  - Sample to model parameter map: `InfDimMCMC.mcmcSampToParamMap`
  - Parameter to solution map: `InfDimMCMC.mcmcForwardMap`
  - Solution to observation map: `InfDimMCMC.mcmcObsMap`
  - Observation to potential map: `InfDimMCMC.mcmcPotMap`
  - Gradient of the sample to model parameter map: `InfDimMCMC.mcmcGradSampToParamMap`
  - Gradient of the parameter to solution map: `InfDimMCMC.mcmcGradForwardMap`
  - Gradient of the solution to observation map: `InfDimMCMC.mcmcGradObsMap`
  - Gradient of the observation to potential map: `InfDimMCMC.mcmcGradPotMap`

  For many problems one or more of these maps may be unnecessary, so the observation to potential map is the only one that must be defined; the sample to parameter map, the parameter to solution map, and the solution to observation map are each assumed to be the identity unless otherwise specified. 

3. Run the sampler with `mcmcRun()`

A simple example is therefore (see `test/normal2D.jl`):

```
using InfDimMCMC, Distributions

m = mcmcProb();
m.nsamp=10;                        #samples
m.nburn=0;                         #burnin
meanPrior =  ones(2);
m.prior = MvNormal(meanPrior,1.0); #prior

#define potential (negative log likelihood) and its gradient
meanLlh   = -ones(2);
InfDimMCMC.mcmcPotMap(s) = -logpdf( MvNormal(meanLlh,1.0), s.samp );
InfDimMCMC.mcmcGradPotMap(s) = (s.samp-meanLlh);

#Set the sampler
smplr = "hmc|0.5|2"; #HMC sampler with step size 0.5 and 2 internal steps
mcmcSetSampler(m,smplr);

#Draw and initial sample from the prior and run
s0 = rand(m.prior);
samples, obs, lpdfs, ar = mcmcRun(m,s0);
```
