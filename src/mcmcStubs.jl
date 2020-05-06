#define some function stubs

#by default we assume that the following are the identity:
# The sample to parameter map
mcmcSampToParamMap     = (s -> s.samp);
mcmcGradSampToParamMap = (s -> I);
# The forward map
mcmcForwardMap         = (s -> s.param);
mcmcGradForwardMap     = (s -> I);
# The observation map
mcmcObsMap             = (s -> s.sol);
mcmcGradObsMap         = (s -> I);


#The following throw an error if called without being redefined:
# The potential map and its gradient
mcmcPotMap(s) = error("You must define the potential (negative log-likelihood) function mcmcPotMap(s::mcmcSample).");
mcmcGradPotMap(s) = error("You must define the gradient of the potential (negative log-likelihood) function mcmcGradPotMap(s::mcmcSample).");
## The sampler
#mcmcSampler(m,s;verbose=0) = error("You must set the MCMC sampler mcmcSampler(), either manually or by using mcmcSetSampler().");
