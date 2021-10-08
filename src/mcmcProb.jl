#Define data type for MCMC problem
mutable struct mcmcProb
  nsamp::Int64     #Number of (retained) samples
  nburn::Int64     #Number of burnin (discarded) samples

  computeGradients::Bool  #whether to compute gradients as part of mcmcFillSample

#  cur::mcmcSample  #Current sample
  #samp0            #Initial sample
  #data             #Data (y)

#  #The next two attributes define the map from a sample to the problem
#  #parameters, and its gradient. For many problems, this map may be the
#  #identity. However, in some cases the user may want to tweak the definition
#  #of a "sample" to be different from the (physical) problem parameters, in
#  #which case these two attributes can be used without touching, e.g.,
#  #gradPotMap.
#  sampToParamMap
#  gradSampToParamMap
#
#  #These two attributes define the map from a given parameter value to its 
#  #associated solution (e.g., of a PDE). For some problems, this may
#  #be the map from a parameter directly to observations, in which case the 
#  #observation map (see below) can be set to the identity.
#  forwardMap       #Map from parameter to solution
#  gradForwardMap   #Gradient of map from parameter to solution
#
#  #This attribute defines the map from a solution to its associated 
#  #observations.
#  obsMap           #Map from solution to observations
#  gradObsMap       #Gradient of map from solution to observations
#
#  #These two attributes define the map from observations to the potential
#  #(negative log-likelihood) - see Dashti, Stuart, 2017 for the definition.
#  #llhMap           #Log-likelihood
#  potMap           #Potential (Phi, negative log-likelihood)
#  gradPotMap       #Gradient of potential (Phi, negative log-likelihood)

  prior            #Prior measure

  #proposalKernel   #Proposal kernel
  #acceptReject     #Accept/reject method

  mcmc::Dict       #Dictionary describing the MCMC method
  step             #Function to run one MCMC step

  mcmcProb() = new()
  #mcmcSample(samp,sol,obs,pot,gradPot) = new(samp,sol,obs,pot,gradPot)
  #function mcmcSample(samp)
  #  s         = mcmcSample();
  #  s.samp    = samp;
  #  #s.sol     = adSolve(s.samp);
  #  #s.obs     = obs(s.sol);
  #  #s.pot     = phi(s.obs);
  #  #s.gradPot = adGradPhi(s);

  #  return s;
  #end
end

