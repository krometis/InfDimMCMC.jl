#Methods for the multi-proposal preconditioned Crank-Nicolson (pCN) MCMC method 
#See:
# Cotter et al, 2013
# Dashti-Stuart, 2017, Algorithm 3


#Proposal kernel for pCN MCMC method
function mcmcCandidatePcn(cur::mcmcSample,m::mcmcProb,beta::Float64)
  return mcmcCandidatePcn(cur.samp,m,beta);
end

function mcmcCandidatePcn(sample::AbstractArray,m::mcmcProb,beta::Float64)
  ##vk = sqrt(1-beta^2)*cur.samp + beta*mu0tovk(mu0samp());
  ##For Gaussian not centered at 0, need to convert uk to a mu0 sample, pull out the mean, do the update, re-add the mean, and convert to vk
  #xij = sqrt(1-beta^2)*vktomu0(cur.samp) + beta*mu0samp() - (beta + sqrt(1-beta^2))*mu0_mean;
  #return mu0tovk(xij + mu0_mean);
  
  mu0=m.prior;
  mu0_mean=mean(mu0);
  xi = rand(mu0);

  return ( sqrt(1-beta^2)*(samp-mu0_mean) + beta*xi + (1-beta)*mu0_mean );
end

#run a step of the sampler
function stepMcPcn(cur::mcmcSample, m::mcmcProb, beta::Float64, nProp::Int64; verbose=0,recompute=true)
  #preallocate proposals and potentials
  proposals = zeros(length(cur.samp),nProp+1);
  pots      = zeros(nProp);

  #copy current value into last entry
  proposals[:,end] = cur.samp;
  pots[end]        = cur.pot;

  can = mcmcSample(); #candidate
  
  #midpoint of proposals
  mdpt = mcmcCandidatePcn(cur,m,beta);

  for p=1:nProp
    can.samp    = mcmcCandidatePcn(mdpt,m,beta);
    mcmcFillSample( can, m ; computeGradients=false);
    proposals[:,p] = can.samp;
    pots[p]        = can.pot;
  end

  #acceptance probabilities 
  #sumPots = sum(exp.(-pots));
  #acceptProb = exp.(-pots) ./ sumPots;
  #(unnormalized) acceptance probabilities 
  weights = exp.(-pots);

  #pick a sample based on acceptance probabilities
  choice = sample(1:length(weights),Weights(weights));

  can.samp = samples[choice];
  if recompute
    #recompute everything - inefficient if expensive
    mcmcFillSample( can, m ; computeGradients=false);
  else
    #assign what we saved
    can.pot = pots[choice];
  end

  #return new sample, accepted index, and potentials
  return can, choice, pots;
end
