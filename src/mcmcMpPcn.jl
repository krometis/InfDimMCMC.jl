#Methods for the multi-proposal preconditioned Crank-Nicolson (pCN) MCMC method 
#See:
# Cotter et al, 2013
# Dashti-Stuart, 2017, Algorithm 3


#run a step of the sampler
function stepMpPcn(cur::mcmcSample, m::mcmcProb, beta::Float64, nProp::Int64; verbose=3,recompute=true)
  #preallocate proposals and potentials
  proposals = zeros(length(cur.samp),nProp+1);
  pots      = zeros(nProp+1);

  #copy current value into last entry
  proposals[:,end] = cur.samp;
  pots[end]        = cur.pot;

  can = mcmcSample(); #candidate
  
  #midpoint of proposals
  mdpt = mcmcCandidatePcn(cur,m,beta);

  for p=1:nProp
    can.samp = mcmcCandidatePcn(mdpt,m,beta);
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

  (verbose>2) && @printf("mp-pCN: Selected sample %d out of %d with potential %10.6f (potential range was %10.6f to %10.6f\n", choice, length(pots), pots[choice], minimum(pots), maximum(pots));

  can.samp = proposals[:,choice];
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
