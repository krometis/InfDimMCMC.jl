#Methods for the "bubble bath" preconditioned Crank-Nicolson (pCN) MCMC method 
#See:
# Cotter et al, 2013
# Dashti-Stuart, 2017, Algorithm 3


#stepBbPcn(): Run a step of the sampler
#Parameters:
# cur         Structure containing current sample
# m           Structure describing MCMC problem (need the prior from this)
# beta        "step size" parameter (see Dashti-Stuart)
# nProp       Number of proposals (one more than this including the current sample)
# verbose     How much information to print about each step (optional)
# recompute   Whether to recompute all information about the selected sample at the end (optional)
#
function stepBbPcn(cur::mcmcSample, m::mcmcProb, beta::Float64, nProp::Int64; verbose=3,recompute=true)
  #preallocate proposals and potentials
  proposals = zeros(length(cur.samp),nProp+1);
  pots      = zeros(nProp+1);

  #copy current value into last entry
  proposals[:,end] = cur.samp;
  pots[end]        = cur.pot;

  # #midpoint of proposals
  # mdpt = mcmcCandidatePcn(cur,m,beta);

  #create all proposals
  can = mcmcSample(); #candidate
  for p=1:nProp
    #can.samp = mcmcCandidatePcn(mdpt,m,beta);          #sample
    can.samp = mcmcCandidatePcn(cur.samp,m,beta);      #sample
    mcmcFillSample( can, m ; computeGradients=false);  #compute potential
    proposals[:,p] = can.samp; #add sample to list
    pots[p]        = can.pot;  #add potential to list
  end

  #acceptance probabilities 
  #sumPots = sum(exp.(-pots));
  #acceptProb = exp.(-pots) ./ sumPots;
  #(unnormalized) acceptance probabilities 
  wgts = exp.(-pots);

  #pick a sample based on acceptance probabilities
  # #quick test of the sample() algorithm:
  # using StatsBase
  # pots = 100.0 * rand(10); w = exp.(-pots); probs = w ./ sum(w);
  # n=Int(1e6); choices=zeros(n); for i=1:n; choices[i] = sample(1:length(w),Weights(w)); end
  # aProbs = [ sum(choices .== i)/n for i=1:length(pots) ];
  # display([probs aProbs abs.(probs-aProbs)])
  choice = sample(1:length(wgts),Weights(wgts));
  accept = ( choice != length(pots) );

  (verbose>2) && @printf("bb-pCN: Selected sample %d out of %d with potential %10.6f (potential range was %10.6f to %10.6f)\n", choice, length(pots), pots[choice], minimum(pots), maximum(pots));

  can.samp = proposals[:,choice];
  if recompute
    #recompute everything - inefficient if expensive
    mcmcFillSample( can, m ; computeGradients=false);
  else
    #assign what we saved
    can.pot = pots[choice];
  end

  #return new sample, accept, and potentials
  return can, accept, pots;
end
