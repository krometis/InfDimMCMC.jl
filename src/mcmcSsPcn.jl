#Methods for the "slingshot" preconditioned Crank-Nicolson (pCN) MCMC method 
#See:
# Cotter et al, 2013
# Dashti-Stuart, 2017, Algorithm 3


function mcmcCandidateSsPcn(sample::AbstractArray,m::mcmcProb,beta::Float64)
  ##vk = sqrt(1-beta^2)*cur.samp + beta*mu0tovk(mu0samp());
  ##For Gaussian not centered at 0, need to convert uk to a mu0 sample, pull out the mean, do the update, re-add the mean, and convert to vk
  #xij = sqrt(1-beta^2)*vktomu0(cur.samp) + beta*mu0samp() - (beta + sqrt(1-beta^2))*mu0_mean;
  #return mu0tovk(xij + mu0_mean);
  
  mu0=m.prior;
  mu0_mean=mean(mu0);
  xi = rand(mu0);
  xi_logpdf = logpdf(mu0,xi);

  prop = ( sqrt(1-beta^2)*(sample-mu0_mean) + beta*xi + (1-beta)*mu0_mean );
  return prop, -xi_logpdf;
end

#stepSsPcn(): Run a step of the sampler
#Parameters:
# cur         Structure containing current sample
# m           Structure describing MCMC problem (need the prior from this)
# beta        "step size" parameter (see Dashti-Stuart)
# nProp       Number of proposals (one more than this including the current sample)
# verbose     How much information to print about each step (optional)
# recompute   Whether to recompute all information about the selected sample at the end (optional)
#
function stepSsPcn(cur::mcmcSample, m::mcmcProb, beta::Float64, nProp::Int64; verbose=3,recompute=true)
  #preallocate proposals and potentials
  proposals = zeros(length(cur.samp),nProp+1);
  pots      = zeros(nProp+1);
  propPots  = zeros(nProp+1);

  #copy current value into last entry
  proposals[:,end] = cur.samp;
  pots[end]        = cur.pot;
  propPots[end]    = -logpdf(m.prior,zeros(length(cur.samp)));

  #create all proposals
  can = mcmcSample(); #candidate
  for p=1:nProp
    can.samp,propPots[p] = mcmcCandidateSsPcn(cur.samp,m,beta);      #sample
    mcmcFillSample( can, m ; computeGradients=false);  #compute potential
    proposals[:,p] = can.samp; #add sample to list
    pots[p]        = can.pot;  #add potential to list
    (verbose>3) && @printf("  ss-pCN: Sample %d out of %d has potential %g and proposal potential %g (weight=%g)\n", p, length(pots), pots[p], propPots[p], exp(propPots[p]-pots[p]));
  end
  (verbose>3) && @printf("  ss-pCN: Sample %d out of %d has potential %g and proposal potential %g (weight=%g)\n", length(pots), length(pots), pots[end], propPots[end], exp(propPots[end]-pots[end]));

  #acceptance probabilities 
  #(unnormalized) acceptance probabilities 
  wgts = exp.(propPots-pots);

  #pick a sample based on acceptance probabilities
  choice = sample(1:length(wgts),Weights(wgts));
  accept = ( choice != length(pots) );

  (verbose>2) && @printf("ss-pCN: Selected sample %d out of %d with potential %10.6f (potential range was %10.6f to %10.6f)\n", choice, length(pots), pots[choice], minimum(pots), maximum(pots));

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
