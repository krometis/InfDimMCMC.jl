#Methods for the Independence Sampler MCMC method 
#See Dashti-Stuart, 2017, Algorithm 2


#Candidate generator for Independence Sampler MCMC method
#Sample from prior
function mcmcCandidateIS(m::mcmcProb)
  #return mu0tovk(mu0samp());
  return rand(m.prior);
end

#Acceptance probability for Independence Sampler MCMC method
function mcmcAcceptProbIS(cur::mcmcSample, can::mcmcSample)
  return min( 1, exp(cur.pot - can.pot) );
end
function mcmcAcceptProbIS(curPot::Float64, canPot::Float64)
  return min( 1, exp(curPot - canPot) );
end

#run a step of the sampler
function stepIS(cur::mcmcSample, m::mcmcProb; verbose=0)
  can         = mcmcSample();
  can.samp    = mcmcCandidateIS(m);
  mcmcFillSample( can, m; computeGradients=false);
  #can.sol     = m.forwardMap(can.samp);
  #can.obs     = m.obsMap(can.sol);
  #can.pot     = m.potMap(can.obs);
  acceptProb  = mcmcAcceptProbIS(cur,can);

  accept = mhAcceptReject(acceptProb,verbose=verbose);

  #return new sample, accept(t/f), candidate
  if (accept)
    return can, true, can;
  else
    return cur, false, can; 
  end
end

