#Methods for the preconditioned Crank-Nicolson (pCN) MCMC method 
#See:
# Cotter et al, 2013
# Dashti-Stuart, 2017, Algorithm 3


#Proposal kernel for pCN MCMC method
function mcmcCandidatePcn(cur::mcmcSample,m::mcmcProb,beta::Float64)
  ##vk = sqrt(1-beta^2)*cur.samp + beta*mu0tovk(mu0samp());
  ##For Gaussian not centered at 0, need to convert uk to a mu0 sample, pull out the mean, do the update, re-add the mean, and convert to vk
  #xij = sqrt(1-beta^2)*vktomu0(cur.samp) + beta*mu0samp() - (beta + sqrt(1-beta^2))*mu0_mean;
  #return mu0tovk(xij + mu0_mean);
  
  mu0=m.prior;
  mu0_mean=mean(mu0);
  xi = rand(mu0);

  return ( sqrt(1-beta^2)*(cur.samp-mu0_mean) + beta*xi + (1-beta)*mu0_mean );
end

#Acceptance probability for pCN MCMC method
function mcmcAcceptProbPcn(cur::mcmcSample, can::mcmcSample)
  return min( 1, exp(cur.pot - can.pot) );
end
function mcmcAcceptProbPCN(curPot::Float64, canPot::Float64)
  return min( 1, exp(curPot - canPot) );
end

#run a step of the sampler
function stepPcn(cur::mcmcSample, m::mcmcProb, beta::Float64; verbose=0)
  can         = mcmcSample();
  can.samp    = mcmcCandidatePcn(cur,m,beta);
  mcmcFillSample( can, m ; computeGradients=false);

  acceptProb  = mcmcAcceptProbPcn(cur,can);

  accept = mhAcceptReject(acceptProb,verbose=verbose);

  #return new sample, accept(t/f), candidate
  if (accept)
    return can, true, can;
  else
    return cur, false, can; 
  end
end

#function mcmcPcnFillSample(s::mcmcSample, m::mcmcProb)
#  s.param = m.sampToParamMap(s);
#  s.sol   = m.forwardMap(s);
#  s.obs   = m.obsMap(s);
#  s.pot   = m.potMap(s);
#end
