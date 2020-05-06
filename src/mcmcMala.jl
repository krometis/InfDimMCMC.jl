#Methods for Metropolis-Adjusted Langevin (MALA) MCMC method
#See Beskos, Girolami, et al, Algorithm 2.2

#Candidate generator for MALA MCMC method
function mcmcCandidateMala(cur::mcmcSample, m::mcmcProb, h::Float64)
  p = (4-h)/(4+h); #should be rho but p is simpler

  #C = mu0cov();
  ##some mess here to recenter to mean 0
  #canSamp = p*(vktomu0(cur.samp)-mu0_mean) + sqrt(1-p^2)*( mu0samp() - mu0_mean - 0.5*sqrt(h)*C*vktomu0(cur.gradPot) );
  #return mu0tovk( canSamp + mu0_mean );
 
  mu0 = m.prior;
  C = cov(mu0);
  mu0_mean = mean(mu0);
  xi = rand(mu0);

  #some mess here to recenter to mean 0
  canSamp = p*(cur.samp-mu0_mean) + sqrt(1-p^2)*( xi - mu0_mean - 0.5*sqrt(h)*C*cur.gradPot );
  return canSamp + mu0_mean;
end

#Acceptance probability for MALA MCMC method
function mcmcAcceptProbMala(cur::mcmcSample, can::mcmcSample, m::mcmcProb, h::Float64)
  p = (4-h)/(4+h); #should be rho but p is simpler

  # #\kappa(u,u')
  # mu0InvNormSqVk(v) = ( vktomu0(v)'*mu0cov()*vktomu0(v) )[1];
  # logKUV = -cur.pot - 0.125*h*mu0InvNormSqVk(cur.gradPot) - 0.5*sqrt(h/(1-p^2))*dot(vktomu0(cur.gradPot),vktomu0(can.samp-mu0_mean - p*(cur.samp-mu0_mean)));
  # logKVU = -can.pot - 0.125*h*mu0InvNormSqVk(can.gradPot) - 0.5*sqrt(h/(1-p^2))*dot(vktomu0(can.gradPot),vktomu0(cur.samp-mu0_mean - p*(can.samp-mu0_mean)));

  # return min( 1, exp(logKVU-logKUV) );

  mu0 = m.prior;
  C = cov(mu0);
  mu0_mean = mean(mu0);

  #\kappa(u,u')
  mu0InvNormSq(v) = dot( v,C*v );
  logKUV = -cur.pot - 0.125*h*mu0InvNormSq(cur.gradPot) - 0.5*sqrt(h/(1-p^2))*dot(cur.gradPot,can.samp-mu0_mean - p*(cur.samp-mu0_mean));
  logKVU = -can.pot - 0.125*h*mu0InvNormSq(can.gradPot) - 0.5*sqrt(h/(1-p^2))*dot(can.gradPot,cur.samp-mu0_mean - p*(can.samp-mu0_mean));

  return min( 1, exp(logKVU-logKUV) );
end


#run a step of the sampler
function stepMala(cur::mcmcSample, m::mcmcProb, h::Float64; verbose=0)
  can         = mcmcSample();
  can.samp    = mcmcCandidateMala(cur,m,h);
  mcmcFillSample(can, m; computeGradients=true);
  #can.sol     = m.forwardMap(can.samp);
  #can.obs     = m.obsMap(can.sol);
  #can.pot     = m.potMap(can.obs);
  #can.gradPot = m.gradPotMap(can);
  acceptProb  = mcmcAcceptProbMala(cur,can,m,h);

  accept = mhAcceptReject(acceptProb,verbose=verbose);

  #return new sample, accept(t/f), candidate
  if (accept)
    return can, true, can;
  else
    return cur, false, can; 
  end
end
