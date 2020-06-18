#Methods for Modified Hamiltonian MCMC (HMC) method
# see Glatt-Holtz, Krometis, Mondaini (add reference in future)
#
function stepLolHmc(cur::mcmcSample, m::mcmcProb, hmc_e::Float64, hmc_tau::Float64; verbose=2, ve=(rand(m.prior)-mean(m.prior)))
  ue   = copy(cur);            #ue is an mcmcSample

  mu0  = m.prior;
  mu0_mean = mean(mu0);
#  ve   = rand(mu0) - mu0_mean; #recenter to mean 0
  C    = cov(mu0);
  Cinv = invcov(mu0);

  mu0NormSq(v) = dot(v-mu0_mean,Cinv*(v-mu0_mean));
  H(u,v) = u.pot + 0.5*mu0NormSq(v+mu0_mean) + 0.5*mu0NormSq(u.samp);
  H0 = H(cur,ve);

  nsteps = round(Int64,hmc_tau / hmc_e);

  cose = cos(hmc_e); sine = sin(hmc_e);
  deltaH = 0.0;

  #leapfrog step
  for i=1:nsteps
  
    #save values from previous step
    u0RC = ue.samp - mu0_mean;  #RC for recentered
    u0GradPot = ue.gradPot;
    v0 = ve;

    #first half step of (12)
    vm = v0 - 0.5*hmc_e*C*u0GradPot;

    #stormer-verlet integrator of (11)
    ueRC  =  cose*u0RC + sine*vm;
    vp    = -sine*u0RC + cose*vm;
    
    #fill out sample (need gradPot for next step)
    ue.samp    = ueRC + mu0_mean;
    #mcmcFillSample(ue,m;computeGradients=true);
    #compute the parameter associated with the sample
    ue.param = mcmcSampToParamMap(ue);
    #for lolmc, don't need to compute the solution at each step
    ##solution
    #ue.sol   = mcmcForwardMap(ue);
    ##observation
    #ue.obs   = mcmcObsMap(ue);
    ##potential
    #ue.pot   = mcmcPotMap(ue);
    ##logpdfs
    #ue.prLogPdf = logpdf(m.prior,ue.samp);    #prior
    #ue.llLogPdf = -ue.pot;                    #likelihood
    #ue.postLogPdf = ue.prLogPdf + ue.llLogPdf; #posterior
    #sample to parameter gradient
    ue.gradSampToParam = mcmcGradSampToParamMap(ue);
    #forward map gradient
    ue.gradSol = mcmcGradForwardMap(ue);
    #observation gradient
    ue.gradObs = mcmcGradObsMap(ue);
    #potential gradient has
    ue.gradPot = mcmcGradPotMap(ue);

    #second half step of (12)
    ve = vp - 0.5*hmc_e*C*ue.gradPot;

    #add to deltaH (we need to know leapfrog values for this)
    deltaH += -0.5*hmc_e*( dot(v0,u0GradPot) + dot(ve,ue.gradPot) );

    Hf = H(ue,ve);
    (verbose>3) && @printf("  Step %3d: H0 = %10.6f, Hf = %10.6f, Hf-H0 = %10.6f\n", i, H0, Hf, Hf-H0);
  end

  #fill out remainder of sample (skipped above)
  #solution
  ue.sol   = mcmcForwardMap(ue);
  #observation
  ue.obs   = mcmcObsMap(ue);
  #potential
  ue.pot   = mcmcPotMap(ue);
  #logpdfs
  ue.prLogPdf = logpdf(m.prior,ue.samp);    #prior
  ue.llLogPdf = -ue.pot;                    #likelihood
  ue.postLogPdf = ue.prLogPdf + ue.llLogPdf; #posterior

  mu0InvNormSq(v) = dot(v,C*v);
  deltaH += ue.pot - cur.pot - 0.125*hmc_e^2*( mu0InvNormSq(ue.gradPot) - mu0InvNormSq(cur.gradPot) ); 
  
  Hf = H(ue,ve);

  (verbose>2) && @printf("H0 = %10.6f, Hf = %10.6f, Hf-H0 = %10.6f, deltaH = %10.6f\n", H0, Hf, Hf-H0, deltaH);

  acceptProb = min( 1, exp(-deltaH) );

  accept = mhAcceptReject(acceptProb,verbose=verbose);

  #return new sample, accept(t/f), candidate
  if (accept)
    return ue, true, ue;
  else
    return cur, false, ue; 
  end
end

