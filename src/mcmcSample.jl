#Define data type for MCMC Sample
mutable struct mcmcSample
  samp             #sample
  param            #parameter associated with this sample
  sol              #forward map applied to this sample
  obs              #observation associated with this sample (G)
  pot::Float64     #potential (Phi, negative log-likelihood)

  gradSampToParam  #gradient of sample to parameter map
  gradSol          #gradient of forward map
  gradObs          #gradient of observation map
  gradPot          #gradient of the potential

  prLogPdf         #logpdf of the prior
  llLogPdf         #logpdf of the likelihood
  postLogPdf         #logpdf of the posterior

  mcmcSample() = new()
  mcmcSample(samp,sol,obs,pot,gradPot) = new(samp,sol,obs,pot,gradPot)
  function mcmcSample(samp)
    s         = mcmcSample();
    s.samp    = samp;
    #s.sol     = adSolve(s.samp);
    #s.obs     = obs(s.sol);
    #s.pot     = phi(s.obs);
    #s.gradPot = adGradPhi(s);

    return s;
  end
end
import Base.copy
function copy(s1::mcmcSample)
  s2 = mcmcSample();
  s2.samp    = copy(s1.samp);
  s2.sol     = copy(s1.sol);
  s2.obs     = copy(s1.obs);
  s2.pot     = copy(s1.pot);
  s2.gradPot = copy(s1.gradPot);
  return s2;
end


