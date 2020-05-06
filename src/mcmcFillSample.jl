function mcmcFillSample(s::mcmcSample, m::mcmcProb; computeGradients=m.computeGradients)
  #compute the parameter associated with the sample
  s.param = mcmcSampToParamMap(s);
  #solution
  s.sol   = mcmcForwardMap(s);
  #observation
  s.obs   = mcmcObsMap(s);
  #potential
  s.pot   = mcmcPotMap(s);

  #logpdfs
  s.prLogPdf = logpdf(m.prior,s.samp);    #prior
  s.llLogPdf = -s.pot;                    #likelihood
  s.postLogPdf = s.prLogPdf + s.llLogPdf; #posterior

  #gradients
  if computeGradients
    #sample to parameter gradient
    s.gradSampToParam = mcmcGradSampToParamMap(s);
    #forward map gradient
    s.gradSol = mcmcGradForwardMap(s);
    #observation gradient
    s.gradObs = mcmcGradObsMap(s);
    #potential gradient has
    s.gradPot = mcmcGradPotMap(s);
  end
end
