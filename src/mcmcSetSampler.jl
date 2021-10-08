function mcmcSetSampler(mcmcP::mcmcProb,d::Dict)
  if d["mcmc"] == :hmc 
    println("Setting sampler to HMC (step size=$(d["eps"]), integration time=$(d["tau"]))");
    mcmcP.step = ( (cur,m;verbose=0) -> stepHmc(cur, m, d["eps"],d["tau"]; verbose=verbose) );
    mcmcP.mcmc = d;
    mcmcP.computeGradients = true;

  elseif d["mcmc"] == :lolhmc 
    println("Setting sampler to LOL-HMC (step size=$(d["eps"]), integration time=$(d["tau"]))");
    mcmcP.step = ( (cur,m;verbose=0) -> stepLolHmc(cur, m, d["eps"],d["tau"]; verbose=verbose) );
    mcmcP.mcmc = d;
    mcmcP.computeGradients = true;

  elseif d["mcmc"] == :mala
    println("Setting sampler to MALA (step size=$(d["h"]))");
    mcmcP.step = ( (cur,m;verbose=0) -> stepMala(cur, m, d["h"]; verbose=verbose) );
    mcmcP.mcmc = d;
    mcmcP.computeGradients = true;
    
  elseif d["mcmc"] == :pcn
    println("Setting sampler to pCN (step size=$(d["beta"]))");
    mcmcP.step = ( (cur,m;verbose=0) -> stepPcn(cur, m, d["beta"]; verbose=verbose) );
    mcmcP.mcmc = d;
    mcmcP.computeGradients = false;

  elseif d["mcmc"] == :is
    println("Setting sampler to IS");
    mcmcP.step = stepIS;
    mcmcP.mcmc = d;
    mcmcP.computeGradients = false;
  else
    error("Unrecognized/unsupported MCMC method: $(d["mcmc"])");
  end
end

function mcmcSetSampler(mcmcP::mcmcProb,s::String)
  mcmcSetSampler(mcmcP,mcmcParse(s));
end
