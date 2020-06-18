function mcmcParse(mcmc::String; delim="|")
  mcmcArray = split(mcmc,delim);
  mcmcType = lowercase(mcmcArray[1]);

  if mcmcType == "hmc"
    if length(mcmcArray) != 3
      error("mcmc type $(mcmcType) requires two arguments (step size and # of steps)");
    end
    hmc_e   = eval(Meta.parse(mcmcArray[2]));
    hmc_I   = eval(Meta.parse(mcmcArray[3]));
    hmc_tau = hmc_e * hmc_I;
    d = Dict("mcmc" => :hmc, "eps"=>hmc_e, "tau"=>hmc_tau);

  elseif mcmcType == "lolhmc"
    if length(mcmcArray) != 3
      error("mcmc type $(mcmcType) requires two arguments (step size and # of steps)");
    end
    hmc_e   = eval(Meta.parse(mcmcArray[2]));
    hmc_I   = eval(Meta.parse(mcmcArray[3]));
    hmc_tau = hmc_e * hmc_I;
    d = Dict("mcmc" => :lolhmc, "eps"=>hmc_e, "tau"=>hmc_tau);

  elseif mcmcType == "mala"
    if length(mcmcArray) != 2
      error("mcmc type $(mcmcType) requires one arguments (step size)");
    end
    mala_h = eval(Meta.parse(mcmcArray[2]));
    d = Dict("mcmc" => :mala, "h"=>mala_h);
  
  elseif mcmcType == "pcn"
    if length(mcmcArray) != 2
      error("mcmc type $(mcmcType) requires one arguments (step size)");
    end
    pcn_beta = eval(Meta.parse(mcmcArray[2]));
    d = Dict("mcmc" => :pcn, "beta"=>pcn_beta);
  
  elseif mcmcType == "ind" || mcmcType == "is"
    if length(mcmcArray) != 1
      error("mcmc type $(mcmcType) takes no arguments");
    end
    d = Dict("mcmc" => :is);

  else
    error("mcmc type $(mcmcType) is not supported.");
  end

  return d;
end
