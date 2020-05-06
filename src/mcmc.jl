## METROPOLIS-HASTINGS FUNCTIONS ##
#See Dashti-Stuart, Section 5


#mhAcceptReject() computes the acceptance probability (Dashti's a)
# inputs:
#  curPot      potential (neg log likelihood) for current sample
#  canPot      potential (neg log likelihood) for candidate sample
#  mhverbose   print information about decision-making (boolean, default=false)
#
function mhAcceptReject(curPot,canPot;mhverbose=false)
  a = min( 1, exp(curPot - canPot) );
  mhverbose && @printf "pot(can) = %.4f;  pot(cur) = %.4f;  a = %.4f; " canPot curPot  a

  compare = rand(); #get a random value between 0 & 1
  if (compare <= a)
    mhverbose && @printf("%.4f <= a (%.4f). candidate accepted.\n",compare,a);# sample is now %.4f\n" compare can
    return true;
  else
    mhverbose && @printf("%.4f >  a (%.4f). candidate rejected.\n",compare,a);# staying with %.4f\n" compare cur
    return false;
  end
end
#mhAcceptReject() computes the acceptance probability (Dashti's a)
# inputs:
#  acceptProb  acceptance probability
#  mhverbose   print information about decision-making (boolean, default=false)
#
# outputs: true if accepted, false if rejected
#
function mhAcceptReject(acceptProb;mhverbose=false)
  #mhverbose && @printf "pot(can) = %.4f;  pot(cur) = %.4f;  a = %.4f; " canPot curPot  a

  compare = rand(); #get a random value between 0 & 1
  if (compare <= acceptProb)
    mhverbose && @printf("%.4f <= a (%.4f). candidate accepted.\n", compare, acceptProb);
    return true;
  else
    mhverbose && @printf("%.4f >  a (%.4f). candidate rejected.\n", compare, acceptProb);
    return false;
  end
end


#mhStep() runs one mh step and returns the next sample
# inputs:
#  cur         current sample
#  y           data y
#  meth        MCMC method
#  mhverbose   print information about decision-making (boolean, default=false)
#
function mhStep(cur::mcmcSample; y=y, meth::AbstractString="pCN", 
    mhverbose=false)

  #start building new sample
  can = mcmcSample();

  if (lowercase(meth)=="ind" || lowercase(meth)=="independence")
    can.samp    = mcmcCandidateIS();
    can.sol     = adSolve(can.samp);
    can.obs     = obs(can.sol);
    can.pot     = phi(can.obs);
    acceptProb  = mcmcAcceptProbIS(cur,can);
  
  elseif (lowercase(meth)=="pcn")
    can.samp    = mcmcCandidatePCN(cur, pcn_beta);
    can.sol     = adSolve(can.samp);
    can.obs     = obs(can.sol);
    can.pot     = phi(can.obs);
    acceptProb  = mcmcAcceptProbPCN(cur,can);

  elseif (lowercase(meth)=="mala" || lowercase(meth)=="pcnl")
    can.samp    = mcmcCandidateMALA(cur, mala_h);
    can.sol     = adSolve(can.samp);
    can.obs     = obs(can.sol);
    can.pot     = phi(can.obs);
    can.gradPot = adGradPhi(can);
    acceptProb  = mcmcAcceptProbMALA(cur,can,mala_h);

  elseif (lowercase(meth)=="hmc")
    can, acceptProb   = mcmcHMC(cur, hmc_e, hmc_tau);

  else
    error("method $method not implemented in mhStep()\n");
  end
  
  #accept or reject the new sample
  accept = mhAcceptReject(acceptProb,mhverbose=mhverbose);
  if (accept)
    return can, true, can.pot;
  else
    return cur, false, can.pot; #return candidate potential so we can see how realistic the candidates have been
  end
end


