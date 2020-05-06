#mhAcceptReject() computes the acceptance probability (Dashti's a)
# inputs:
#  acceptProb  acceptance probability
#  verbose   print information about decision-making if > 0
#
# outputs: true if accepted, false if rejected
#
function mhAcceptReject(acceptProb;verbose=0)
  #mhverbose && @printf "pot(can) = %.4f;  pot(cur) = %.4f;  a = %.4f; " canPot curPot  a

  compare = rand(); #get a random value between 0 & 1
  if (compare <= acceptProb)
    (verbose>1) && @printf("%.4f <= a (%.4f). candidate accepted.\n", compare, acceptProb);
    return true;
  else
    (verbose>1) && @printf("%.4f >  a (%.4f). candidate rejected.\n", compare, acceptProb);
    return false;
  end
end
