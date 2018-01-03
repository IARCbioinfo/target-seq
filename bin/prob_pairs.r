# probability to observe p pairs from a draw of k samples in a set of N pairs
prob_pair=function(p,N,k) {
  if (p==0) {
    i=0:(k-1);prod(1-i/(2*N-i))
  } else {
    i1=0:(p-1)
    i2=0:(k-p-1)
    prod(choose(k-2*i1,2))*prod(2*N-2*i2)/factorial(p)/choose(2*N,k)/factorial(k)
  }
}

# probability to observe p or more pairs from a draw of k samples in a set of N pairs
pvalue_pair=function(p,N,k) {
  sum(unlist(lapply(p:(floor(k/2)+1),prob_pair,N,k)))
}
