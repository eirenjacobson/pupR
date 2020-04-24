#include <TMB.hpp>
#include <iostream>


//Logit Transform
template <class Type>
Type logitt(Type x){
  return log(x/(Type(1)-x));
}


template <class vector<Type>>
vector<Type> conv(vector<Type>& x, vector<Type>& h)
{

  int j;
  int m = x.size();
  int n = h.size();

  vector<Type> lowerlimit(2);
  vector<Type> upperlimit(2);


  vector<Type> y(n+m-1);

  for(int k=0;k<(m+n-1);k++)
  {
    y(k) = 0;
    lowerlimit(0) = 1;
    lowerlimit(1) = k+1-n;
    upperlimit(0) = k;
    upperlimit(1) = m;

    for(int j=max(lowerlimit);j<min(upperlimit);j++)
    {
      y(k) = y(k) + x(j)*h(k-j+1);
    }
  }

  return y;

}

// Simple linear regression.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_VECTOR(x);
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(logSigma);

  vector<Type> result = conv(x,h);
  REPORT(result);

  ADREPORT(exp(2*logSigma));
  Type nll = -sum(dnorm(Y, a+b*x, exp(logSigma), true));
  return nll;
}

