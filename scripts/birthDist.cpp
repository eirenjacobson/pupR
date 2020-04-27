#include <TMB.hpp>
#include <iostream>

//Bound between 0 and 1
template <class Type>
Type ilogit(Type x){
  return Type(1.0)/(Type(1.0)+exp(-x));
}

//Bound between -1 and 1
template <class Type>
Type bound(Type x){
  return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);
}

//Logit Transform
template <class Type>
Type logitt(Type x){
  return log(x/(Type(1)-x));
}

//posfun-function
template <class Type>
Type posfun(Type x, Type eps){
  if (x>=eps)
  {
    return x;
  }
  else
  {
    Type y=1.0-x/eps;
    return eps*(1./(1.+y+y*y));
  }
}

template <class Type>
vector<Type> conv(vector<Type>& x, vector<Type>& h)
{

  //int j;
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

    for(int j=CppAD::Integer(max(lowerlimit));j<=CppAD::Integer(min(upperlimit));j++)
    {
      y(k-1) = y(k-1) + x(j-1)*h(k-j);
    }
  }

  return y;

}


template<class Type>
Type objective_function<Type>::operator() ()
{

  //Input data
  DATA_INTEGER(ndays);
  DATA_INTEGER(nstages);
  DATA_SCALAR(spacing);
  DATA_VECTOR(days);
  DATA_MATRIX(staging);
  DATA_MATRIX(stageprop);
  //DATA_INTEGER(nphi);
  DATA_VECTOR(phi1);
  DATA_VECTOR(phi2);
  //DATA_INTEGER(ncphi);
  DATA_VECTOR(cphi1);
  DATA_VECTOR(cphi2);
  DATA_VECTOR(cphi3);
  //DATA_INTEGER(ntm1);
  DATA_VECTOR(tm1);
  //DATA_INTEGER(ntn1);
  DATA_VECTOR(tn1);
  //DATA_INTEGER(ntn2);
  DATA_VECTOR(tn2);
  //DATA_VECTOR(ntn3);
  DATA_VECTOR(tn3);
  //DATA_INTEGER(nttot);
  DATA_INTEGER(datoind);
  DATA_VECTOR(ttot);


  //Catch level in future projections
  //Estimated parameters
  PARAMETER(pmub);
  PARAMETER(logsigmab);


  Type sigmab = exp(logsigmab);
  Type mub = Type(1)+(exp(pmub)/(Type(1.0)+exp(pmub)))*Type(29.0); //Bounded between 1-30

  int nphi = phi1.size();
  int ncphi = cphi1.size();
  int ntm1 = tm1.size();
  //int ntn1 = tn1.size();
  //int ntn2 = tn2.size();
  //int ntn3 = tn3.size();
  int nttot = ttot.size();



  matrix<Type> S(nstages,ndays);
  matrix<Type> P(nstages,ndays);
  matrix<Type> SPmatr(nstages,ndays);
  array<Type> Nout(nttot,4);

  vector<Type> m1(ntm1);
  vector<Type> m2(ntm1+nphi-1);
  vector<Type> m3(ntm1+2*nphi-2);
  vector<Type> n1(ntm1+ncphi-1);
  vector<Type> n2(ntm1+nphi+ncphi-2);
  vector<Type> n3(ntm1+2*nphi+ncphi-3);

  vector<Type> nn1(nttot);
  vector<Type> nn2(nttot);
  vector<Type> nn3(nttot);
  vector<Type> nn(nttot);

    //vector Nk(1,ndays)
    //vector positions(1,ndays)
    vector<Type> SPvect(ndays);
    vector<Type> b(2);
    //vector m2(1,nttot)
    //number ind
    //number counter

    Type PropEst;
    Type SN;
    Type S1;
    Type S2;
    Type S3;
    //Type tmp;


    //int i;
    //int k;
    int ind;
    int counter;
    Type pi=3.14159265358979323844;

    for(int i=0;i<ntm1;i++)
    {
      m1(i) = Type(1.0)/sqrt(Type(2.0)*pi*sigmab*sigmab)*exp(-0.5*((tm1(i)-mub)*(tm1(i)-mub))/(sigmab*sigmab));
    }

    m2 = conv(m1,phi1);
    m2 = m2*spacing;

    m3 = conv(m2,phi2);
    m3 = m3*spacing;

    n1 = conv(m1,cphi1);
    n1 = n1*spacing;

    n2 = conv(m2,cphi2);
    n2 = n2*spacing;

    n3 = conv(m3,cphi3);
    n3 = n3*spacing;

    ind = 0;
    while (ttot(ind) < tn1(0))
    {
      ind = ind+1;
    }
    ind = ind - 1;
    counter = ind;
    for(int k=0;k< n1.size();k++)
    {
      nn1(counter) = n1(k);
      counter = counter + 1;
    }

    ind = 1;
    while (ttot(ind) < tn2(0))
    {
      ind = ind+1;
    }
    ind = ind - 1;
    counter = ind;
    for(int k=0;k< n2.size();k++)
    {
      nn2(counter) = n2(k);
      counter = counter + 1;
    }

    ind = 1;
    while (ttot(ind) < tn3(0))
    {
      ind = ind+1;
    }
    ind = ind - 1;
    counter = ind;
    for(int k=0;k<n3.size();k++)
    {
      nn3(counter) = n3(k);
      counter = counter + 1;
    }

    nn = nn1+nn2+nn3;

    PropEst = nn(datoind);

    for(int k=0;k<ndays;k++)
    {
      ind = 1;
      while (ttot(ind) < days(k))
      {
        ind = ind+1;
      }
      ind = ind - 1;

      P(0,k) = nn1(ind)/nn(ind);
      P(1,k) = nn2(ind)/nn(ind);
      P(2,k) = nn3(ind)/nn(ind);

    }

    S = staging.transpose();

    //cout<<P<<endl;
    //cout<<stageprop<<endl;

    matrix<Type> Plog = log(P.array());
    SPmatr = S.array()*Plog.array();
    SPvect = SPmatr.colwise().sum();

  Type nll = Type(0.0);

  vector<Type> indv;
  for(int k=0;k<ndays;k++)
  {
    indv = staging.row(k);
    ind = CppAD::Integer(indv.sum());

    /*Type sumstaging;
    sumstaging = Type(0.0);
    for(int l=0;l<nstages ;l++){
      sumstaging = sumstaging + staging(k,l);
    }
    ind = CppAD::Integer(sumstaging);
    */

    counter = 1;
    SN = 0;
    for(int i=0;i<ind;i++)
    {
      SN = SN + log(counter);
      counter = counter + 1;
    }

    if (S(0,k) == 0) S1 = 0; else
    {
      //ind = (S(1,k));
      counter = 1;
      S1 = 0;
      for(int i=0;i<staging(k,0);i++)
      {
        S1 = S1 + log(counter);
        counter = counter + 1;
      }
    }

    if (S(1,k) == 0) S2 = 0; else
    {
      counter = 1;
      S2 = 0;
      for(int i=0;i<staging(k,1);i++)
      {
        S2 = S2 + log(counter);
        counter = counter + 1;
      }
    }

    if (S(2,k) == 0) S3 = 0; else
    {
      counter = 1;
      S3 = 0;
      for(int i=0;i<staging(k,2);i++)
      {
        S3 = S3 + log(counter);
        counter = counter + 1;
      }
    }

    nll = nll + (SN - S1 - S2 - S3 + SPvect(k));
  }

  for(int i=0;i<nttot;i++)
  {
    Nout(i,0) = nn1(i);
    Nout(i,1) = nn2(i);
    Nout(i,2) = nn3(i);
    Nout(i,3) = nn(i);
  }

  nll = -nll;

  ADREPORT(PropEst);
  ADREPORT(mub);
  ADREPORT(sigmab)
  REPORT(Nout);

  return nll;

}
