#' Calculate the modelled fitting to the staging data
#'
#' Calculate the modelled fittet proportions to the observed staging data given a birth distribution
#' @param mub Mean of the birth distribution
#' @param sigmab Sigma for det birth distribution
#' @param data The data set with photo counts
#' @param type Type of correction. Use 1 for separate correction for each reader (assumes two readers) and 2 for the same correction for both readers
#' @return
#' @keywords
#' @export
#' @examples
#' propFitFun()

propFitFun <- function(mub = NA,
                       sigmab = NA,
                       data = NA,
                       timeDF = NA)
{

  # Input parameters
  spacing = data$spacing

  tm1 = data$tm1
  phi1 = data$phi1
  phi2 = data$phi2
  phi3 = data$phi3
  cphi1 = data$cphi1
  cphi2 = data$cphi2
  cphi3 = data$cphi3

  tau = timeDF$tau
  tm1 = timeDF$tm1
  tm2 = timeDF$tm2
  tm3 = timeDF$tm3
  tn1 = timeDF$tn1
  tn2 = timeDF$tn2
  tn3 = timeDF$tn3
  t_min = timeDF$t_min
  t_max = timeDF$t_max
  t_tot = timeDF$t_tot

  nn2 = array(0,length(t_tot));
  nn1 = nn2;
  nn3 = nn2;

  m1 = dnorm(tm1,mub,sigmab)
  m2 = convolve(m1,rev(phi1),type = "op")*spacing;
  m2 = m2 / (trapz(1:length(m2),m2)*spacing);
  m3 = convolve(m2,rev(phi2),type = "op")*spacing;
  m3 = m3 / (trapz(1:length(m3),m3)*spacing);


  n1 = convolve(m1,rev(cphi1),type = "op")*spacing;
  n2 = convolve(m2,rev(cphi2),type = "op")*spacing;
  n3 = convolve(m3,rev(cphi3),type = "op")*spacing;

  pos = min(which(t_tot>tn1[1]));
  nn1[pos:(pos+length(n1)-1)] = n1;
  pos = min(which(t_tot>=tn2[1]));
  nn2[pos:(pos+length(n2)-1)] = n2;
  pos = min(which(t_tot>=tn3[1]));
  nn3[pos:(pos+length(n3)-1)] = n3;

  nn = nn1+nn2+nn3;


  returnList = list()
  returnList$nn1 = nn1
  returnList$nn2 = nn2
  returnList$nn3 = nn3
  returnList$nn = nn

  }


############################


  tau = seq(0,10,by = spacing)

  phi_1 = dgamma(tau[-1],kappa,1/rho_1);
  phi_2 = dgamma(tau[-1],kappa,1/rho_2);
  phi_3 = dgamma(tau[-1],kappa,1/rho_3);


  S = t(S)

  tau = seq(spacing,30,by = spacing)
  t = seq(0,15,by = 1/24)
  tm1 = seq(spacing,40,by = spacing)

  #lognormal birth distr
  # m1 = plnorm(tm1,mu_b,sigma_b);
  #Gamma birth distr
  # m1 = (tm1-mu_b).^(kappa_b-1).*exp(-(tm1-mu_b)/rho_b)/rho_b^kappa_b/gamma(kappa_b);
  # %Gaussian birth distr
  #m1 = 1/sqrt(2*pi*sigma_b^2)*exp( -1/2/sigma_b^2*(tm1-mu_b).^2 );
  m1 = dnorm(tm1,mu_b,sigma_b)


  m2 = convolve(m1,rev(phi_1),type = "op")*spacing;
  m2 = m2 / (trapz(1:length(m2),m2)*spacing);
  tm2 = seq(tau[2]+tm1[1],tau[length(tau)]+tm1[length(tm1)],length = length(m2));
  m3 = convolve(m2,rev(phi_2),type = "op")*spacing;
  m3 = m3 / (trapz(1:length(m3),m3)*spacing);
  tm3 = seq(tau[2]+tm2[1],tau[length(tau)]+tm2[length(tm2)],length = length(m3));

  cphi_1 = 1-pgamma(tau,kappa,1/rho_1);
  cphi_2 = 1-pgamma(tau,kappa,1/rho_2);
  cphi_3 = 1-pgamma(tau,kappa,1/rho_3);

  n1 = convolve(m1,rev(cphi_1),type = "op")*spacing;
  tn1 = seq(tm1[1]+tau[1],tm1[length(tm1)]+tau[length(tau)], length = length(n1));

  n2 = convolve(m2,rev(cphi_2),type = "op")*spacing;
  tn2 = seq(tm2[1]+tau[1],tm2[length(tm2)]+tau[length(tau)], length = length(n2));

  n3 = convolve(m3,rev(cphi_3),type = "op")*spacing;
  tn3 = seq(tm3[1]+tau[1],tm3[length(tm3)]+tau[length(tau)], length = length(n3));

  t_min = min(c(tn1[1],tn2[1],tn3[1]))-1;
  t_max = max(c(tn1[length(tn1)],tn2[length(tn2)],tn3[length(tn3)]))+1;
  t_tot = seq(t_min,t_max,by = spacing);
  nn2 = array(0,length(t_tot));
  nn1 = nn2;
  nn3 = nn2;

  pos = min(which(t_tot>tn1[1]));
  nn1[pos:(pos+length(n1)-1)] = n1;
  pos = min(which(t_tot>=tn2[1]));
  nn2[pos:(pos+length(n2)-1)] = n2;
  pos = min(which(t_tot>=tn3[1]));
  nn3[pos:(pos+length(n3)-1)] = n3;

  nn = nn1+nn2+nn3;

  P = matrix(NA,3,length(positions))
  P[1,] = nn1[positions]/nn[positions];
  P[2,] = nn2[positions]/nn[positions];
  P[3,] = nn3[positions]/nn[positions];

  SPmatr = S*log(P);

