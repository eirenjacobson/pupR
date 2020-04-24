
m = 3 #x
n = 2 #h

lowerlimit = rep(0,2)
upperlimit = rep(0,2)

y = rep(NA,n+m-1)

# for(k in 1:(m+n-1)){
#   cat("k = ",k,"\n")
#   lowerlimit[1]=1
#   lowerlimit[2]=k+1-n
#   upperlimit[1]=k
#   upperlimit[2]=m
#   for(j in max(lowerlimit):min(upperlimit)){
#     cat("    j = ",j,"\n")
#     cat("    y(",k,")\n")
#     cat("    x(",j,")\n")
#     cat("    h(",k-j+1,")\n")
#   }
#   cat("\n ---------------- \n")
# }


for(k in 1:(m+n-1)){
  cat("k = ",k,"\n")
  lowerlimit[1]=1
  lowerlimit[2]=k+1-n
  upperlimit[1]=k
  upperlimit[2]=m
  for(j in max(lowerlimit):min(upperlimit)){
    cat("    j = ",j,"\n")
    cat("    y(",k-1,")\n")
    cat("    x(",j-1,")\n")
    cat("    h(",k-j+1-1,")\n")
  }
  cat("\n ---------------- \n")
}

# for(k in 0:(m+n-2)){
#   cat("k = ",k,"\n")
#   lowerlimit[1]=0
#   lowerlimit[2]=k-n
#   upperlimit[1]=k-1
#   upperlimit[2]=m-1
#   for(j in max(lowerlimit):(min(upperlimit)-1)){
#     cat("    j = ",j,"\n")
#     cat("    y:",k,"\n")
#     cat("    x:",j,"\n")
#     cat("    h:",k-j+1,"\n")
#   }
#   cat("\n ---------------- \n")
# }

  # for(k=1;k<=(int) (m+n-1);k++)
  # {
  #   y(k) = 0;
  #   lowerlimit(1) = 1;
  #   lowerlimit(2) = k+1-n;
  #   upperlimit(1) = k;
  #   upperlimit(2) = m;
  #
  #   for(j=max(lowerlimit);j<=(int) min(upperlimit);j++)
  #   {
  #     y(k) = y(k) + x(j)*h(k-j+1);
  #   }
  # }
  #
