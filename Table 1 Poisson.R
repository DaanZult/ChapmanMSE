set.seed(42)
options(scipen=99)

r = 1000000
m11=20
n11=rpois(reps,m11)
mean(1/n11)


1/mean(n11)
1/mean(n11) - mean((n11-mean(n11))/mean(n11)^2)
1/mean(n11) - mean((n11-mean(n11))/mean(n11)^2) + mean((n11-mean(n11))^2/mean(n11)^3) 
1/mean(n11) - mean((n11-mean(n11))/mean(n11)^2) + mean((n11-mean(n11))^2/mean(n11)^3) - mean((n11-mean(n11))^3/mean(n11)^4) 
1/mean(n11) - mean((n11-mean(n11))/mean(n11)^2) + mean((n11-mean(n11))^2/mean(n11)^3) - mean((n11-mean(n11))^3/mean(n11)^4) + mean((n11-mean(n11))^4/mean(n11)^5)

1/mean(n11) - mean(1/n11)
1/mean(n11) - mean((n11-mean(n11))/mean(n11)^2) - mean(1/n11)
1/mean(n11) - mean((n11-mean(n11))/mean(n11)^2) + mean((n11-mean(n11))^2/mean(n11)^3) - mean(1/n11)
1/mean(n11) - mean((n11-mean(n11))/mean(n11)^2) + mean((n11-mean(n11))^2/mean(n11)^3) - mean((n11-mean(n11))^3/mean(n11)^4) - mean(1/n11)
1/mean(n11) - mean((n11-mean(n11))/mean(n11)^2) + mean((n11-mean(n11))^2/mean(n11)^3) - mean((n11-mean(n11))^3/mean(n11)^4) + mean((n11-mean(n11))^4/mean(n11)^5) - mean(1/n11)



mean(1/(n11+1)) 
mean(1/(n11+1)) + mean(1/((n11+1)*(n11+2)))
mean(1/(n11+1)) + mean(1/((n11+1)*(n11+2))) + mean(2/((n11+1)*(n11+2)*(n11+3)))
mean(1/(n11+1)) + mean(1/((n11+1)*(n11+2))) + mean(2/((n11+1)*(n11+2)*(n11+3))) + mean(6/((n11+1)*(n11+2)*(n11+3)*(n11+4)))
mean(1/(n11+1)) + mean(1/((n11+1)*(n11+2))) + mean(2/((n11+1)*(n11+2)*(n11+3))) + mean(6/((n11+1)*(n11+2)*(n11+3)*(n11+4))) + mean(24/((n11+1)*(n11+2)*(n11+3)*(n11+4)*(n11+5)))

mean(1/(n11+1)) - mean(1/n11)
mean(1/(n11+1)) + mean(1/((n11+1)*(n11+2))) - mean(1/n11)
mean(1/(n11+1)) + mean(1/((n11+1)*(n11+2))) + mean(2/((n11+1)*(n11+2)*(n11+3))) - mean(1/n11)
mean(1/(n11+1)) + mean(1/((n11+1)*(n11+2))) + mean(2/((n11+1)*(n11+2)*(n11+3))) + mean(6/((n11+1)*(n11+2)*(n11+3)*(n11+4))) - mean(1/n11)
mean(1/(n11+1)) + mean(1/((n11+1)*(n11+2))) + mean(2/((n11+1)*(n11+2)*(n11+3))) + mean(6/((n11+1)*(n11+2)*(n11+3)*(n11+4))) + mean(24/((n11+1)*(n11+2)*(n11+3)*(n11+4)*(n11+5))) - mean(1/n11)

