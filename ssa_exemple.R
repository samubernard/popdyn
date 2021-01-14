# Script SSA
# n1: cellules souches
# n2: cellules matures
# 1 mort de CS: mu1*n1
# 2 mort de M: mu2*n2
# 3 replication de CS: beta*n1*(1-n1/K)*(n1 <= K)
# 4 differenciation de CS -> M: f0*theta^h/(theta^h + n2^h)*n1

# initialisation
n1 <- 10
n2 <- 5

# parametres (jours)
mu1 <- 0.01
mu2 <- 0.1
beta <- 1.0
K <- 100
f0 <- 0.8 # f0 < beta
theta <- 50
h <- 2 # entre 1 et 4

tfinal <- 100.0
t <- 0
nbrIter <- 0
ntot <- n1+n2

while ( (t[nbrIter+1] < tfinal) && (nbrIter < 10000) && ntot > 0)
{
  nbrIter <- nbrIter + 1
  liste_taux <- c(mu1*n1[nbrIter], 
                  mu2*n2[nbrIter], 
                  beta*n1[nbrIter]*(1-n1[nbrIter]/K)*(n1[nbrIter] <= K),
                  f0*theta^h/(theta^h + n2[nbrIter]^h)*n1[nbrIter]) 
  somme_taux <- sum(liste_taux)
  tau <- rexp(1,rate = somme_taux)
  c <- runif(1, min = 0, max = 1)
  s <- cumsum(liste_taux)/somme_taux
  index_evenement <- min(which(c < s))
  if (index_evenement == 1) # mort n1
  {
    n1 <- c(n1,n1[nbrIter] - 1)
    n2 <- c(n2,n2[nbrIter])
  }
  else if ( index_evenement == 2 )  # mort n2
  {
    n2 <- c(n2,n2[nbrIter] - 1)
    n1 <- c(n1,n1[nbrIter])
  }
  else if ( index_evenement == 3 ) # replication
  {
    n1 <- c(n1,n1[nbrIter] + 1)
    n2 <- c(n2,n2[nbrIter])
  }
  else if ( index_evenement == 4 ) # differenciation
  {
    n1 <- c(n1,n1[nbrIter] - 1)
    n2 <- c(n2,n2[nbrIter] + 1)    
  }
  t <- c(t,t[nbrIter] + tau)
  ntot <- n1[nbrIter+1]+n2[nbrIter+1]
  
}
plot(t,n1, type = 'l', xlim = c(0,100), ylim = c(0,150), 
     xlab = "temps", ylab = "n1, n2")
lines(t,n2, col = "red")
