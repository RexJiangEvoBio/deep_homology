# Let lineages undergo prior expression divergence before the onset of selection
# Representing state 1 divergence in 2-staged simulations

setwd("/your_directory/")

n_regulator=50

# Random divergence of regulator expression
n_rep=20 # Number of replicate lineages (species)
sigma_x=0.1 # Evolutionary rate (assumed to be equal for all)
T_0=100 # Time
X_ances_new=matrix(0,nrow=n_rep,ncol=n_regulator)
for(i in 1:n_rep){
  delta=rnorm(n_regulator,mean=0,sd=sigma_x*sqrt(T_0)) # Amount of evolutionary change (log scale)
  X_ances_new[i,]=X_ances*exp(delta) # Expression at the end of stage 1
}
write.table(X_ances_new,file="X_ances_diverged.txt",sep="\t") # Write file
