# Generate a data matrix that contains the effect of each regulator on each phenotypic trait

setwd("/your_directory/")

n_regulator=50 # The number of regulators available
n_dim=1000 # The total number of phenotypic dimensions (not all under selection)
effect_all=matrix(0,nrow=n_dim,ncol=n_regulator) # Data matrix containing all regulators' phenotypic effects per unit expression
# Effects of regulators on a random trait
effect_size_base_1=rexp(n_regulator/2,rate=2) # Sample absolute values of effects from an exponential distribution
effect_size_base_2=-effect_size_base_1
effect_size_base=c(effect_size_base_1,effect_size_base_2) # Sum to zero such that there is no mutational bias at phenotypic level
for(n in 1:n_dim){
  effect_all[n,]=signif(sample(effect_size_base,n_regulator),5) # Shuffle the same vector for each trait, such that every trait has the same amount of mutational variance
}
write.table(effect_all,file="data_effect_all.txt",sep="\t") # Save data for future use
