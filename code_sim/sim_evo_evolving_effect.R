# Simulate the origin of a novel character identity under a sequential-fixation regime
# Mutations affecting regulator expression & mutations affecting cis-trans interactions are both allowed

setwd("/your_directory/")

# Function to calculate fitness (Gaussian function of a multivariate phenotype with w_max=1)
fitness_calc<-function(z,opt,omega) {
  Lambda<-solve(omega)
  D2<-t(z-opt)%*%Lambda%*%(z-opt)
  w<-exp(-D2/2)
  return(as.numeric(w))
}

# Function to calculate fixation probability
# Parameters: ancestral fitness (wa), mutant fitness (wm), effective population size (Ne)
fix_prob_calc <- function(wa,wm,Ne){
  if(wa>0){
    s=(wm/wa)-1 # Coefficient of selection
    if(s==0){
      p=1/(2*Ne) # Neutral mutation
    }else{
      p=(1-exp(-2*s))/(1-exp(-4*Ne*s))
    }
  }else{ # The ancestral fitness is close to 0 (close enough to be recognized as zero by R)
    if(wm==0){ # Both ancestral and mutant fitness are close to 0
      p=1/(2*Ne) # The mutation is considered neutral
    }else{ # Ancestral fitness is close to 0 while mutant fitness isn't
      p=1 # The mutation is considered strongly beneficial
    }
  }
  return(p)
}

# Function to simulate evolution
# Arguments: ancestral expression levels of regulators (X_ances, a vector), per-unit-expression phenotypic effect of regulators (effect_all, a matrix), target state (z_opt, a vector), adaptive landscape (omega, a covariance matrix), 
# Arguments (continued): effective population size (Ne), simulation time (n_step), mutation rate per genome per time step (U), mutation effect size SD (sigma)
evo_sim_plus <- function(X_ances,effect_all_ances,z_opt,omega,Ne,n_step,U,sigma,U_cis,sigma_cis){
  n_dim=nrow(effect_all_ances);n_regulator=ncol(effect_all_ances)
  X=X_ances # 'Initialize' X
  effect_all=effect_all_ances # 'Initialize' effect matrix
  n_mut=rpois(1,lambda=2*Ne*(U+U_cis)*n_step) # Total number of mutations to occur
  if(n_mut>0){
    z=effect_all%*%log(X_ances)
    w=fitness_calc(z,z_opt,omega)
    type_mut=rbinom(n=n_mut,size=1,prob=U/(U+U_cis)) # Whether each mutation affects regulator expression or cis-trans regulation
    for(n in 1:n_mut){
      if(type_mut[n]==1){
        X_mut=X
        affect=sample(1:length(X),1) # Each mutation affects a single gene's expression
        effect=rnorm(1,mean=0,sd=sigma) # Effect of mutation (log scale)
        X_mut[affect]=X_mut[affect]*exp(effect) # Mutant expression level
        z_mut=effect_all%*%log(X_mut) # Mutant phenotype
        w_mut=fitness_calc(z_mut,z_opt,omega) # Mutant fitness
        fp=fix_prob_calc(w,w_mut,Ne) # Fixation probability
        fix=rbinom(n=1,size=1,prob=fp) # Whether the mutations fixes
        if(fix==1){
          X=X_mut;z=z_mut;w=w_mut # Pass the mutant's genotype & phenotype down to next step if the mutation is fixed
        }
      }else{
        effect_all_mut=effect_all
        affect=c(sample(1:n_dim,1),sample(1:n_regulator,1)) # Each mutation affects a single gene's effect on a single trait
        effect=rnorm(1,mean=0,sd=sigma) # Effect of mutation (additive)
        effect_all_mut[affect[1],affect[2]]=effect_all_mut[affect[1],affect[2]]+effect
        z_mut=effect_all_mut%*%log(X) # Mutant phenotype
        w_mut=fitness_calc(z_mut,z_opt,omega) # Mutant fitness
        fp=fix_prob_calc(w,w_mut,Ne) # Fixation probability
        fix=rbinom(n=1,size=1,prob=fp) # Whether the mutations fixes
        if(fix==1){
          effect_all=effect_all_mut;z=z_mut;w=w_mut
        }
      }
    }
  }
  return(list(X,effect_all))
}

# Read pre-generated matrices
effect_all<-read.table("data_effect_all.txt",sep="\t");effect_all=as.matrix(effect_all)
n_regulator=ncol(effect_all);n_dim=nrow(effect_all)

# Set ancestral state
X_ances=rep(1,n_regulator) # An ancestral state where everything is lowly expressed
z_ances=effect_all%*%X_ances # Calculate ancestral character state

# Simulation-specific parameters
n_selection=5 # The number of phenotypic dimensions under selection
# Read-pre-existing list of traits under selection
fn0=paste("selection_list_",n_selection,".txt",sep="")
sub_selection<-read.table(fn0,sep="\t")
sub_selection=sub_selection[,1]
sub_selection
effect_sub=effect_all[sub_selection,];z_ances_sub=z_ances[sub_selection] # Extract rows corresponding to traits under selection

# Optimal phenotype
z_div_scale=2 # Distance between the optimum for each trait to the ancestral state (assumed to be the same for all traits for now)
z_opt=z_ances[sub_selection]+z_div_scale # Optimum
omega=matrix(0,nrow=n_selection,ncol=n_selection);diag(omega)=100 # Adaptive landscape shape (no correlational selection for now)

# Set parameters for simulation function
Ne=1000 # Effective population size
n_step=2e5 # Duration of simulation
U=n_regulator*5e-6 # Total mutation rate
sigma=0.1 # SD of mutational effect on each regulator's expression (assumed to be the same for all regulators)
U_cis=n_selection*n_regulator*1e-7;sigma_cis=.1 # Rate and effect size distribution of mutations affecting phenotypic effects of regulators (each mutation affects one regulator's effect on one trait)

# Simulate
n_rep=20 # Number of replicate lineages (species) per parameter combination
X_out_all=matrix(nrow=n_rep,ncol=n_regulator)
z_out_all=matrix(nrow=n_rep,ncol=n_selection)
fitness_out_all=rep(0,n_rep)
for(n in 1:n_rep){
  X_out=evo_sim(X_ances,effect_sub,z_opt,omega,Ne,n_step,U,sigma)
  z_out=effect_sub%*%log(X_out)
  fitness_out=fitness_calc(z_out,z_opt,omega)
  X_out_all[n,]=X_out
  z_out_all[n,]=z_out
  fitness_out_all[n]=fitness_out
  cat(paste(n,"\t",fitness_out,"\n",sep=""))
}

#fn0=paste("selection_list_",n_selection,".txt",sep="") # The set of traits under selection
fn1=paste("sim_out_X_",n_selection,".txt",sep="") # Regulator expression programs at the end
fn2=paste("sim_out_z_",n_selection,".txt",sep="") # Phenotypic states at the end
fn3=paste("sim_out_fitness_",n_selection,".txt",sep="") # Fitness at the end
#write.table(sub_selection,file=fn0,sep="\t")
write.table(X_out_all,file=fn1,sep="\t")
write.table(z_out_all,file=fn2,sep="\t")
write.table(fitness_out_all,file=fn3,sep="\t")




