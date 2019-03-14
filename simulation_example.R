library(bipartite)
connectance=.1	# True connectance (proportion of possible links that actually occur)
n_plants=50		# Number of plants (or other group A)
n_gallers=50 	# Number of gallers (or other group B)
n_samples=30 # Number of sampling days


# Create a random network with desired size and connectance
Make_A<-function(n_plants,n_gallers,connectance){
	# Randomly assign links among species 
	links=c(rep(1,connectance*n_plants*n_gallers),rep(0,(1-connectance)*n_plants*n_gallers))
	# Create a matrix with links randomly assigned to give our desired connectance
	A=matrix(nrow=n_plants,ncol=n_gallers,data=0)
	# At this stage, all plants and gallers must have at least one interaction
	zero_r=length(which(rowSums(A)==0))
	zero_c=length(which(colSums(A)==0))
	while(zero_r+zero_c>0){
		A=matrix(nrow=n_plants,ncol=n_gallers,data=sample(links,replace=FALSE))
		zero_r=length(which(rowSums(A)==0))
		zero_c=length(which(colSums(A)==0))
	}	
	return(A)
}


# Create a network based on A that incorporates process uncertainty
# Assume the probability of an interaction occurring during a given sampling day,
# given the two species co-occur, is uniformly distributed
# Assume the probability of two species co-occurring is also uniformly distributed.
# Assume the probability of each interaction is independent
Make_B<-function(A,n_samples){
	B=matrix(nrow=nrow(A),ncol=ncol(A))
	for(plant in 1:nrow(A)){
		for(galler in 1:ncol(A)){
			process=runif(1) # Choose a random process uncertainty for each interaction
			cooccur=round(runif(n=1,min=0,max=n_samples))
			# Assuming the interaction was feasible, how often would it be observed based on the number of co-occurrences?
			n_local_occurrances=rbinom(1,cooccur,process)
			# If the interaction was feasible and occurred at least once, it is included in B
			if(A[plant,galler]==1 && n_local_occurrances>0){
				B[plant,galler]=1
			} else {
				B[plant,galler]=0
			}
		}
	}
	return(B)
}

# Create a network based on B that incorporates detection uncertainty
# Assume that a few interactions are easily detectable but many are difficult to detect
# Detection probability might therefore follow a beta distribution with parameters 2 and 10

# Assume that the number of co-occurrences for each species pair is uniformly distributed
# Assume that the probability of detecting an interaction, for a given number of co-occurrences, is beta distributed
# The overall detection uncertainty is a binomial process depending on the number of co-occurrences and detection probability for each interaction

Make_C<-function(B,n_samples){
	C=matrix(nrow=nrow(B),ncol=ncol(B))
	for(plant in 1:nrow(B)){
		for(galler in 1:ncol(B)){
			cooccur=round(runif(n=1,min=0,max=n_samples)) # Number of co-occurrences is a random integer
			detection=rbeta(1,2,10) # Intrinsic probability of detecting an interaction is uniformly distributed
			detections_if_occurring=rbinom(1,cooccur,detection) # Is the interaction detectable?
			# If the interaction occurred during sampling and is detectable in at least one sample, it is included in C
			if(B[plant,galler]==1 && detections_if_occurring>0){
				C[plant,galler]=1
			} else {
				C[plant,galler]=0
			}
		}
	}
	return(C)
}


test_network<-function(net){
	# Ignore species not observed interacting at least once
	trim_rows=net[which(rowSums(net)>0),]
	trimmed=trim_rows[,which(colSums(trim_rows)>0)]
	# Calculate properties
	n_plants=nrow(trimmed)
	n_gallers=ncol(trimmed)
	C=sum(trimmed)/length(trimmed)
	mean_L_per_p=mean(rowSums(trimmed))
	mean_L_per_g=mean(colSums(trimmed))
	NODF=nestedness(trimmed,null.models=FALSE)$temperature

	properties=c(n_plants,n_gallers,C,mean_L_per_p,mean_L_per_g,NODF)
	return(properties)
}

# Some matrices to store network properties
# A_props=matrix(nrow=10,ncol=2)
# colnames(A_props)=c("L_per_plant","L_per_galler")
# B_props=matrix(nrow=100,ncol=5)
# colnames(B_props)=c("N_plants","N_gallers","C","L_per_plant","L_per_galler")
# C_props=matrix(nrow=1000,ncol=5)
# colnames(C_props)=c("N_plants","N_gallers","C","L_per_plant","L_per_galler")
A_props=matrix(nrow=100,ncol=3)
colnames(A_props)=c("L_per_plant","L_per_galler","NODF")
B_props=matrix(nrow=10000,ncol=6)
colnames(B_props)=c("N_plants","N_gallers","C","L_per_plant","L_per_galler","NODF")
C_props=matrix(nrow=1000000,ncol=6)
colnames(C_props)=c("N_plants","N_gallers","C","L_per_plant","L_per_galler","NODF")
for(a in 1:100){
	print(a)
	A=Make_A(n_plants,n_gallers,connectance)
	properties=test_network(A)
	A_props[a,]=properties[4:6]
	# A_props[a,]=properties[4:5]
	for(b in 1:100){
		B=Make_B(A,n_samples)
		properties=test_network(B)
		rownum=b+(a-1)*100
		B_props[rownum,]=properties
		for(c in 1:100){
			C=Make_C(B,n_samples)
			properties=test_network(C)
			rownum2=c+(b-1)*100+(a-1)*10000
			C_props[rownum2,]=properties
		}
	}
}

# save.image('network_properties.Rdata')
