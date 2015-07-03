# Function to simulate an SIR model on a patch space 
ddv.run <- function(
		beta=0.3,
    predrop=0.35,
		B=0.2,
		Tg=2.6,
		D_R=365*4,
		alpha=5,
		pow=10,
    t_alpha=364*10,
    t_drop=354*5,
		offset=15,
		dxdy = 10,
		nrow=3,
		ncol=10,
		noreals = 1,
		tmax = 364*20,
		dt = 0.5,
		dtrep = 1,
		N = 1000000
) {
  
  # dtrep must be a multiple of dt
  if (dtrep %% dt > 1e-10) stop("dtrep must be a multiple of dt")
	
	# Setup state variables and incidence array
  nts <- floor(tmax/dt)
  ntrep <- min(max(nts*dt),floor(tmax/dtrep))
	inc <- array(dim=c(nrow,ncol,ntrep,noreals))
	st_sus <- array(dim=c(nrow,ncol))
	st_inf <- array(dim=c(nrow,ncol))
	st_rec <- array(dim=c(nrow,ncol))
  st_N <- array(dim=c(nrow,ncol))
	ev_inf <- array(dim=c(nrow,ncol))
	ev_rec <- array(dim=c(nrow,ncol))
	ev_bec <- array(dim=c(nrow,ncol))
  rt_inf <- array(dim=c(nrow,ncol))
	p_inf <- array(dim=c(nrow,ncol))
  per_seed <- array(dim=c(nrow,ncol))
  inc_tmp <- array(dim=c(nrow,ncol))
  
  # Setup distances and kernel values
  mat_dist <- array(dim=c(nrow,ncol,nrow,ncol))
  for (i in 1:nrow) {
    for (j in 1:ncol) {
      for (k in 1:nrow) {
        for (l in 1:ncol) {
          mat_dist[i,j,k,l] <- dxdy*sqrt((i-k)^2+(j-l)^2)   
        }
      }
    }
  }
  mat_kappa <- 1/(1+(mat_dist/offset)^pow)
  
  p_rec <- 1-exp(-1/Tg)
  p_bec <- 1-exp(-1/D_R)
	
	# preconditions for the reals loop
	inc[] <- 0
	
	# Initiate realizations
	count_reals <- 1
  while (count_reals <= noreals) {
		
		# Initialize the state variables for this realization
		st_sus[] <- round(N / (nrow*ncol))
		st_inf[] <- 0
		st_rec[] <- 0
    st_N <- st_sus + st_inf + st_rec
    inc_tmp[] <- 0
		currenttime <- 0
    per_seed[] <- 0
    per_seed[1,1] <- alpha / st_N[1,1] / 364
		
    # Define a lambda function
		func_lambda <- function(x){sum(st_inf*x/st_N)}
    
    # Start the time loop
    count_times <- 1
    count_reps <- 1
    while (count_times <= nts) {
      
      # Generate the force of infection for each population
      # at the current time
      beta_t <- ifelse(
        count_times * dt > t_drop,
        beta * (1+ B * sin(2*pi/364*count_times*dt)),
        predrop * beta * (1+ B * sin(2*pi/364*count_times*dt))
      )
      rt_inf <- beta_t * apply(mat_kappa,c(1,2),func_lambda) + per_seed
      for (i in 1:nrow) {
        for (j in 1:ncol) {
          # rt_inf[i,j] <- beta_t * sum(st_inf * mat_kappa[i,j,,]) / st_N[i,j] + per_seed[i,j]          
        }
      }
           
      # Calculate the probability of an infection event
      p_inf[] <- 1-exp(-rt_inf * dt)
      
      # Generate infection events
      # if (min(c(st_sus,st_inf,st_rec)) < -0.1) browser()
      if (sum(is.na(c(st_sus,st_inf,st_rec))) > 0) browser()
      for (i in 1:nrow) {
        ev_inf[i,] <- rbinom(ncol,st_sus[i,],p_inf[i,])
        ev_rec[i,] <- rbinom(ncol,st_inf[i,],p_rec)
        ev_bec[i,] <- rbinom(ncol,st_rec[i,],p_bec)
      }
      
      # Apply the actual events
      st_sus  <- st_sus   - ev_inf + ev_bec
      st_inf  <- st_inf   + ev_inf - ev_rec
      st_rec  <- st_rec   + ev_rec - ev_bec
      inc_tmp <- inc_tmp  + ev_inf 
      
      # Save to reported incidence
      if (abs(count_times*dt - count_reps*dtrep) < 1e-10) {
        inc[,,count_reps,count_reals] <- inc[,,count_reps,count_reals] + inc_tmp
        inc_tmp[]   <- 0
        count_reps  <- count_reps + 1
      }
      
      # Close the timestep loop
      count_times <- count_times +1      

    }
    
		# Close the realization loop
		count_reals <- count_reals + 1
	}
	
	list(inc=inc,dtrep=dtrep,dxdy=dxdy,N=st_sus+st_inf+st_rec)
	
}

ddv.plot.ts <- function(x,log="") {
  
  # Decompse the 
  plot(
    x$dtrep*(1:(dim(x$inc)[3]))/364,
    # (apply(x$inc[,,,1],c(3),sum)+1)/x$dtrep/sum(x$N),
    (apply(x$inc[,,,1],c(3),sum)+1)/x$dtrep/sum(x$N)*100,
    type="l",
    log=log,
    xlab="Time",
    ylab="Incidence (+1, log scale)")
  
}