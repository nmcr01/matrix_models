
    # Simple matrix projection model for SIR process
    # Assume 3 states as in standard SIR model: Susceptible, Infected, Removed
    # There is no immigration or emigration and no mortality.  We assume that 
    # it is possible to define a time interval such that meaningful 
    # probabilities can be assumed for the state transitions.
    
    # Start by listing the probabiliites.  Note that this step is unnecessary for simple
    # models where we can just directly write the matrix, but it's good practice to
    # develop because it's easier to track the numerical values of probabilities in 
    # a list than in the block of code for the matrix.
    
    p_is <-0.25
    p_ss <-1-p_is
    p_ri <-0.3
    p_ii <-1-p_ri
    p_rr <-1
    
    # Now define a matrix to take the transition probabiliites
    
    A <-matrix(nrow=3, ncol=3, 
               data=c(p_ss,p_is,0,
                      0,p_ii,p_ri,
                       0,0,p_rr))
    # print A to make sure it has been constucted correctly
    A
    # Check the column entries of each column of A sum to 1
    # There are prettier, more efficient ways to implement
    # this for larger matrices.
    
    sum(A[,1])
    sum(A[,2])
    sum(A[,3])
    
    # Now we need a matrix to hold the output.  Declaring it
    # in advance is more efficient than trying to dynamically
    # add columns during calculation. Each column of the output
    # is a time step, so the width of this matrix defines the
    # total duration of the projection.
    
    projection <-matrix(nrow=3, ncol=21)
    
    # Now we define the initial states and assign them to the first columnn of
    # the output matrix
    
    initial <-c(10000,0,0)
    projection[,1]<-initial
    
    # The numerical projection of the SIR model is now performed using a 
    # "for" loop to iterate a matrix multiplication
    
    for (i in seq(2,21,1)) {
      projection[,i]<-A%*%projection[,i-1]
    }
    
    # Note that the iteration starts at time 2 and works by
    # calculating the "current" value from the immediate past
    # value.  That is, compared to the way we would typically
    # write the projection as a projection of the future, we are
    # calculating it with the pair of time index values slid back
    # by one step.
    
    # Make a vector to hold the timesteps and generate a plot of
    # the numerical output
    
    t<-seq(1,21,1)
    graph1<-plot(t,projection[1,], ty="l", lwd=3, col="darkblue",
                 xlab="time", ylab="individuals")
     lines(t,projection[2,], lwd=3, lty=2, col="red")
     lines(t,projection[3,], lwd=3, lty=3, col="black")
     legend("topright", inset=0.015, legend=c("S","I","R"),
     lty=c(1,2,3), col=c("darkblue","red","black"), lwd=c(3,3,3))
