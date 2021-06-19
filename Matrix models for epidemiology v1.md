## Introduction
As the title suggests, in this chapter we will take a look how matrices can be used to build dynamic models of epidemics.  In contrast to most of the the approaches discussed in the book, which are based on the classical approach of treating time as a continuous variable, the matrix models we will be working with in this chapter consider discrete events occurring along a series of time steps that can be thought of as an ordered index.  In addition, most of the models do not include (in any direct way) density-dependence or non-linear components.  Despite these limitations we will see that matrix models are  able to capture a range of useful settings in which some kind of quantitative framework is needed for understanding disease dynamics. Indeed, we will see that when combined with the process of drawing a state transition diagram for the system we want to study, matrix models offer a powerful and adaptable way for non-experts to get involved in the process of model building, because of the homology between the representations of the system in the state diagram and the projection matrix.   Using this approach, modeling starts with a simple picture of the system being studied, and adds the mathematical and computational aspects in a logical way. 

The state diagram is simply a picture of the states in the system and how they are connected by possible transitions one state to another (_e.g._ whether susceptible, healthy individuals can become exposed to a pathogen).  The projection matrix takes the set of states and the set of transitions and captures them in a formal structure that guarantees only the transitions shown in the diagram occur when we calculate numerical results.  The rules for translating the diagram to the matrix are easy to learn, facilitating the process of model building.  We labour these points because to achieve real-world impact there is no question that epidemiology has to be accessible to the audiences whose problems we want to solve.  Matrix models will sometimes not offer the appropriate match between biology and mathematics needed to solve complex problems, but they certainly provide an access point for non-experts that may greatly help the exchange of ideas between stakeholders and scientists.  Simple "toy" models implemented in matrix format may provide a focus for discussion that will help to identify the specifications for something more complex to be implemented in a different format.

If you are not familiar with matrix models, and this brief introduction whets your appetite to learn more, you are in luck.  There is an outstanding single reference text that will likely tell you everything (and much more) that you will want to know about matrix models.  That text is Hal Caswell's [Matrix Population Models: Construction, Analysis, and Interpretation](https://global.oup.com/ushe/product/matrix-population-models-9780878931217?q=Caswell&cc=us&lang=en). The chapter on matrix models in Henry Stevens', [A Primer of Ecology with R](https://www.springer.com/gp/book/9780387898810) from the Use R! series gives a useful, short introduction in the context of the R language.  In keeping with the rest of this book we will be implementing all of our examples in R, but in the spirit of highlighting accessibility for matrix modeling we should note that popular spreadsheet programs include a basic function for matrix multiplication, including free, cloud-based applications such as Google Sheets, so for those who want to try modeling but don't want to learn any coding the possibility of matrix modeling still exists.

## A simple Susceptible-Infected-Removed (SIR) model

### The State Diagram 
In figure 1 you can see what will probably be a familiar diagram; three bubbles (or _nodes_) labelled _S_, _I_, and _R_.  These represent three possible _states_ that an individual can have in many different types of disease epidemic.  Note that the states are drawn in sequence from left to right and are connected by arrows (or directed _edges_) that indicate a flow from one state to the next. The susceptible (_S_) state is usually taken to be synonymous with healthy, uninfected, or not exposed. Individuals in this state are available to become infected.  In the simplest _SIR_ models (such as the one we're dealing with here) susceptible (healthy) individuals move directly to being infected (_I_) without any intermediate exposed state or state of latent infection.  Infected individuals can become removed (_R_) from the epidemic.  The _R_ can be taken to mean _recovered_ in some situations.  Unfortunately, there are situations where removed has unhappy connotations, but it is a more general term, indicating that individuals (for whatever reason) are no longer participating in the spread of the disease so we will continue to think of _R_ as removed.  Note that in what has been said so far we have referred to individuals moving among states.  Thinking in terms of individuals obviously has a natural interpretation in human disease epidemiology, but it can also naturally apply to leaves, plants, fields or any other organizational unit.  The important point here is that the model applies to some kind of discrete unit and various kinds of events result in these units moving among a defined set of mutually exclusive states.

The question that naturally arises now is how, in the model, does the change from one state to another occur?  To answer that question, recall that we are thinking in terms of discrete time steps of equal length.  At the beginning of the disease process, we can ask ourselves, given a time step of duration $\delta$_t_ what is the probability that a susceptible (_S_) individual will become infected (_I_)? A similar question can be asked about the probability that an infected individual (_I_) will be removed (_R_).  The two probabilities are represented by the _edges_ in the state transition diagram. The edges originating from a node represent all the possible events that can be experienced by individuals in that state during a time step.  Since the edges represent probabilities, and the events they represent are mutually exclusive and exhaustive they must sum to 1.

Unless all susceptible individuals become infected and all infected individuals are removed each time step there must some individuals of each type who **do not** change state.  For completeness, if we want to represent the chance that some individuals do not change state in a time step, we add "self loops" to the states, as shown in Figure 2.  It is often the case in the process of model building that the value given to the self loop is found after we have specified the values on all the other edges.  We know from the law of total probability that the value on the self loop must be 1 minus the sum of the values on the other edges. This attribute of the modeling approach also gives us mechanism for making sure that we are accounting for all possible fates for individuals in each state, when we translate the state diagram to its corresponding projection matrix (see below).

**Quick recap** The state diagram represents the possible states for individuals in the population we are modeling.  We consider fixed, discrete time intervals and think about the dynamic processes in terms of the probability that individuals will change state in a single time step.  Choices we have to make are what the individuals are (what scale are we modeling as the fundamental unit of observation?) and what period of time is represented by the time step.

### The Projection Matrix
The state diagram is a static visual model of the system we are studying.  If we want to look at the dynamics of the system we have to translate the diagram into a suitable format to perform calculations or conduct algebraic analysis.  To achieve the desired translation we are going to exploit the fact that for any set of $n$ states there are $n^{2}$ possible transitions (assuming we count forward and backward transitions separately and allow for self-loops).  Thus, whatever transitions actually occur in reality, we can capture them in an $n \times n$ square matrix in which the rows and columns of the matrix represent the set of $n$ states.  In the current example there are three states _S_, _I_, and _R_, so the set of all possible transitions can be represented in matrix format using a 3 by 3 square matrix.  We can identify the matrix by any letter or symbol we choose.  In some conventions matrices are named with upper case English letters written in bold typeface.  We follow those conventions here.  For the simple _SIR_ model we might have something such as$$
\mathbf{A=}\begin{pmatrix}
(1-pIS) & 0 & 0 \\
pIS & (1-pRI) & 0 \\
0 & pRI & 1 \\
\end{pmatrix}$$
in which the rows and columns follow the same order as the states in the _SIR_ model.  We are going to use the matrix as a multiplier to recursively change a vector in which the individual elements are the numbers (or proportions) of individuals in each of the three states, _S_, _I_, and _R_.  To do this we will use standard matrix multiplication.  You might remember from classes on linear algebra that matrix multiplication follows strict rules in which the number of columns in the matrix on the left of the product must be equal to the number of rows on the right of the product.

In theory we're free to choose which orientation we use to set up the projection matrix and the vector of states. Different disciplines tend to follow one arrangement or the other.  In population ecology and demography the convention is to write the projection matrix as shown above, so that the **columns** show the possible fates for individuals in each state according to the transitions indicated in the state diagram.  For example, the two non-zero entries in the first column show that individuals in _S_ can remain in _S_, with probability (1-_pIS_) or they can move to _I_, with probability _pIS_.  Since transition from _S_ directly to _R_ is not possible, the entry in the third row of the first column is a 0 and there is no corresponding edge in the state diagram.

Recalling that we mentioned formatting the model in this way provides a convenient error check on our working, we note that the values in each column of the matrix should sum to one.  This is an easy calculation to add to the program code or spreadsheet during the process of model building.

### Projecting the state vector
Before we look at the matrix equation we will make a slight detour to consider the analogous case of population growth in which a population increases (or decreases) by a constant proportion each time step.  We might represent such a process using a simple equation such as:$$
n_{t+1}=\rho \cdot n_{t}$$
in which $\rho$ is a constant.  When $\rho >1$ the population will grow exponentially and when $\rho<1$ the population will shrink exponentially with time.  In a loose sense the projection matrix can be thought of as a multi-value replacement for $\rho$ in a difference equation describing the simultaneous change in several states which together comprise the whole population:$$
\begin{pmatrix}
S_{t+1}\\
I_{t+1}\\
R_{t+1}\\
\end{pmatrix} =
\begin{pmatrix}
(1-pIS) & 0 & 0 \\
pIS & (1-pRI) & 0 \\
0 & pRI & 1 \\
\end{pmatrix} \cdot
\begin{pmatrix}
S_{t}\\
I_{t}\\
R_{t}\\
\end{pmatrix}$$

Or, in compact form:$$
\mathbf{n}_{t+1}=\mathbf{A} \cdot \mathbf{n}_{t}$$

Note that the entry in the third row, third column of $\mathbf{A}$ is 1, corresponding to the idea, captured in the state diagram, that individuals who reach _R_ remain there. Nodes such as _R_ in the _SIR_ model, from which there is no exiting transition are known as _absorbing states_ because they capture (or absorb) individuals flowing into them.  If a system has one or more absorbing states given enough time all of the individuals in the population will end up in the absorbing state(s).

Once the values for the transition probabilities are known, it is a relatively trivial task in R to produce a numerical projection for a vector specifying the numbers of individuals (or proportions of the population) in each state at each time.  The calculations only need features included in the base R installation, the most important being the matrix multiplication operator indicated by the %*% symbol combination.  The calculation makes use of a _for_ loop to recursively calculate successive instances of the state vector from the previous entry and the projection matrix.

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

