# Dynamic Control Law
This method is a nonlinear control algorithm that solves optimal control problem offline through solving an algebraic inquality. Its performance is much better than LQR and SDRE because this dyanmic control retain the nonlinearity of the problem. There are additional control parameters that can be optimised, and Monte Carlo and genetic algorithm have both been appled to optimise these parameters in the hopes of obtaining the optimal conroller performance.  

## Files
**Para_Opt_MonteCarlo.m** - using monte carlo simulations to generate random samples within a range to find optimal control parameters  
**Para_Opt_GA.m** - using genetic algorithm to find the optimal control parameters 
