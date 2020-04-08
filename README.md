# Viraly - Virality concepts demo

**Introduction**

This is a simple command line program that simulate an epidemic over user given parameters under 4 different models:

* permanent infection, infinite population: exponential growth
* permanent infection, finite population: logistic growth
* temporary infection with duration T, finite population: epidemic curve
* temporary infection with guassian duration of avg T and stdev L, finite population: epidemic curve

**Parameters and their meaning:**
```
h     # average number of contacts per unit of time
p     # probability of transmission during a contact
T     # average duration of infections
L     # standard deviation of the normal distribution of the infection duration
h1    # average number of contacts per unit of time under contention
p1    # probability of transmission during a contact under contention
tint  # simulation time with initial parameters (i.e., before contention)
tmax  # total time
M     # population size
N0    # initial number of infections
DR    # death rate
```
The parameters must be given via command line in the order listed above as quoted comma-separated list.

**Examples:**
```
python3 viraly.py "4,0.1145,15,3,2,0.02,120,120,10276617,4,0.03"
python3 viraly.py "4,0.1145,15,3,2,0.02,24 ,120,10276617,4,0.03"
```

The first example simulates a free epidemic for 120 days, whereas the second example simulates an epidemic for 120 days with sudden change of h and p (contention) at time t=24.

**Outputs:**

* plot with active cases, new cases, recoveries an deaths
* plot with acumulated cases and acumulated deaths
* plot with the usual SIR variables: Susceptible, Infected and Removed (Recovered or Dead)
* plot with comparison of models: exponential, logistic and epidemic (with fixed recovery time) and epidemic2 (with gaussian recovery time)

**Example outputs:**

Output for model 4 with a sudden parameter change (contention) at t=24:
![Output for model4 with a parameter change shock at t=24](https://github.com/ghomem/viraly/blob/master/images/example_t24_shock.png)

Same as above but as a free epidemic (no parameter change):
![Same as above but as a free epidemic (no parameter change)[https://github.com/ghomem/viraly/blob/master/images/example_no_shock.png)

SIR plot for the case above:
![SIR plot)[https://github.com/ghomem/viraly/blob/master/images/example_no_SIR.png)

Comparison of models for the case above:
![Comparison)[https://github.com/ghomem/viraly/blob/master/images/example_no_comp.png)

**Disclaimer:**

This is an experiment related to the math of virality (whether it means marketing, ideas, content or something else) and should not be used for decisions related to real world public health situations.
