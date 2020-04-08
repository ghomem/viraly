# Viraly - Virality concepts demo

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

The first example simulates a free epidemic for 120 days, whereas the second example simulates an epidemic for 120 days with sudden change of h and p (contenton) at time t=24.

**Outputs:**

* plot with active cases, new cases, recoveries an deaths
* plot with acumulated cases and acumulated deaths
* plot with comparison of models: exponential, logistic and epidemic (with fixed recovery time) and epidemic2 (with gaussian recovery time)

**Disclaimer:**

This is an experiment related to the math of virality (whether it means marketing, ideas, content or something else) and should not be used for decisions related to real world public health situations.
