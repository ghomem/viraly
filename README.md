# Viraly - Virality concepts demo

**Introduction**

This is a simple command line program that simulates an epidemic in 4 different models:

1. permanent infection, infinite population: exponential growth
2. permanent infection, finite population: logistic growth
3. temporary infection with duration T, finite population: epidemic curve
4. temporary infection with guassian duration of avg T and stdev L, finite population: epidemic curve

An epidemic is the propagation of something (a disease, a phrase, a brand, an idea, ...) over a population by means of interactions between elements. A technical summary document can be found [here](https://github.com/ghomem/viraly/blob/master/doc/viral-summary-github.pdf). There is a bokeh based web frontend on the /web folder an instance of which is live [here](https://lo.gic.li/viral).

**Requirements**

Scipy, Numpy and Matplotlib

**Usage**

```
python3 viraly.py "h,p,T,L,I,h2,p2,tint,tmax,M,N0,DR"
python3 viraly.py "h,p,T,L,I,h2,p2,tint,tmax,M,N0,DR,progressive,ttime"
```

**Parameters and their meaning**
```
h           # average number of contacts per unit of time
p           # probability of transmission during a contact
T           # average duration of infections
L           # standard deviation of the normal distribution of the infection duration
I           # incubation time
h2          # average number of contacts per unit of time under contention
p2          # probability of transmission during a contact under contention
tint        # simulation time with initial parameters (i.e., before contention)
tmax        # total simulation time
M           # population size
N0          # initial number of infections
DR          # death rate
progressive # [optional] whether or not the change of parameters at time tint should be progressive
ttime       # [optional] the parameters transition time (if progressive == True)
```
The parameters must be given via command line in the order listed above as quoted comma-separated list.

Notes:
* h and p are the parameters that drive propagation
* they  are presented as independent for physical intuition purposes but only the product hp matters in practice
* the simulation includes two phases:  **\[ 0, tint \[** with parameters (h,p) and **\[ tint, tmax \[**  with parameters (h2,p2)
* if tint == tmax the simulation reduces to a single phase with parameters (h,p)
* if progressive == True the transition between phases is done using a linear variation (h,p) -> (h2,p2)
* the Basic Reproduction Number (R0) is given by hpT
* for model 3 the infections remain constant if hp = 1/T and decrease to zero if hp < 1/T
* for model 4 hp needs to be slightly lower for the situations above to occur

**Examples**
```
python3 viraly.py "4,0.1145,15,3,1,2,0.02,120,120,10276617,4,0.03"
python3 viraly.py "4,0.1145,15,3,1,2,0.02,24 ,120,10276617,4,0.03"
python3 viraly.py "4,0.1145,15,3,1,2,0.02,24 ,120,10276617,4,0.03,True,7"
```

The first example simulates a free epidemic for 120 days, whereas the second example simulates an epidemic for 120 days with sudden change of h and p (contention) at time t=24. The third example is equal to the second but with a linear change of parameters over the course of 7 days.

**Outputs**

* plot with active cases, new cases, recoveries an deaths
* plot with acumulated cases and acumulated deaths
* plot with the usual SIR variables: Susceptible, Infected and Removed (Recovered or Dead)
* plot with comparison of models: exponential, logistic and epidemic (with fixed recovery time) and epidemic2 (with gaussian recovery time)
* plot with evolution of the effective reproduction number R(t) over time
* console output for the preferred model with Active Cases, New Cases, Removals, Susceptibles and R(t), plus misc stats (peak, totals, ...)

**Configuration**

The boolean global variable PREFER_MOD4 controls whether or not model4 is the preferred model for console output and plots. If PREFER_MOD4 is False the preferred model is model3.

**Example outputs**

Example 1: output for model 4 with a sudden parameter change (contention) at t=24 such that h<sub>2</sub>p<sub>2</sub>T < 1:
![Output for model4 with a parameter change shock at t=24 hpT < 1](https://github.com/ghomem/viraly/blob/master/images/example_t24_shock.png)

Example 2: output for model 4 with a sudden parameter change (contention) at t=24 such that h<sub>2</sub>p<sub>2</sub>T > 1:
![Output for model4 with a parameter change shock at t=24 htT > 1](https://github.com/ghomem/viraly/blob/master/images/example_t24_shock_larger_R.png)

Example 3: same as above but as a free epidemic (no parameter change):
![Same as above but as a free epidemic (no parameter change)](https://github.com/ghomem/viraly/blob/master/images/example_no_shock.png)

Example 4: SIR plot for the case above:
![SIR plot](https://github.com/ghomem/viraly/blob/master/images/example_no_shock_SIR.png)

Example 5: comparison of models for the case above:
![Comparison](https://github.com/ghomem/viraly/blob/master/images/example_no_shock_comp.png)

Example 6: output for model 3 with critical choice of parameters hp = 1/T so that the epidemic is in the limit of propagation:
![Critical](https://github.com/ghomem/viraly/blob/master/images/example_no_shock_prop_stall_model3.png)

Example 7: same as above for model 4 where due to the gaussian recovery hp needs to be slightly lower:
![Critical](https://github.com/ghomem/viraly/blob/master/images/example_no_shock_prop_stall_model4.png)

Example 8: output of early erratication case due to sub critical hp for model 4:
![Critical](https://github.com/ghomem/viraly/blob/master/images/example_no_shock_erradication_model4.png)

**Bonus section**

There are extra parameters not documented here which can be checkend on the source code. With those parameters a 3 stage simulation can be performed and the script can be integrated with an external application as demonstrated in this [web frontend](https://lo.gic.li/viral).

**Disclaimer**

This is an experiment related to the math of virality and should not be used for decisions related to real world public health situations.
