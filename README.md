# viraly
Virality concept demo

Examples:

```
python3 viraly.py "4.1,0.1,15,3,2,0.02,120,120,10276617,4,0.03"
python3 viraly.py "4.1,0.1,15,3,2,0.02,24 ,120,10276617,4,0.03"
```

Parameters and their meaning:

```
h     # average number of contacts per unit of time
p     # probability of transmission during a contact
T     # average duration of infection
L     # standard deviation of the normal distribution
h1    # average number of contacts per unit of time under contention
p1    # probability of transmission during a contact under contention
tint  # time with initial parameters (i.e., before contention)
tmax  # total time
M     # population size
N0    # initial number of infections
DR    # death rate
```
