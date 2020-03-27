#!/usr/bin/python3

import os
import sys
import scipy.stats
import numpy
import matplotlib.pyplot as plt 

# misc parameters
E_OK  = 0
E_ERR = 1

SEP = ';'
YLABEL_STR = 'Aware / Recovered / Inactive'
STATS_STR  = 'Total transmissions, Total infections, recoveries, inactivations'

#STATS_STR = 'Total cases, recoveries, deaths'
#YLABEL    = 'Active / Recovered / Dead'

### functions ###

def print_usage ():
    basename = os.path.basename(sys.argv[0])
    print()
    print( 'Usage:\n\npython3 ' + basename + ' \"h,p,T,L,h1,p1,tint,tmax,M,N0,DR\"\n')

# model 1 - permanent infection, infinite population
def get_next_model1 ( current, h, p ):

    return min ( current*(1 + h*p), M )

# model 2 - permanent infection, finite population correction
def get_next_model2 ( current, h, p ):

    if current < 0:
        return 0

    if current == 0:
        return 0

    correction = ( 1 - current / M )
    return min ( current*(1 + h*p*correction), M)

# model 3 - temporary infection of fixed duration, finite population correction

def get_older_model3 ( time, history ):

    delta = time - T
    if delta < 0:
        return 0

    return history [ delta ]

def get_next_model3 ( current, h, p, time, history, m ):

    outgoing = get_older_model3 ( time, history )
    # the correction here is different becase current does not include outgoers...
    # we need to use the share of the population available for infection
    correction = max(( 1 - (M-m)/M ),0)
    # new cases
    nc = current*h*p*correction
    # new current - it sometimes goes negative by a very small value
    current_updated = max(current + nc - outgoing,0)

    return current_updated, nc, outgoing

# model 4 - temporary infection with gaussian duration of parameters T and L, finite population corrections

def get_fraction ( center, stdev, t1, t2 ):

    n1 = scipy.stats.norm.cdf( t1, center, stdev )
    n2 = scipy.stats.norm.cdf( t2, center, stdev )

    return n2 - n1

def test_fraction ():

    # test "standard normal" 
    # which denotes the normal distribution with zero mean and unit variance

    print ( get_fraction ( 0, 1, 0 , 1) ) # centered in 0, stdev 1, interval [0,1]
    print ( get_fraction ( 0, 1, -1, 0) ) # centered in 0, stdev 1, interval [-1,0]
    print ( get_fraction ( 0, 1, -1, 1) ) # centered in 0, stdev 1, interval [-1,1]

    # comparison to this table 
    # https://en.wikipedia.org/wiki/Standard_normal_table#Cumulative_from_mean_(0_to_Z)

    print ( get_fraction ( 0, 1, 0 , 0.09) ) # centered in 0, stdev 1, interval [0,0.09], result 0.03586
    print ( get_fraction ( 0, 1, 0 , 0.19) ) # centered in 0, stdev 1, interval [0,0.19], result 0.07535
    print ( get_fraction ( 0, 1, 0 , 0.29) ) # centered in 0, stdev 1, interval [0,0.29], result 0.11409

def get_older_model4 ( time, history ):

    aux_t    = 0
    aux_n    = 0
    count    = 0

    for n in history:
        aux_n = n*get_fraction(T + aux_t, L, time-1, time )
        #print('debug', time, aux_t, n, aux_n)
        count = count + aux_n
        aux_t = aux_t + 1

    #print('debug out 4', count)
    return count

def get_next_model4 ( current, h, p, time, history, m ):

    # we get the outgoing cases (recoveries) from the gaussian 
    outgoing = get_older_model4 ( time, history )
    return min (max (0, get_next_model2 ( current, h, p ) - outgoing ), m), outgoing, 0 # FIXME all broken

# plotting

def plot_data ( data1, data2, data3, data4, label_str, labely, allplots ):

    interval = len(data1)
    x = numpy.linspace(0, interval, interval)
    plt.xkcd() # <3 <3 <3
    plt.title( label_str)
    plt.xlabel('Time (days)')
    plt.ylabel(labely)
    plt.plot(x,data1)

    if allplots:
        if len(data2) > 0:
            plt.plot(x,data2)
        if len(data3) > 0:
            plt.plot(x,data3)
        if len(data4) > 0:
            plt.plot(x,data4)

    plt.show()
    plt.clf()


### main ###

# parse input

if len(sys.argv) < 2:
    print_usage()
    exit(E_OK)

### simulation parameters ###

# we accept CLI arguments for h,p,T,L,h1,p1,tint,tmax,M,N0,DR in a single string
#
# example: 
# python3 viraly.py "4.574,0.1,15,3,1,0.02,24,60,10000000,4,0.03"

myparams_str  = sys.argv[1]
myparams_list = myparams_str.split(',')

if len(myparams_list) < 11:
    print_usage()
    exit(E_OK)

h     = float(myparams_list[0])  # average number of contacts per unit of time
p     = float(myparams_list[1])  # probability of transmission during a contact
T     = int  (myparams_list[2])  # average duration of infection
L     = int  (myparams_list[3])  # standard deviation of the normal distribution
h1    = float(myparams_list[4])  # average number of contacts per unit of time under contention
p1    = float(myparams_list[5])  # probability of transmission during a contact under contention
tint  = int  (myparams_list[6])  # time with initial parameters (i.e., before contention)
tmax  = int  (myparams_list[7])  # total time
M     = float(myparams_list[8])  # population size
N0    = float(myparams_list[9])  # initial number of infections
DR    = float(myparams_list[10]) # death rate

# simulation

n1 = N0
n2 = N0
n3 = N0
n4 = N0

# history of active numbers
n1_history = [ N0 ]
n2_history = [ N0 ]
n3_history = [ N0 ]
n4_history = [ N0 ]

# history of outgoing numbers
o3_history = [ 4 ]
o4_history = [ 4 ]

# history of new cases
nc3_history = [ 4 ]
nc4_history = [ 4 ]

# currently available population

m3 = M
m4 = M

# initial situation
print (0, SEP, n1, SEP, n2, SEP, n3)

for t in range (1, tint):
    n1 = get_next_model1 (n1, h, p) 
    n2 = get_next_model2 (n2, h, p) 
    n3, nc3, o3 = get_next_model3 (n3, h, p, t, nc3_history, m3)
    n4, nc4, o4 = get_next_model4 (n4, h, p, t, nc4_history, m4)

    # new cases that appeared at time t
    nc3_history.append(nc3)
    nc4_history.append(nc4)
   
    # cases that went out at time t 
    o3_history.append(o3)
    o4_history.append(o4)

    # number of active cases at time t
    n1_history.append(n1)
    n2_history.append(n2)
    n3_history.append(n3)
    n4_history.append(n4)

    # neither the outgoing nor the infected are available targets for new infections
    # but infected are still infecting causing new infections
    m3 = max(m3 - nc3 - o3,0)
    m4 = m4 - nc4 - o4

    #print (t, SEP, n1, SEP, n2, SEP, n3, SEP, n4, SEP, o3, SEP, nc3, SEP, m3) #FIXME remove o3
    print (t, SEP, n1, SEP, n2, SEP, n3, SEP, o3, SEP, nc3, SEP, m3) #FIXME remove o3

# change in parameters
for t in range (tint, tmax):
    n1 = get_next_model1 (n1, h1, p1) 
    n2 = get_next_model2 (n2, h1, p1) 
    n3, o3 = get_next_model3 (n3, h1, p1, t, nc3_history, m3)
    n4, o4 = get_next_model4 (n4, h1, p1, t, nc4_history, m4)

    nc3 = max( n3 - n3_history[-1], 0)
    nc4 = max( n4 - n4_history[-1], 0)  
    nc3_history.append(nc3)
    nc4_history.append(nc4)
    
    o3_history.append(o3)
    o4_history.append(o4)

    n1_history.append(n1)
    n2_history.append(n2)
    n3_history.append(n3)
    n4_history.append(n4)

    m3 = m3 - nc3 - o3
    m4 = m4 - nc4 - o4

    #print (t, SEP, n1, SEP, n2, SEP, n3, SEP, n4, SEP, o3, SEP, nc3, SEP, m3) #FIXME remove o3
    print (t, SEP, n1, SEP, n2, SEP, n3, SEP, o3, SEP, nc3, SEP, m3) #FIXME remove o3

# deaths vs recoveries

d3_history = numpy.array(o3_history) * DR
r3_history = numpy.array(o3_history) * (1-DR)
d4_history = numpy.array(o4_history) * DR
r4_history = numpy.array(o4_history) * (1-DR)

# choose which model in use from here on

n_history  = n3_history
nc_history = nc3_history
d_history  = d3_history
r_history  = r3_history

# calculate some statistics

print ('Maximum value')
print (numpy.argmax(n_history), ' ' , numpy.amax(n_history))

print (STATS_STR)
t_transmissions = numpy.array(nc_history).sum()
t_infections   =  t_transmissions + N0
print ( t_transmissions, t_infections, numpy.array(r_history).sum(), numpy.array(d_history).sum() )

# produce a plot

# technical string that labels the plot with the simulation parameters
tech_str = 'h={h}, p={p}, T={T}, L={L}, h1={h1}, p1={p1}, tint={tint}, tmax={tmax}, M={M}, N0={N0}, DR={DR}'.format(h=h, p=p, T=T, L=L, h1=h1,p1=p1, tint=tint, tmax=tmax, M=M, N0=N0, DR=DR)

plot_data( n_history, nc_history, r_history, d_history, tech_str, YLABEL_STR, True )

