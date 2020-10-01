#!/usr/bin/python3

import os
import sys
import scipy.stats
import numpy
import json
import matplotlib.pyplot as plt 
from distutils.util import strtobool
from collections import deque

# misc parameters
E_OK  = 0
E_ERR = 1

SEP = ';'
YLABEL_STR = 'Count'
STATS_STR  = 'transmissions, infections, recoveries, deaths'

# do we prefer model4 over model3?
PREFER_MOD4 = False

# shall we output the other models to the terminal?
OUTPUT_ALL = False

### functions ###

def print_usage ():
    basename = os.path.basename(sys.argv[0])
    print()
    print( 'Usage:\n\npython3 ' + basename + ' \"h,p,T,L,I,h2,p2,tint,tmax,M,N0,DR\"')
    print( 'python3 ' + basename + ' \"h,p,T,L,I,h2,p2,tint,tmax,M,N0,DR,progressive,ttime\"\n')

# model 1 - permanent infection, infinite population

def get_next_model1 ( current, h, p, M ):

    return min ( current*(1 + h*p), M )

# model 2 - permanent infection, finite population correction

def get_next_model2 ( current, h, p, M ):

    if current < 0:
        return 0

    if current == 0:
        return 0

    correction = ( 1 - current / M )
    return min ( current*(1 + h*p*correction), M)

# model 3 - temporary infection of fixed duration, finite population correction

def get_older_model3 ( time, history, T ):

    delta = time - T
    if delta < 0:
        return 0

    return history [ delta ]

# model 4 - temporary infection with gaussian duration of parameters T and L, finite population corrections

# note: this model leaves a residue of infections that do not disappear which is noticeable if L is high compared to T
# it could be fixed by complicating the code, but the physical situation does not make much sense
# 99% (so to say) of the normal should be to the right of t=0, or the physical model is not good

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

def get_older_model4 ( time, history, M, T, L ):

    aux_n    = 0
    aux_nc   = 0
    count    = 0

    if L == 0:
        L_effective = 1
    else:
        L_effective = L

    for j in range (0, time):
        aux_nc = history[j]
        aux_n = aux_nc*get_fraction( j + T, L_effective, time-1, time )
        count = count + aux_n

    return count

# common to models 3 and 4

def get_next_model34 ( current, h, p, time, nc_history, m, M, T, L, gaussian = False ):

    # we get the outgoing cases (recoveries, deaths) from the gaussian
    # outgoers are computed from the history of new cases either with
    # a batch recovery after T units of time (gaussian = False) or with
    # a recovery spread over moments controlled by normal distribution of
    # parameters T and L

    if gaussian:
        outgoing = get_older_model4 ( time, nc_history, M, T, L )
    else:
        outgoing = get_older_model3 ( time, nc_history, T )

    # the correction here is different becase current does not include outgoers...
    # we need to use the share of the population available for infection
    correction = max(( 1 - (M-m)/M ),0)
    # new cases
    nc = current*h*p*correction
    # Rt - attempt at estimating
    rt = h*p*T*correction

    return nc, outgoing, rt

# old helper function to model parameters evolution over time (supports only two stages)

def get_parameters_old ( h, p, h2, p2, t, tint, progressive = False, delta = 14 ):

    # check transtition time
    if progressive == False:
        ttime = 0
    else:
        ttime = delta

    # normal cases
    if t < tint:
        return h, p

    if t >= (tint + ttime):
        return h2, p2

    # if we are in the transition period
    if progressive == False:
        return h2, p2
    else:
        delta_t  = t-(tint+ttime)
        p_h2     = h2 + ((h2-h)/ttime)*delta_t
        p_p2     = p2 + ((p2-p)/ttime)*delta_t
        return p_h2, p_p2

# helper function to model parameters evolution over time

def get_parameters ( h, p, h2, p2, t, tint, progressive = False, delta = 14, h3 = 0, p3 = 0, tint2 = 0, delta2 = 0 ):

    if progressive == False:
        ttime  = 0
        ttime2 = 0
    else:
        ttime  = delta
        ttime2 = delta2

    # free phase
    if t < tint:
        #print ("debug:", t, h, p, h*p, 'free phase')
        return h, p

    # first transition
    if t in range(tint, tint + ttime):
        if progressive == False:
            return h2, p2
        else:
            delta_t  = t-(tint+ttime)
            p_h2     = h2 + ((h2-h)/ttime)*delta_t
            p_p2     = p2 + ((p2-p)/ttime)*delta_t
            #print ("debug:", t, p_h2, p_p2, p_h2*p_p2, 'first transition')
            return p_h2, p_p2

    # contention phase
    if tint2 > 0 and  t in range(tint + ttime, tint2 ):
        #print ("debug:", t, h2, p2, h2*p2, 'contention')
        return h2, p2

    # second transition
    if tint2 > 0 and t in range(tint2, tint2 + ttime2):
        if progressive == False:
            #print ("debug:", t, h3, p3, h3*p3, 'second free')
            return h3, p3
        else:
            delta_t  = t-(tint2+ttime2)
            p_h3     = h3 + ((h3-h2)/ttime2)*delta_t
            p_p3     = p3 + ((p3-p2)/ttime2)*delta_t
            #print ("debug:", t, p_h3, p_p3, p_h3*p_p3, 'second transition')
            return p_h3, p_p3

    # either second free phase or contention still
    if t >= tint2 + ttime2 :
        if tint2 > 0:
            #print ("debug:", t, h3, p3, h3*p3, 'second free')
            return h3, p3
        else:
            #print ("debug:", t, h2, p2, h2*p2, 'contention')
            return h2,p2

# plotting

# data is a list of lists of data
# labels is a list of labels

def plot_multiple ( data, labels, title_str, labely, legend_loc = "upper right" , block_execution = True):

    interval = len(data[0])
    x = numpy.linspace(0, interval, interval)

    plt.xkcd() # <3 <3 <3
    plt.figure() # necessary to make the plots separate
    plt.title( title_str)
    plt.xlabel('Time (days)')
    plt.ylabel(labely)

    for t in range (0, len(data)):
        plt_data  = data[t]
        plt_label = labels[t]
        if len(plt_data) > 0: 
            plt.plot(x, plt_data, label=plt_label )

    # fine tune legend location
    plt.legend(loc=legend_loc)

    # add a 10% head room for the y axis
    bottom, top = plt.ylim()
    plt.ylim((bottom, top*1.1))

    return plt

# print stuff to the terminal

def print_output ( t, x1, x2, x3_data, x4_data , prefer_x4 = False, output_all = False, output_silent = False ):

    if output_silent:
        return

    if prefer_x4:
        x_data = x4_data
    else:
        x_data = x3_data

    if output_all:
        print (t, SEP, x1, SEP, x2, SEP, x_data[0], SEP, x_data[1], SEP, x_data[2], SEP, x_data[3], SEP, x_data[4])
    else:
        print (t, SEP, x_data[0], SEP, x_data[1], SEP, x_data[2], SEP, x_data[3], SEP, x_data[4])

# main simulation function

def run_simulation ( h, p, T, L, I, h2, p2, tint, tmax, M, N0, DR, progressive, ttime, h3, p3, tint2, ttime2, silent, prefer_mod4 = PREFER_MOD4 ):

    # initial infections
    n1 = N0
    n2 = N0
    n3 = N0
    n4 = N0

    R0 = h*p*T

    # history of active numbers
    n1_history = [ N0 ]
    n2_history = [ N0 ]
    n3_history = [ N0 ]
    n4_history = [ N0 ]

    # fifos for cases in incubation
    incubator3 = deque([0]*(I-1))
    incubator4 = deque([0]*(I-1)) 

    # history of outgoing numbers
    o3_history = [ 0 ]
    o4_history = [ 0 ]

    # history of new cases
    nc3_history = [ N0 ]
    nc4_history = [ N0 ]

    # currently available population
    m3 = M - N0
    m4 = M - N0

    # history of available population
    m3_history = [ m3 ]
    m4_history = [ m4 ]

    n3_data = [ n3, N0, 0, M, R0 ]
    n4_data = [ n4, N0, 0, M, R0 ]

    # Rt history

    rt3_history = [ R0 ]
    rt4_history = [ R0 ]

    # stored parameters because h and p change over time
    sh = h
    sp = p

    # initial situation
    print_output (0, n1, n2, n3_data, n4_data, prefer_mod4, OUTPUT_ALL, silent )

    # we simulate tmax days, but the result contains the extra initial condition day at position 0
    for t in range (1, tmax + 1):
        # get new cases for the dummy models (new cases = active cases as there are no outgoers here)
        n1 = get_next_model1 (n1, h, p, M)
        n2 = get_next_model2 (n2, h, p, M)
        # get new cases, outgoing and rt3 for the two models that matter
        nc3i, o3, rt3 = get_next_model34 (n3, h, p, t, nc3_history, m3, M, T, L, False)
        nc4i, o4, rt4 = get_next_model34 (n4, h, p, t, nc4_history, m4, M, T, L, True)
        # update simulation parameters over time
        h, p = get_parameters( h,p, h2, p2, t, tint, progressive, ttime, h3, p3, tint2, ttime2)

        # but nc3i and nc4i go for incubation still and we need to fetch the ones that are ready to infect
        # in the SIER model incubation cases are called "Exposed" - they are infected but not infectious
        # note: we need to append before popping to support the case where the incubation time is 1 
        # which recovers the tried and tested behaviour we had before introducing this parameter
        incubator3.appendleft(nc3i)
        incubator4.appendleft(nc4i)
        nc3 = incubator3.pop()
        nc4 = incubator4.pop()

        # new current - it sometimes goes negative by a very small value
        n3 = max(n3 + nc3 - o3,0)
        n4 = max(n4 + nc4 - o4,0)

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

        # neither the outgoing nor the exposed (i.e. in incubation) are available targets for new infections
        # but the infected are still causing new infections
        # note: we remove the cases for the susceptibles pool as soon as they are exposed (nc3i instead of nc3, etc)
        m3 = max(m3 - nc3i, 0)
        m4 = max(m4 - nc4i, 0)

        m3_history.append(m3)
        m4_history.append(m4)

        rt3_history.append(rt3)
        rt4_history.append(rt4)

        n3_data = [ n3, nc3, o3, m3, rt3 ]
        n4_data = [ n4, nc4, o4, m4, rt4 ]
        print_output (t, n1, n2, n3_data, n4_data, prefer_mod4, OUTPUT_ALL, silent )

    # deaths vs recoveries

    d3_history = numpy.array(o3_history) * DR
    r3_history = numpy.array(o3_history) * (1-DR)
    d4_history = numpy.array(o4_history) * DR
    r4_history = numpy.array(o4_history) * (1-DR)

    # choose which epidemic model in use from here on

    if prefer_mod4:
        n_final    = n4
        m_final    = m4
        n_history  = n4_history
        nc_history = nc4_history
        d_history  = d4_history
        r_history  = r4_history
        o_history  = o4_history
        m_history  = m4_history
        rt_history = rt4_history
    else:
        n_final    = n3
        m_final    = m3
        n_history  = n3_history
        nc_history = nc3_history
        d_history  = d3_history
        r_history  = r3_history
        o_history  = o3_history
        m_history  = m3_history
        rt_history = rt3_history

    # calculate and print some statistics

    t_transmissions = numpy.array(nc_history).sum()
    t_infections    = t_transmissions + N0
    t_inactivations = numpy.array(d_history).sum()
    t_recoveries    = numpy.array(r_history).sum()
    t_removals      = numpy.array(o_history).sum()

    # prepare some acumulated data

    j = 1
    na_history = []
    da_history = []
    ra_history = []

    for value in n_history:
        na = numpy.array(nc_history[0:j]).sum()
        da = numpy.array(d_history[0:j]).sum()
        ra = numpy.array(r_history[0:j]).sum()
        na_history.append(na)
        da_history.append(da)
        ra_history.append(ra)
        j=j+1

    # plot time

    if not silent:
        print ('Maximum value')
        print (numpy.argmax(n_history), ' ' , numpy.amax(n_history))

        print('Totals:')
        print (STATS_STR)
        print ( t_transmissions, t_infections, t_recoveries, t_inactivations )

        # technical string that labels the plot with the simulation parameters

        tech_str = 'h={h}, p={p}, T={T}, L={L}, h2={h2}, p2={p2}, tint={tint}, tmax={tmax}, M={M}, N0={N0}, DR={DR} progressive={progressive} ttime={ttime}'.format(h=sh, p=sp, T=T, L=L, h2=h2,p2=p2, tint=tint, tmax=tmax, M=M, N0=N0, DR=DR, progressive=progressive, ttime=ttime)

        # produce a complete plot for the chosen epidemic model

        mydata   = [ n_history,      nc_history,   r_history,    d_history ]
        mylabels = [ 'Active cases', 'New Cases',  'Recoveries', 'Deaths'  ]

        plt1 = plot_multiple( mydata, mylabels, tech_str, YLABEL_STR, "upper right" )
        # plot acumulated cases and acumulated deaths

        mydata   = [ na_history,          da_history ]
        mylabels = [ 'Acumulated cases', 'Acumulated deaths' ]

        plt2 = plot_multiple( mydata, mylabels, tech_str, YLABEL_STR, "upper left" )

        # typical SIR plot with Susceptible, Infected and Removed (Recovered or Dead)

        # population history: should be constant, can be added to the plot just to check consistency
        po_history = numpy.array(m_history) +  numpy.array(n_history) + numpy.array(ra_history) + numpy.array(da_history)

        mydata   = [ m_history,     n_history,  ra_history,  da_history ]
        mylabels = [ 'Susceptible', 'Infected', 'Recovered', 'Dead'     ]

        plt3 = plot_multiple( mydata, mylabels, tech_str, YLABEL_STR, "upper left" )

        # compare epidemic model with simple exponential and logisic models

        mydata   = [ n1_history,      n2_history,   n3_history, n4_history  ]
        mylabels = [ 'Exponential',   'Logistic',  'Epidemic',  'Epidemic2' ]

        plt4 = plot_multiple( mydata, mylabels, tech_str, YLABEL_STR, "upper left" )

        # plot Rt

        mydata   = [ rt_history ]
        mylabels = [ 'R(t)'     ]

        plt4 = plot_multiple( mydata, mylabels, tech_str, YLABEL_STR, "upper left" )

        plt1.show(block = True)
    else:
        # the list cast is only to uniformized because some of the elements were converted to numpy arrays
        dataset = [ n_history, nc_history, list(r_history), list(d_history), m_history, n_history, ra_history, da_history, rt_history, na_history ]

        return dataset

# optimized version only to be used by the web interface:
# runs model 4 and is silent

def run_simulation_web ( h, p, T, L, I, h2, p2, tint, tmax, M, N0, DR, progressive, ttime, h3, p3, tint2, ttime2, silent = True, prefer_mod4 = PREFER_MOD4 ):

    n4 = N0
    R0 = h*p*T

    # history of active numbers
    n4_history = [ N0 ]

    # fifos for cases in incubation
    incubator4 = deque([0]*(I-1))

    # history of outgoing numbers
    o4_history = [ 0 ]

    # history of new cases
    nc4_history = [ N0 ]

    # currently available population
    m4 = M - N0

    # history of available population
    m4_history = [ m4 ]

    n4_data = [ n4, N0, 0, M, R0 ]

    # Rt history
    rt4_history = [ R0 ]

    # stored parameters because h and p change over time
    sh = h
    sp = p

    # we simulate tmax days, but the result contains the extra initial condition day at position 0
    for t in range (1, tmax + 1):

        # get new cases, outgoing and rt
        nc4i, o4, rt4 = get_next_model34 (n4, h, p, t, nc4_history, m4, M, T, L, prefer_mod4)
        # update simulation parameters over time
        h, p = get_parameters( h,p, h2, p2, t, tint, progressive, ttime, h3, p3, tint2, ttime2)

        # but nc3i and nc4i go for incubation still and we need to fetch the ones that are ready to infect
        # in the SIER model incubation cases are called "Exposed" - they are infected but not infectious
        # note: we need to append before popping to support the case where the incubation time is 1
        # which recovers the tried and tested behaviour we had before introducing this parameter
        incubator4.appendleft(nc4i)
        nc4 = incubator4.pop()

        # new current - it sometimes goes negative by a very small value
        n4 = max(n4 + nc4 - o4,0)

        # new cases that appeared at time t
        nc4_history.append(nc4)

        # cases that went out at time t
        o4_history.append(o4)

        # number of active cases at time t
        n4_history.append(n4)

        # neither the outgoing nor the exposed (i.e. in incubation) are available targets for new infections
        # but the infected are still causing new infections
        # note: we remove the cases for the susceptibles pool as soon as they are exposed (nc3i instead of nc3, etc)
        m4 = max(m4 - nc4i, 0)
        m4_history.append(m4)

        rt4_history.append(rt4)

        n4_data = [ n4, nc4, o4, m4, rt4 ]

    # deaths vs recoveries

    d4_history = numpy.array(o4_history) * DR
    r4_history = numpy.array(o4_history) * (1-DR)

    n_final    = n4
    m_final    = m4
    n_history  = n4_history
    nc_history = nc4_history
    d_history  = d4_history
    r_history  = r4_history
    o_history  = o4_history
    m_history  = m4_history
    rt_history = rt4_history

    # calculate and print some statistics

    t_transmissions = numpy.array(nc_history).sum()
    t_infections    = t_transmissions + N0
    t_inactivations = numpy.array(d_history).sum()
    t_recoveries    = numpy.array(r_history).sum()
    t_removals      = numpy.array(o_history).sum()

    # prepare some acumulated data

    j = 1
    na_history = []
    da_history = []
    ra_history = []

    for value in n_history:
        na = numpy.array(nc_history[0:j]).sum()
        da = numpy.array(d_history[0:j]).sum()
        ra = numpy.array(r_history[0:j]).sum()
        na_history.append(na)
        da_history.append(da)
        ra_history.append(ra)
        j=j+1

    # the list cast is only to uniformized because some of the elements were converted to numpy arrays
    dataset = [ n_history, nc_history, list(r_history), list(d_history), m_history, n_history, ra_history, da_history, rt_history, na_history ]

    return dataset

def main():

    # the silent mode is for integration with external tools, it only exports dataset
    silent = False

    # parse input

    if len(sys.argv) < 2:
        print_usage()
        exit(E_OK)

    if len(sys.argv) > 2:
        silent=True

    # simulation parameters

    # We accept all the CLI arguments in a single comma separated string. Check the documentation for examples.

    myparams_str  = sys.argv[1]
    myparams_list = myparams_str.split(',')

    if len(myparams_list) < 12:
        print_usage()
        exit(E_OK)

    h     = float(myparams_list[0])  # average number of contacts per unit of time
    p     = float(myparams_list[1])  # probability of transmission during a contact
    T     = int  (myparams_list[2])  # average duration of infection
    L     = int  (myparams_list[3])  # standard deviation of the normal distribution
    I     = int  (myparams_list[4])  # incubation time
    h2    = float(myparams_list[5])  # average number of contacts per unit of time under contention
    p2    = float(myparams_list[6])  # probability of transmission during a contact under contention
    tint  = int  (myparams_list[7])  # time with initial parameters (i.e., before contention)
    tmax  = int  (myparams_list[8])  # total time
    M     = float(myparams_list[9])  # population size
    N0    = float(myparams_list[10]) # initial number of infections
    DR    = float(myparams_list[11]) # death rate

    if len(myparams_list) > 12:
        progressive = strtobool(myparams_list[12])
    else:
        progressive = False

    if len(myparams_list) > 13:
        ttime = int(myparams_list[13])
    else:
        ttime = 0

    # bonus: stage 3
    if len(myparams_list) > 17:
        h3     = float(myparams_list[14])  # average number of contacts per unit of time after contention
        p3     = float(myparams_list[15])  # probability of transmission during a contact after contention
        tint2  = int(myparams_list[16])    # time at which we start the second transition
        ttime2 = int(myparams_list[17])    # x2 -> x3 parameters transition duration
    else:
        h3     = 0
        p3     = 0
        tint2  = 0
        ttime2 = 0

    if tint > tmax:
        print('tint2 must be smaller than', tmax)
        exit(E_ERR)

    if tint2 > 0 and tint2 > tmax:
        print('tint2 must be smaller than', tmax)
        exit(E_ERR)

    if tint2 > 0 and tint2 < tint + ttime:
        print('tint2 must be greater than', tint, '+', ttime)
        exit(E_ERR)

    # this is another bonus for external tools integration which does not break the historical CLI usage:
    #   we allow them model selection to be done as a function of the value of L:
    #     if L = 0 -> model 3
    #     if L > 0 -> model 4
    #
    # NOTE: model4 is much much slower than model3

    # simulation
    if len(myparams_list) > 18:
        if L == 0:
            prefer_mod4 = False
        else:
            prefer_mod4 = True
    else:
        prefer_mod4 = PREFER_MOD4

    dataset = run_simulation ( h, p, T, L, I, h2, p2, tint, tmax, M, N0, DR, progressive, ttime, h3, p3, tint2, ttime2, silent, prefer_mod4 )
    #print(dataset)
 
### Main block ###

if __name__ == "__main__":
    main()

