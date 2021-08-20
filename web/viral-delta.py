''' 
Present an interactive epidemic simulation

'''
import subprocess
import json
import numpy as np

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider, TextInput, BoxAnnotation, HoverTool, Button, Spacer, Div
from bokeh.plotting import figure

# imports from separate  file
from viraly import *

### Configuration

# Risk plot, for testing
ENABLE_RISK = False

# Geometry and visuals
PLOT_TOOLS    ='save,reset,pan,wheel_zoom,box_zoom'
PLOT_HEIGHT   = 300
PLOT_WIDTH    = 500
TEXT_WIDTH    = 300
LMARGIN_WIDTH = 20

PLOT_LINE_WIDTH  = 3
PLOT_LINE_ALPHA  = 0.9

PLOT_LINE_ALPHA_DIFF = 0.7

PLOT_X_LABEL  = 'Days'
PLOT_Y_LABEL  = 'Count'
PLOT_Y_LABEL2 = 'Value'

PLOT_LINE_ACTIVE_COLOR    = 'blue'
PLOT_LINE_NEW_COLOR       = 'red'
PLOT_LINE_RECOVERED_COLOR = 'green'
PLOT_LINE_DEAD_COLOR      = 'black'

# Bureaucracy :-)
APP_DIR    = '/home/gustavo/viral'
CMD_DATA   = APP_DIR + '/viraly.py'
CMD_PYTHON = '/usr/bin/python3'

# Population
POP_MIN   = 1
POP_MAX   = 330
POP_START = 10.2 
POP_STEP  = 0.5

# Initial infections
IIF_MIN   = 0
IIF_MAX   = 10000
IIF_START = 7500
IIF_STEP  = 50

# Infectious period
T_MIN   = 1
T_MAX   = 15
T_START = 6

# and its Standard Deviation
T_STDEV_MIN   = 0
T_STDEV_MAX   = 3
T_STDEV_START = 0

# Latent period
# Not the same as incubation period, and can be shorter:
# https://en.wikipedia.org/wiki/Latent_period_(epidemiology)
L_MIN   = 1
L_MAX   = 20
L_START = 4

# Simulation time
DAYS = 365

# Propagation rate parameters
#
# Note:
#   BETA = hp
#   R0   = hpT
#
BETA_MIN  =  0

H1_MIN   = 0
H1_MAX   = 100
H1_START = 8
H1_STEP  = 0.1

P1_MIN   = 0
P1_MAX   = 100
P1_START = 0.125 * 100
P1_STEP  = 0.01

# uncertainty %
P_DELTA_MIN   = 0
P_DELTA_MAX   = 30
P_DELTA_START = 20
P_DELTA_STEP  = 0.5

DRATE_MIN   = 0.05
DRATE_MAX   = 10
DRATE_START = 0.15
DRATE_STEP  = 0.05

IM_MIN   = 0
IM_MAX   = 100
IM_START = 50
IM_STEP  = 0.5

# for the incidence plot
INCIDENCE_PERIOD = 14

# labels and strings
PAGE_TITLE  ='Simple epidemic simulator'
PLOT_TITLE  ='Active'
PLOT2_TITLE ='New, Recovered'
PLOT3_TITLE ='Rt estimation'
PLOT4_TITLE ='Immunity'
PLOT5_TITLE ='Dead'
PLOT6_TITLE ='Accumulated cases / recoveries'
PLOT7_TITLE ='Accumulated deaths'
PLOT8_TITLE =str(INCIDENCE_PERIOD) + ' day incidence per 100m'
PLOT9_TITLE ='Prevalence'

T_LABEL       = 'Infectious period'
T_STDEV_LABEL = 'NOT IN USE Infectious Period Standard Deviation'
L_LABEL       = 'Latent Period'
POP_LABEL     = 'Population (Millions)'
IIF_LABEL     = 'Initial number of infections'
H1_LABEL      = 'Organic contacts per day'
P1_LABEL      = 'Probability of transmission (%)'
P_DELTA_LABEL = 'Uncertainty on the prob/ of transmission (%)'
DRATE_LABEL   = 'Death rate (%)'
IM_LABEL      = 'Pre immunized (%)'

TEXT_INTRO    = 'Use the mouse for initial selection and cursors for fine tuning:'
TEXT_SUMMARY  = 'Stats:'
TEXT_NOTES    ='<b>Notes:</b><br/>\
              &bull; &beta; = hp<br/>\
              &bull; R<sub>0</sub> = hpT<br/>\
              &bull; Technical info at <a href="https://github.com/ghomem/viraly">github.com/ghomem/viraly</a>'
### End of configuration

### Functions

# the function that we are plotting
def get_data(x, pop, n0, period, period_stdev, latent, d1, d2, tr1, tr2, b1, b2,b3, tmax, dr, prog_change, IM = 0 ):

    h  = 1
    p  = float (b1 / 100) # input is multiplied by 100 for precision on the sliders
    h2 = 1
    p2 = float (b2 / 100) # input is multiplied by 100 for precision on the sliders
    h3 = 1
    p3 = float (b3 / 100) # input is multiplied by 100 for precision on the sliders
    T  = period
    I  = latent
    N0 = n0
    DR = float(dr/100)   # input is in percentage
    M  =  pop*1000000    # input is in millions
    L  = period_stdev
    I0 = round(M*IM/100) # input is in percentage

    tint  = d1
    tint2 = d1 + d2
    progressive = prog_change # bool
    ttime  = tr1
    ttime2 = tr2

    # decide which model to use based on the value of L
    if L == 0:
        prefer_mod4 = False
    else:
        prefer_mod4 = True

    # prepare debug friendly string for CLI troubleshoot
    str_params = '{h},{p},{T},{L},{I},{h2},{p2},{tint},{tmax},{M},{N0},{DR},{progressive},{ttime},{h3},{p3},{tint2},{ttime2},{prefer_mod4},{I0}'.format(h=h, p=p, T=T, L=L, I=I, h2=h2,p2=p2,       \
                                                                                                                                                  tint=tint, tmax=tmax, M=M, N0=N0, DR=DR,           \
                                                                                                                                                  progressive=progressive, ttime=ttime, h3=h3, p3=p3,\
                                                                                                                                                  tint2=tint2, ttime2=ttime2, prefer_mod4=prefer_mod4, I0=I0)
    print(str_params)

    # this function is included from viraly.py
    silent = True
    top_level = run_simulation_web ( h, p, T, L, I, h2, p2, tint, tmax, M, N0, DR, progressive, ttime, h3, p3, tint2, ttime2, silent, prefer_mod4, I0 )

    # dataset from viraly.py: [ n_history, nc_history, list(r_history), list(d_history), m_history, n_history, ra_history, da_history, rt_history, na_history, i_history ]
    # we chop the first element because it is the initial condition (ex: new cases don't make sense there, especially on a second wave simulation )
    n_history  = top_level[0][1:]
    nc_history = top_level[1][1:]
    r_history  = top_level[2][1:]
    d_history  = top_level[3][1:]
    ra_history = top_level[6][1:]
    da_history = top_level[7][1:]
    rt_history = top_level[8][1:]
    na_history = top_level[9][1:]
    i_history = top_level[10][1:]

    # calculate % of initial population which is recovered
    rc_history = list ( numpy.array( ra_history ) * (100/M) )
    # same % calculation for immunity history which comes in absolute numbers
    im_history = list ( numpy.array( i_history )  * (100/M) )

    # calculate standard 14 day incidence per 100 000 people
    ic_history = []
    for j in range (0, len(n_history)):
        if ( j < INCIDENCE_PERIOD ):
            my_incidence = 0
        else:
            my_incidence = numpy.array( nc_history[ j - ( INCIDENCE_PERIOD + 1 ) : j-1 ] ).sum() / ( population.value / 0.1 )
        ic_history.append( my_incidence )

    # prevalence = active cases / population * 100
    pr_history = []
    for j in range (0, len(n_history)):
        my_prevalence = ( n_history[j] / ( population.value * 1000000 ) ) * 100
        pr_history.append ( my_prevalence )

    t_transmissions = int(numpy.array(nc_history).sum())
    t_recoveries    = int(numpy.array(r_history).sum())
    t_deaths        = int(numpy.array(d_history).sum())
    p_transmissions = round( (t_transmissions / M )*100 ,2)

    ar_stats = [ t_transmissions, t_recoveries, t_deaths, p_transmissions ]

    # Active, New, Recovered, Dead, Rt, Immunized + accumulated Cases, Recoveries and Deaths + Stats
    return n_history, nc_history, r_history, d_history, rt_history, rc_history, im_history, na_history, ra_history, da_history, ic_history, pr_history, ar_stats

# make an interval :-)
def mki( xa, x, xb, unit='' ):

    return '<b>' + str(x) + unit + '</b>' + ' (' + str(xa) + ' , ' + str(xb) + ')'

# callback function for updating the data
def update_data(attrname, old, new):

    # Generate the new curve with the slider values
    x = np.linspace(0, DAYS, DAYS)
    y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, ar_stats = get_data(x, population.value, iinfections.value, period.value, period_stdev.value, latent.value, DAYS, 0, 0, 0, h1.value*p1.value, 0, 0, DAYS, drate.value, True, im.value )

    # upper and lower values
    p1_ = p1.value*( 1 + p_delta.value/100 )
    _p1 = p1.value*( 1 - p_delta.value/100 )

    y1_, y2_, y3_, y4_, y5_, y6_, y7_, y8_, y9_, y10_, y11_, y12_, ar_stats_ = get_data(x, population.value, iinfections.value, period.value, period_stdev.value, latent.value, DAYS, 0, 0, 0, h1.value*_p1,       0, 0, DAYS, drate.value, True, im.value )
    _y1, _y2, _y3, _y4, _y5, _y6, _y7, _y8, _y9, _y10, _y11, _y12, _ar_stats = get_data(x, population.value, iinfections.value, period.value, period_stdev.value, latent.value, DAYS, 0, 0, 0, h1.value*p1_,       0, 0, DAYS, drate.value, True, im.value )

    # Only the global variable data sources need to be updated
    source_active.data= dict(x=x, y=y1)
    source_new.data    = dict(x=x, y=y2)
    source_rec.data    = dict(x=x, y=y3)
    source_dead.data   = dict(x=x, y=y4)
    source_rt.data     = dict(x=x, y=y5)
    source_rc.data     = dict(x=x, y=y6)
    source_im.data     = dict(x=x, y=y7)
    source_na.data     = dict(x=x, y=y8)
    source_ra.data     = dict(x=x, y=y9)
    source_da.data     = dict(x=x, y=y10)
    source_ic.data     = dict(x=x, y=y11)
    source_pr.data     = dict(x=x, y=y12)

    source_active_.data = data=dict(x=x, y=y1_)
    source_new_.data    = data=dict(x=x, y=y2_)
    source_rec_.data    = data=dict(x=x, y=y3_)
    source_dead_.data   = data=dict(x=x, y=y4_)
    source_rt_.data     = data=dict(x=x, y=y5_)
    source_rc_.data     = data=dict(x=x, y=y6_)
    source_im_.data     = data=dict(x=x, y=y7_)
    source_na_.data     = data=dict(x=x, y=y8_)
    source_ra_.data     = data=dict(x=x, y=y9_)
    source_da_.data     = data=dict(x=x, y=y10_)
    source_ic_.data     = data=dict(x=x, y=y11_)
    source_pr_.data     = data=dict(x=x, y=y12_)

    _source_active.data = data=dict(x=x, y=_y1)
    _source_new.data    = data=dict(x=x, y=_y2)
    _source_rec.data    = data=dict(x=x, y=_y3)
    _source_dead.data   = data=dict(x=x, y=_y4)
    _source_rt.data     = data=dict(x=x, y=_y5)
    _source_rc.data     = data=dict(x=x, y=_y6)
    _source_im.data     = data=dict(x=x, y=_y7)
    _source_na.data     = data=dict(x=x, y=_y8)
    _source_ra.data     = data=dict(x=x, y=_y9)
    _source_da.data     = data=dict(x=x, y=_y10)
    _source_ic.data     = data=dict(x=x, y=_y11)
    _source_pr.data     = data=dict(x=x, y=_y12)

    # Incidence vs Rt
    source_phase_space.data = dict(x=y5, y=y11)

    beta          = round ( h1.value * p1.value / 100 , 4)
    beta_         = round ( h1.value * p1_ / 100 , 4)
    _beta         = round ( h1.value * _p1 / 100 , 4)

    R0            = round ( beta  * period.value , 4)
    R0_           = round ( beta_ * period.value , 4)
    _R0           = round ( _beta * period.value , 4)

    im_threshold  = max (round ( ( 1 - 1/R0  )*100, 2 ),0)
    im_threshold_ = max (round ( ( 1 - 1/R0_ )*100, 2 ),0)
    _im_threshold = max (round ( ( 1 - 1/_R0 )*100, 2 ),0)

    pre_str       = '&beta;: ' + mki(_beta, beta, beta_) + '<br/>R<sub>0</sub>: ' + mki(_R0, R0, R0_)  + '<br/>Immunity threshold: ' + mki(_im_threshold, im_threshold, im_threshold_, '%')
    extra_str     = ''
    stats_str     = pre_str + '<br/>Transmissions: '   + mki(_ar_stats[0], ar_stats[0], ar_stats_[0]) + \
                              '<br/>Transmissions %: ' + mki(_ar_stats[3], ar_stats[3], ar_stats_[3], '%') + \
                              '<br/>Recoveries: '      +  mki(_ar_stats[1], ar_stats[1], ar_stats_[1]) + \
                              '<br/>Deaths: '          +  mki(_ar_stats[2], ar_stats[2], ar_stats_[2]) + extra_str
    stats.text = stats_str

def reset_data():
    population.value   = POP_START
    iinfections.value  = IIF_START
    period.value       = T_START
    period_stdev.value = T_STDEV_START
    latent.value       = L_START
    h1.value           = H1_START
    p1.value           = P1_START
    p_delta.value      = P_DELTA_START
    drate.value        = DRATE_START
    im.value           = IM_START

    # we seem to need to pass something here because the slider callback needs to have a declaration of 3 parameters
    update_data('xxxx',0,0)

def vaccinate_data():
    R0       = h1.value * (p1.value / 100) * period.value
    im.value  = max ((1 - 1 / (R0)) * 100, 0)

    # we seem to need to pass something here because the slider callback needs to have a declaration of 3 parameters
    update_data('xxxx',0,0)

def vaccinate70_data():
    im.value  = 70

    # we seem to need to pass something here because the slider callback needs to have a declaration of 3 parameters
    update_data('xxxx',0,0)

# set properties common to all the plots
def set_plot_details ( aplot, ahover, alabel = PLOT_Y_LABEL ):
    aplot.toolbar.active_drag    = None
    aplot.toolbar.active_scroll  = None
    aplot.toolbar.active_tap     = None

    # add the hover tool
    aplot.add_tools(ahover)
    aplot.toolbar.active_inspect = ahover

    # control placement / visibility of toolbar
    aplot.toolbar_location       = None

    aplot.xaxis.axis_label = PLOT_X_LABEL
    aplot.yaxis.axis_label = alabel

### Main

# Set up widgets
population  = Slider(title=POP_LABEL, value=POP_START, start=POP_MIN, end=POP_MAX, step=POP_STEP)
iinfections = Slider(title=IIF_LABEL, value=IIF_START, start=IIF_MIN, end=IIF_MAX, step=IIF_STEP)

period       = Slider(title=T_LABEL,       value=T_START,       start=T_MIN,       end=T_MAX,       step=1)
period_stdev = Slider(title=T_STDEV_LABEL, value=T_STDEV_START, start=T_STDEV_MIN, end=T_STDEV_MAX, step=1)

latent = Slider(title=L_LABEL, value=L_START, start=L_MIN, end=L_MAX, step=1)

h1 = Slider(title=H1_LABEL, value=H1_START, start=H1_MIN, end=H1_MAX, step=H1_STEP)
p1 = Slider(title=P1_LABEL, value=P1_START, start=P1_MIN, end=P1_MAX, step=P1_STEP)

p_delta = Slider(title=P_DELTA_LABEL, value=P_DELTA_START, start=P_DELTA_MIN, end=P_DELTA_MAX, step=P_DELTA_STEP)

drate = Slider(title=DRATE_LABEL, value=DRATE_START, start=DRATE_MIN, end=DRATE_MAX, step=DRATE_STEP)

im = Slider(title=IM_LABEL, value=IM_START, start=IM_MIN, end=IM_MAX, step=IM_STEP)

button  = Button(label="Reset",                        button_type="default")
button2 = Button(label="Immunize critical proportion", button_type="default")
button3 = Button(label="Immunize 70%",                 button_type="default")

# text widgets
intro   = Div(text='', width=TEXT_WIDTH)
summary = Div(text='', width=TEXT_WIDTH)
stats   = Div(text='', width=TEXT_WIDTH)
notes   = Div(text='', width=TEXT_WIDTH)

# Assign widgets to the call back function
# updates are on value_throtled because this is too slow for realtime updates
for w in [population, iinfections, period, period_stdev, latent, h1, p1, p_delta, drate, im ]:
    w.on_change('value_throttled', update_data)

# reset button call back
button.on_click(reset_data)

# vaccinate the population
button2.on_click(vaccinate_data)
button3.on_click(vaccinate70_data)

# initial plot

p1_ = p1.value*( 1 + p_delta.value/100 )
_p1 = p1.value*( 1 - p_delta.value/100 )

x = np.linspace(1, DAYS, DAYS)
y1,  y2,  y3,  y4,  y5,  y6,  y7,  y8,  y9,  y10,  y11,  y12,  ar_stats = get_data(x, population.value, iinfections.value, period.value, period_stdev.value, latent.value, DAYS, 0, 0, 0,   h1.value*p1.value, 0, 0, DAYS, drate.value, True, im.value )

# upper and lower values
y1_, y2_, y3_, y4_, y5_, y6_, y7_, y8_, y9_, y10_, y11_, y12_, ar_stats_  = get_data(x, population.value, iinfections.value, period.value, period_stdev.value, latent.value, DAYS, 0, 0, 0, h1.value*_p1,       0, 0, DAYS, drate.value, True, im.value )
_y1, _y2, _y3, _y4, _y5, _y6, _y7, _y8, _y9, _y10, _y11, _y12, _ar_stats_ = get_data(x, population.value, iinfections.value, period.value, period_stdev.value, latent.value, DAYS, 0, 0, 0, h1.value*p1_,       0, 0, DAYS, drate.value, True, im.value )

# Active, New, Recovered, Dead, Rt, % Immunine
source_active = ColumnDataSource(data=dict(x=x, y=y1))
source_new    = ColumnDataSource(data=dict(x=x, y=y2))
source_rec    = ColumnDataSource(data=dict(x=x, y=y3))
source_dead   = ColumnDataSource(data=dict(x=x, y=y4))
source_rt     = ColumnDataSource(data=dict(x=x, y=y5))
source_rc     = ColumnDataSource(data=dict(x=x, y=y6))
source_im     = ColumnDataSource(data=dict(x=x, y=y7))
source_na     = ColumnDataSource(data=dict(x=x, y=y8))
source_ra     = ColumnDataSource(data=dict(x=x, y=y9))
source_da     = ColumnDataSource(data=dict(x=x, y=y10))
source_ic     = ColumnDataSource(data=dict(x=x, y=y11))
source_pr     = ColumnDataSource(data=dict(x=x, y=y12))

source_active_ = ColumnDataSource(data=dict(x=x, y=y1_))
source_new_    = ColumnDataSource(data=dict(x=x, y=y2_))
source_rec_    = ColumnDataSource(data=dict(x=x, y=y3_))
source_dead_   = ColumnDataSource(data=dict(x=x, y=y4_))
source_rt_     = ColumnDataSource(data=dict(x=x, y=y5_))
source_rc_     = ColumnDataSource(data=dict(x=x, y=y6_))
source_im_     = ColumnDataSource(data=dict(x=x, y=y7_))
source_na_     = ColumnDataSource(data=dict(x=x, y=y8_))
source_ra_     = ColumnDataSource(data=dict(x=x, y=y9_))
source_da_     = ColumnDataSource(data=dict(x=x, y=y10_))
source_ic_     = ColumnDataSource(data=dict(x=x, y=y11_))
source_pr_     = ColumnDataSource(data=dict(x=x, y=y12_))

_source_active = ColumnDataSource(data=dict(x=x, y=_y1))
_source_new    = ColumnDataSource(data=dict(x=x, y=_y2))
_source_rec    = ColumnDataSource(data=dict(x=x, y=_y3))
_source_dead   = ColumnDataSource(data=dict(x=x, y=_y4))
_source_rt     = ColumnDataSource(data=dict(x=x, y=_y5))
_source_rc     = ColumnDataSource(data=dict(x=x, y=_y6))
_source_im     = ColumnDataSource(data=dict(x=x, y=_y7))
_source_na     = ColumnDataSource(data=dict(x=x, y=_y8))
_source_ra     = ColumnDataSource(data=dict(x=x, y=_y9))
_source_da     = ColumnDataSource(data=dict(x=x, y=_y10))
_source_ic     = ColumnDataSource(data=dict(x=x, y=_y11))
_source_pr     = ColumnDataSource(data=dict(x=x, y=_y12))

# Incidence vs Rt
source_phase_space = ColumnDataSource(data=dict(x=y5, y=y11))

# plot 1

hover = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL, "@y{0}")], mode="mouse" )
hover.point_policy='snap_to_data'
hover.line_policy='nearest'

plot = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )

m = PLOT_LINE_ALPHA_DIFF

plot.line('x', 'y', source=source_active,  line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA,       line_color=PLOT_LINE_ACTIVE_COLOR, legend_label='Active' )
plot.line('x', 'y', source=source_active_, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_ACTIVE_COLOR, legend_label='Active' )
plot.line('x', 'y', source=_source_active, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_ACTIVE_COLOR, legend_label='Active' )

set_plot_details(plot, hover)

# plot 2

# using mode="mouse" because the mouse mode produces overlapping tooltips when multiple lines are used
hover2 = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL, "@y{0}")], mode="mouse" )
hover2.point_policy='snap_to_data'
hover2.line_policy='nearest'

plot2 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT2_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS],)

m = PLOT_LINE_ALPHA_DIFF

plot2.line('x', 'y', source=source_new,  line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA,       line_color=PLOT_LINE_NEW_COLOR,       legend_label='New cases' )
plot2.line('x', 'y', source=source_new_, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_NEW_COLOR,       legend_label='New cases' )
plot2.line('x', 'y', source=_source_new, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_NEW_COLOR,       legend_label='New cases' )

plot2.line('x', 'y', source=source_rec,  line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA,       line_color=PLOT_LINE_RECOVERED_COLOR, legend_label='Recoveries')
plot2.line('x', 'y', source=source_rec_, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_RECOVERED_COLOR, legend_label='Recoveries')
plot2.line('x', 'y', source=_source_rec, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_RECOVERED_COLOR, legend_label='Recoveries')

set_plot_details(plot2, hover2)

# plot 3

# custom precision
hover3 = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL2, "@y{0.00}")], mode="mouse" )
hover3.point_policy='snap_to_data'
hover3.line_policy='nearest'

plot3 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT3_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )

m = PLOT_LINE_ALPHA_DIFF

plot3.line('x', 'y', source=source_rt,  line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA,       line_color=PLOT_LINE_ACTIVE_COLOR, legend_label='Rt' )
plot3.line('x', 'y', source=source_rt_, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_ACTIVE_COLOR, legend_label='Rt' )
plot3.line('x', 'y', source=_source_rt, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_ACTIVE_COLOR, legend_label='Rt' )

set_plot_details(plot3, hover3, PLOT_Y_LABEL2)

# plot 4

# custom precision
hover4 = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL2, "@y{0.00}")], mode="mouse" )
hover4.point_policy='snap_to_data'
hover4.line_policy='nearest'

plot4 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT4_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )

m = PLOT_LINE_ALPHA_DIFF

plot4.line('x', 'y', source=source_rc,  line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA,       line_color=PLOT_LINE_RECOVERED_COLOR, legend_label='% Recovered' )
plot4.line('x', 'y', source=source_rc_, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_RECOVERED_COLOR, legend_label='% Recovered' )
plot4.line('x', 'y', source=_source_rc, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_RECOVERED_COLOR, legend_label='% Recovered' )

plot4.line('x', 'y', source=source_im,  line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA,       line_color=PLOT_LINE_ACTIVE_COLOR,    legend_label='% Immune' )
plot4.line('x', 'y', source=_source_im, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_ACTIVE_COLOR,    legend_label='% Immune' )
plot4.line('x', 'y', source=source_im_, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_ACTIVE_COLOR,    legend_label='% Immune' )

plot4.legend.location = 'bottom_right'

set_plot_details(plot4, hover4, PLOT_Y_LABEL2)

# plot 5

hover5 = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL, "@y{0}")], mode="mouse" )
hover5.point_policy='snap_to_data'
hover5.line_policy='nearest'

plot5 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT5_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )

m = PLOT_LINE_ALPHA_DIFF

plot5.line('x', 'y', source=source_dead,  line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_DEAD_COLOR, legend_label='Deaths' )
plot5.line('x', 'y', source=source_dead_, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_DEAD_COLOR, legend_label='Deaths' )
plot5.line('x', 'y', source=_source_dead, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_DEAD_COLOR, legend_label='Deaths' )

set_plot_details(plot5, hover5)

# plot 6

# using mode="mouse" because the mouse mode produces overlapping tooltips when multiple lines are used
hover6 = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL, "@y{0}")], mode="mouse" )
hover6.point_policy='snap_to_data'
hover6.line_policy='nearest'

plot6 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT6_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS],)

m = PLOT_LINE_ALPHA_DIFF

plot6.line('x', 'y', source=source_na,  line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA,       line_color=PLOT_LINE_NEW_COLOR,       legend_label='Cases' )
plot6.line('x', 'y', source=_source_na, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_NEW_COLOR,       legend_label='Cases' )
plot6.line('x', 'y', source=source_na_, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_NEW_COLOR,       legend_label='Cases' )


plot6.line('x', 'y', source=source_ra,  line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA,       line_color=PLOT_LINE_RECOVERED_COLOR, legend_label='Recoveries')
plot6.line('x', 'y', source=_source_ra, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_RECOVERED_COLOR, legend_label='Recoveries')
plot6.line('x', 'y', source=source_ra_, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_RECOVERED_COLOR, legend_label='Recoveries')
plot6.legend.location = 'bottom_right'

set_plot_details(plot6, hover6)

# plot 7

hover7 = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL, "@y{0}")], mode="mouse" )
hover7.point_policy='snap_to_data'
hover7.line_policy='nearest'

plot7 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT7_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )

m = PLOT_LINE_ALPHA_DIFF

plot7.line('x', 'y', source= source_da, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA,       line_color=PLOT_LINE_DEAD_COLOR, legend_label='Dead' )
plot7.line('x', 'y', source=_source_da, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_DEAD_COLOR, legend_label='Dead' )
plot7.line('x', 'y', source=source_da_, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_DEAD_COLOR, legend_label='Dead' )

plot7.legend.location = 'bottom_right'

set_plot_details(plot7, hover7)

# plot 8

hover8 = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL2, "@y{0.00}")], mode="mouse" )
hover8.point_policy='snap_to_data'
hover8.line_policy='nearest'

plot8 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT8_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )

m = PLOT_LINE_ALPHA_DIFF

plot8.line('x', 'y', source=source_ic,  line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA,       line_color=PLOT_LINE_NEW_COLOR, legend_label='Incidence' )
plot8.line('x', 'y', source=_source_ic, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_NEW_COLOR, legend_label='Incidence' )
plot8.line('x', 'y', source=source_ic_, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_NEW_COLOR, legend_label='Incidence' )

set_plot_details(plot8, hover8, PLOT_Y_LABEL2)

# plot 9

hover9 = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL2, "@y{0.00}")], mode="mouse" )
hover9.point_policy='snap_to_data'
hover9.line_policy='nearest'

plot9 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT9_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )

m = PLOT_LINE_ALPHA_DIFF

plot9.line('x', 'y', source=source_pr,  line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA,       line_color=PLOT_LINE_NEW_COLOR, legend_label='% Prevalence' )
plot9.line('x', 'y', source=_source_pr, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_NEW_COLOR, legend_label='% Prevalence' )
plot9.line('x', 'y', source=source_pr_, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA*(1-m), line_color=PLOT_LINE_NEW_COLOR, legend_label='% Prevalence' )

set_plot_details(plot9, hover9, PLOT_Y_LABEL2)

# misc text
intro.text    = TEXT_INTRO
summary.text  = TEXT_SUMMARY
notes.text    = TEXT_NOTES

summary.style = { 'font-weight' : 'bold' }

beta = round ( h1.value * p1.value / 100 , 4)
R0   = round ( beta * period.value , 4)

# let's force an update so the stats are written, etc
update_data('xxxx',0,0)

# plot 10 - test plot for evolution of risk

PLOT10_TITLE = 'Risk evolution'
RT_LIM = 1
IC_LIM = 120
PLOT10_ALPHA  = 0.2

hover10 = HoverTool(tooltips=[ ('Rt', "@x{0}"), ('Incidence', "@y{0.00}")], mode="mouse" )
hover10.point_policy='snap_to_data'
hover10.line_policy='nearest'

plot10 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT10_TITLE, tools=PLOT_TOOLS, x_range=[0, R0], )

plot10.line('x', 'y', source=source_phase_space, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_DEAD_COLOR, )

set_plot_details(plot10, hover10, PLOT_Y_LABEL2)

plot10.xaxis.axis_label = 'Rt'
plot10.yaxis.axis_label = '14 day Incidence per 100m'

quadrant3 = BoxAnnotation(left=0, right=RT_LIM, bottom=0, top=IC_LIM, fill_alpha=PLOT10_ALPHA, fill_color='green')
quadrant4 = BoxAnnotation(left=1, right=R0    , bottom=0, top=IC_LIM, fill_alpha=PLOT10_ALPHA, fill_color='orange')
quadrant1 = BoxAnnotation(left=1, right=R0    , bottom=IC_LIM,        fill_alpha=PLOT10_ALPHA, fill_color='red')
quadrant2 = BoxAnnotation(left=0, right=1     , bottom=IC_LIM,        fill_alpha=PLOT10_ALPHA, fill_color='orange')

plot10.add_layout(quadrant1)
plot10.add_layout(quadrant2)
plot10.add_layout(quadrant3)
plot10.add_layout(quadrant4)

# Set up layouts and add to document
notespacer = Spacer(width=TEXT_WIDTH, height=10, width_policy='auto', height_policy='fixed')

# simplified set
inputs = column(intro, population, iinfections, period, latent, h1, p1, p_delta, drate, im, button2, button3, button, summary, stats, notespacer, notes)

curdoc().title = PAGE_TITLE

if ENABLE_RISK:
    last_plot = plot10
else:
    last_plot = plot9

# useful for mobile scrolling on the left side
leftmargin = Spacer(width=LMARGIN_WIDTH, height=400, width_policy='fixed', height_policy='auto')
curdoc().add_root( row(leftmargin,inputs, column(plot, plot2, plot5), column(plot4, plot6, plot7), column(plot3, plot8, last_plot)) )
