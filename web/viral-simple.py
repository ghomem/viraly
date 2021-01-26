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

# Geometry and visuals
PLOT_TOOLS    ='save,reset,pan,wheel_zoom,box_zoom'
PLOT_HEIGHT   = 300
PLOT_WIDTH    = 500
TEXT_WIDTH    = 300
LMARGIN_WIDTH = 20

PLOT_LINE_WIDTH = 3
PLOT_LINE_ALPHA = 0.6

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
IIF_MAX   = 1000
IIF_START = 100

# Infectious period
T_MIN   = 1
T_MAX   = 15
T_START = 5

# and its Standard Deviation
T_STDEV_MIN   = 0
T_STDEV_MAX   = 3
T_STDEV_START = 0

# Latent period
# Not the same as incubation period, and can be shorter:
# https://en.wikipedia.org/wiki/Latent_period_(epidemiology)
L_MIN   = 1
L_MAX   = 20
L_START = 1

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
H1_START = 17
H1_STEP  = 0.1

P1_MIN   = 0
P1_MAX   = 100
P1_START = 0.0165 * 100
P1_STEP  = 0.01

DRATE_MIN   = 0.05
DRATE_MAX   = 10
DRATE_START = 0.25
DRATE_STEP  = 0.05

IM_MIN   = 0
IM_MAX   = 100
IM_START = 0
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
DRATE_LABEL   = 'Death rate (%)'
IM_LABEL      = 'Pre immunized (%)'

TEXT_INTRO    = 'Use the mouse for initial selection and cursors for fine tuning:'
TEXT_SUMMARY  = 'Stats:'
TEXT_NOTES    ='<b>Notes:</b><br/>\
              &bull; &beta; = hp<br/>\
              &bull; R0 = hpT<br/>\
              &bull; Technical info at <a href="https://github.com/ghomem/viraly">github.com/ghomem/viraly</a>'
### End of configuration

### Functions

# the function that we are plotting
def get_data(x, pop, n0, period, period_stdev, latent, d1, d2, tr1, tr2, b1, b2,b3, tmax, dr, prog_change, I0 = 0 ):

    h  = 1
    p  = float (b1 / 100) # input is multiplied by 100 for precision on the sliders
    h2 = 1
    p2 = float (b2 / 100) # input is multiplied by 100 for precision on the sliders
    h3 = 1
    p3 = float (b3 / 100) # input is multiplied by 100 for precision on the sliders
    T  = period
    I  = latent
    N0 = n0
    DR = float(dr/100) # input is in percentage
    M  =  pop*1000000   # input is in millions
    L  = period_stdev

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
    str_params = '{h},{p},{T},{L},{I},{h2},{p2},{tint},{tmax},{M},{N0},{DR},{progressive},{ttime},{h3},{p3},{tint2},{ttime2},{prefer_mod4}, {I0}'.format(h=h, p=p, T=T, L=L, I=I, h2=h2,p2=p2,       \
                                                                                                                                                  tint=tint, tmax=tmax, M=M, N0=N0, DR=DR,           \
                                                                                                                                                  progressive=progressive, ttime=ttime, h3=h3, p3=p3,\
                                                                                                                                                  tint2=tint2, ttime2=ttime2, prefer_mod4=prefer_mod4, I0 = I0)
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

    ar_stats = [ t_transmissions, t_recoveries, t_deaths ]

    # Active, New, Recovered, Dead, Rt, Immunized + accumulated Cases, Recoveries and Deaths + Stats
    return n_history, nc_history, r_history, d_history, rt_history, rc_history, im_history, na_history, ra_history, da_history, ic_history, pr_history, ar_stats

# callback function dor updating the data
def update_data(attrname, old, new):

    # Generate the new curve with the slider values
    x = np.linspace(0, DAYS, DAYS)
    y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, ar_stats = get_data(x, population.value, iinfections.value, period.value, period_stdev.value, latent.value, DAYS, 0, 0, 0, h1.value*p1.value, 0, 0, DAYS, drate.value, True, im.value )

    # Only the global variable data sources need to be updated
    source_active.data = dict(x=x, y=y1)
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

    beta          = round ( h1.value * p1.value / 100 , 4)
    R0            = round ( beta * period.value , 4)
    pre_str       = 'Beta: ' + str(beta) + '<br/>R0: ' + str(R0) 
    extra_str     = ''
    stats_str     = pre_str + '<br/>Transmissions: ' + str(ar_stats[0]) + '<br/>Recoveries: ' + str(ar_stats[1]) + '<br/>Deaths: ' + str(ar_stats[2]) + extra_str
    stats.text = stats_str

def reset_data():
    population.value   = POP_START
    iinfections.value  = IIF_START
    period.value       = T_START
    period_stdev.value = T_STDEV_START
    latent.value       = L_START
    h1.value           = H1_START
    p1.value           = P1_START
    drate.value        = DRATE_START
    im.value           = IM_START

    # we seem to need to pass something here because the slider callback needs to have a declaration of 3 parameters
    update_data('xxxx',0,0)

### Main

# Set up widgets
population  = Slider(title=POP_LABEL, value=POP_START, start=POP_MIN, end=POP_MAX, step=POP_STEP)
iinfections = Slider(title=IIF_LABEL, value=IIF_START, start=IIF_MIN, end=IIF_MAX, step=1)

period       = Slider(title=T_LABEL,       value=T_START,       start=T_MIN,       end=T_MAX,       step=1)
period_stdev = Slider(title=T_STDEV_LABEL, value=T_STDEV_START, start=T_STDEV_MIN, end=T_STDEV_MAX, step=1)

latent = Slider(title=L_LABEL, value=L_START, start=L_MIN, end=L_MAX, step=1)

h1 = Slider(title=H1_LABEL, value=H1_START, start=H1_MIN, end=H1_MAX, step=H1_STEP)
p1 = Slider(title=P1_LABEL, value=P1_START, start=P1_MIN, end=P1_MAX, step=P1_STEP)

drate = Slider(title=DRATE_LABEL, value=DRATE_START, start=DRATE_MIN, end=DRATE_MAX, step=DRATE_STEP)

im = Slider(title=IM_LABEL, value=IM_START, start=IM_MIN, end=IM_MAX, step=IM_STEP)

button = Button(label="Reset", button_type="default")

# text widgets
intro   = Div(text='', width=TEXT_WIDTH)
summary = Div(text='', width=TEXT_WIDTH)
stats   = Div(text='', width=TEXT_WIDTH)
notes   = Div(text='', width=TEXT_WIDTH)

# Assign widgets to the call back function
# updates are on value_throtled because this is too slow for realtime updates
for w in [population, iinfections, period, period_stdev, latent, h1, p1, drate, ]:
    w.on_change('value_throttled', update_data)

# reset button call back
button.on_click(reset_data)

# initial plot
x = np.linspace(1, DAYS, DAYS)
y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, ar_stats = get_data(x, population.value, iinfections.value, period.value, period_stdev.value, latent.value, DAYS, 0, 0, 0, h1.value*p1.value, 0, 0, DAYS, drate.value, True, im.value )

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

# plot 1

hover = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL, "@y{0}")], mode="vline" )
hover.point_policy='snap_to_data'
hover.line_policy='nearest'

plot = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )
plot.xaxis.axis_label = PLOT_X_LABEL
plot.yaxis.axis_label = PLOT_Y_LABEL
plot.add_tools(hover)
plot.toolbar.active_inspect = None

plot.line('x', 'y', source=source_active, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_ACTIVE_COLOR, legend_label='Active' )

# plot 2

# using mode="mouse" because the vline mode produces overlapping tooltips when multiple lines are used
hover2 = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL, "@y{0}")], mode="mouse" )
hover2.point_policy='snap_to_data'
hover2.line_policy='nearest'

plot2 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT2_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS],)
plot2.xaxis.axis_label = PLOT_X_LABEL
plot2.yaxis.axis_label = PLOT_Y_LABEL

plot2.line('x', 'y', source=source_new, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_NEW_COLOR,       legend_label='New cases' )
plot2.line('x', 'y', source=source_rec, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_RECOVERED_COLOR, legend_label='Recoveries')
plot2.add_tools(hover2)
plot.toolbar.active_inspect = None

# plot 3

# custom precision
hover3 = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL2, "@y{0.00}")], mode="vline" )
hover3.point_policy='snap_to_data'
hover3.line_policy='nearest'

plot3 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT3_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )
plot3.xaxis.axis_label = PLOT_X_LABEL
plot3.yaxis.axis_label = PLOT_Y_LABEL2
plot3.add_tools(hover3)
plot3.toolbar.active_inspect = None

plot3.line('x', 'y', source=source_rt, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_ACTIVE_COLOR, legend_label='Rt' )

# plot 4

# custom precision
hover4 = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL2, "@y{0.00}")], mode="vline" )
hover4.point_policy='snap_to_data'
hover4.line_policy='nearest'

plot4 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT4_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )
plot4.xaxis.axis_label = PLOT_X_LABEL
plot4.yaxis.axis_label = PLOT_Y_LABEL2
plot4.add_tools(hover4)
plot4.toolbar.active_inspect = None

plot4.line('x', 'y', source=source_rc, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_ACTIVE_COLOR, legend_label='% Recovered' )
plot4.line('x', 'y', source=source_im, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_ACTIVE_COLOR, legend_label='% Immune' )
plot4.legend.location = 'bottom_right'

# plot 5

hover5 = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL, "@y{0}")], mode="vline" )
hover5.point_policy='snap_to_data'
hover5.line_policy='nearest'

plot5 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT5_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )
plot5.xaxis.axis_label = PLOT_X_LABEL
plot5.yaxis.axis_label = PLOT_Y_LABEL
plot5.add_tools(hover5)
plot5.toolbar.active_inspect = None

plot5.line('x', 'y', source=source_dead, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_DEAD_COLOR, legend_label='Deaths' )

# plot 6

# using mode="mouse" because the vline mode produces overlapping tooltips when multiple lines are used
hover6 = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL, "@y{0}")], mode="mouse" )
hover6.point_policy='snap_to_data'
hover6.line_policy='nearest'

plot6 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT6_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS],)
plot6.xaxis.axis_label = PLOT_X_LABEL
plot6.yaxis.axis_label = PLOT_Y_LABEL

plot6.line('x', 'y', source=source_na, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_NEW_COLOR,       legend_label='Cases' )
plot6.line('x', 'y', source=source_ra, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_RECOVERED_COLOR, legend_label='Recoveries')
plot6.legend.location = 'bottom_right'
plot6.add_tools(hover6)
plot.toolbar.active_inspect = None

hover7 = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL, "@y{0}")], mode="vline" )
hover7.point_policy='snap_to_data'
hover7.line_policy='nearest'

plot7 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT7_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )
plot7.xaxis.axis_label = PLOT_X_LABEL
plot7.yaxis.axis_label = PLOT_Y_LABEL
plot7.add_tools(hover7)
plot7.toolbar.active_inspect = None

plot7.line('x', 'y', source=source_da, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_DEAD_COLOR, legend_label='Dead' )
plot7.legend.location = 'bottom_right'

hover8 = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL2, "@y{0.00}")], mode="vline" )
hover8.point_policy='snap_to_data'
hover8.line_policy='nearest'

plot8 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT8_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )
plot8.xaxis.axis_label = PLOT_X_LABEL
plot8.yaxis.axis_label = PLOT_Y_LABEL2
plot8.add_tools(hover8)
plot8.toolbar.active_inspect = None

plot8.line('x', 'y', source=source_ic, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_NEW_COLOR, legend_label='Incidence' )

hover9 = HoverTool(tooltips=[ (PLOT_X_LABEL, "@x{0}"), (PLOT_Y_LABEL2, "@y{0.00}")], mode="vline" )
hover9.point_policy='snap_to_data'
hover9.line_policy='nearest'

plot9 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT9_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )
plot9.xaxis.axis_label = PLOT_X_LABEL
plot9.yaxis.axis_label = PLOT_Y_LABEL2
plot9.add_tools(hover8)
plot9.toolbar.active_inspect = None

plot9.line('x', 'y', source=source_pr, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_NEW_COLOR, legend_label='% Prevalence' )

# misc text
intro.text    = TEXT_INTRO
summary.text  = TEXT_SUMMARY
summary.style = { 'font-weight' : 'bold' }

beta          = round ( h1.value * p1.value / 100 , 4)
R0            = round ( beta * period.value , 4)
pre_str       = 'Beta: ' + str(beta) + '<br/>R0: ' + str(R0)
extra_str     = ''
stats_str     = pre_str + '<br/>Transmissions: ' + str(ar_stats[0]) + '<br/>Recoveries: ' + str(ar_stats[1]) + '<br/>Deaths: ' + str(ar_stats[2]) + extra_str
stats.text = stats_str
notes.text    = TEXT_NOTES

# Set up layouts and add to document
notespacer = Spacer(width=TEXT_WIDTH, height=10, width_policy='auto', height_policy='fixed')

# simplified set for the marketing simulation
inputs = column(intro, population, iinfections, period, h1, p1, drate, im, button, summary, stats, notespacer, notes)

curdoc().title = PAGE_TITLE

# useful for mobile scrolling on the left side
leftmargin = Spacer(width=LMARGIN_WIDTH, height=400, width_policy='fixed', height_policy='auto')
curdoc().add_root( row(leftmargin,inputs, column(plot, plot2, plot5), column(plot4, plot6, plot7), column(plot3, plot8, plot9)) )
