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
IIF_MAX   = 50
IIF_START = 4

# Infectious period
T_MIN   = 7
T_MAX   = 50
T_START = 19

# and its Standard Deviation
T_STDEV_MIN   = 0
T_STDEV_MAX   = 3
T_STDEV_START = 2

# Incubation period
I_MIN   = 1
I_MAX   = 20
I_START = 1

# Stage durations 
DUR_MIN  = 1
DUR1_MAX = 30
DUR2_MAX = 60

DUR1_START = 20
DUR2_START = 55

# Transition durations
TRA_MIN  = 0
TRA1_MAX = 45
TRA2_MAX = 60

TRA1_START = 18
TRA2_START = 40

# Simulation time
DAYS = 150 # DUR1_MAX + DUR2_MAX + DUR3 < DAYS 

# Propagation rate parameters
#
# Note:
#   BETA = hp
#   R0   = hpT
#
BETA_MIN  =  0

BETA1_MAX   = 0.6  * 10
BETA1_START = 0.46 * 10
BETA1_STEP  = 0.01

BETA2_MAX   = 0.1   * 10
BETA2_START = 0.042 * 10
BETA2_STEP  = 0.01

BETA3_MAX   = 0.6   * 10
BETA3_START = 0.057 * 10
BETA3_STEP  = 0.01

DRATE_MIN   = 0
DRATE_MAX   = 50
DRATE_START = 0.85
DRATE_STEP  = 0.05

# labels and strings
PAGE_TITLE  ='3 stage epidemic simulator'
PLOT_TITLE  ='Active'
PLOT2_TITLE ='New, Recovered'
PLOT3_TITLE ='Rt estimation'
PLOT4_TITLE ='Immunity'
PLOT5_TITLE ='Dead'
PLOT6_TITLE ='Accumulated cases / recoveries'
PLOT7_TITLE ='Accumulated deaths'

T_LABEL       = 'Infectious Period'
T_STDEV_LABEL = 'Infectious Period Standard Deviation'
I_LABEL       = 'Incubation Period'
POP_LABEL     = 'Population (Millions)'
IIF_LABEL     = 'Initial number of infections'
DUR1_LABEL    = 'Free expansion duration'
DUR2_LABEL    = 'Confinement duration (including transition)'
TRA1_LABEL    = 'Transition to confinement duration'
TRA2_LABEL    = 'Transition to deconfinement duration'
BETA1_LABEL   = 'Beta during free expansion (x10)'
BETA2_LABEL   = 'Beta during confinement (x10)'
BETA3_LABEL   = 'Beta after confinement (x10)'
DRATE_LABEL   = 'Death rate (%)'

TEXT_INTRO    = 'Use the mouse for initial selection and cursors for fine tuning:'
TEXT_SUMMARY  = 'Stats:'
TEXT_NOTES    ='<b>Notes:</b><br/>\
              &bull; Confinement phase is highlighted in red.<br/>\
              &bull; Adjacent highlights represent transition periods.<br/>\
              &bull; Pre/post confinement phases appear in white.<br/>\
              &bull; For each phase &beta; = hp.<br/>\
              &bull; More info at <a href="https://github.com/ghomem/viraly">github.com/ghomem/viraly</a>'
### End of configuration

### Functions

# the function that we are plotting
def get_data(x, pop, n0, period, period_stdev, incubation, d1, d2, tr1, tr2, b1, b2,b3, tmax, dr, prog_change ):

    h  = 1
    p  = float (b1 / 10) # input is multiplied by 10 for precision on the sliders
    h2 = 1
    p2 = float (b2 / 10) # input is multiplied by 10 for precision on the sliders
    h3 = 1
    p3 = float (b3 / 10) # input is multiplied by 10 for precision on the sliders
    T  = period
    I  = incubation
    N0 = n0
    DR = float(dr/100) # input is in percentage
    M  =  pop*1000000   # input is in millions
    L  = period_stdev

    tint  = d1
    tint2 = d1 + d2
    progressive = prog_change # bool
    ttime  = tr1
    ttime2 = tr2

    str_params = '{h},{p},{T},{L},{I},{h2},{p2},{tint},{tmax},{M},{N0},{DR},{progressive},{ttime},{h3},{p3},{tint2},{ttime2}'.format(h=h, p=p, T=T, L=L, I=I, h2=h2,p2=p2, \
                                                                                                                          tint=tint, tmax=tmax, M=M, N0=N0, DR=DR, \
                                                                                                                          progressive=progressive, ttime=ttime, h3=h3, p3=p3, tint2=tint2, ttime2=ttime2)
    print(str_params)

    # decide which model to use based on the value of L
    if L == 0:
        prefer_mod4 = False
    else:
        prefer_mod4 = True

    # this function is included from viraly.py
    top_level = run_simulation ( h, p, T, L, I, h2, p2, tint, tmax, M, N0, DR, progressive, ttime, h3, p3, tint2, ttime2, True, prefer_mod4 )

    # dataset from viraly.py: [ n_history, nc_history, list(r_history), list(d_history), m_history, n_history, ra_history, da_history, rt_history, na_history ]
    n_history  = top_level[0]
    nc_history = top_level[1]
    r_history  = top_level[2]
    d_history  = top_level[3]
    ra_history = top_level[6]
    da_history = top_level[7]
    rt_history = top_level[8]
    na_history = top_level[9]

    # calculate % of initial population which is immunized
    im_history = list ( numpy.array( ra_history ) * (100/M) )

    t_transmissions = int(numpy.array(nc_history).sum())
    t_recoveries    = int(numpy.array(r_history).sum())
    t_deaths        = int(numpy.array(d_history).sum())

    ar_stats = [ t_transmissions, t_recoveries, t_deaths ]

    # Active, New, Recovered, Dead, Rt, Immunized
    return n_history, nc_history, r_history, d_history, rt_history, im_history, na_history, ra_history, da_history, ar_stats

# callback function dor updating the data
def update_data(attrname, old, new):

    # Generate the new curve with the slider values
    x = np.linspace(0, DAYS, DAYS)
    y1, y2, y3, y4, y5, y6, y7, y8, y9, ar_stats = get_data(x, population.value, iinfections.value, period.value, period_stdev.value, incubation.value, duration1.value, duration2.value, transition1.value, transition2.value, beta1.value, beta2.value, beta3.value, DAYS, drate.value, True )

    # Only the global variable data sources need to be updated
    source_active.data = dict(x=x, y=y1)
    source_new.data    = dict(x=x, y=y2)
    source_rec.data    = dict(x=x, y=y3)
    source_dead.data   = dict(x=x, y=y4)
    source_rt.data     = dict(x=x, y=y5)
    source_im.data     = dict(x=x, y=y6)
    source_na.data     = dict(x=x, y=y7)
    source_ra.data     = dict(x=x, y=y8)
    source_da.data     = dict(x=x, y=y9)

    transition1_begin = duration1.value
    transition1_end   = transition1_begin + transition1.value
    confinement_begin = transition1_end
    confinement_end   = confinement_begin + ( duration2.value - transition1.value )
    transition2_begin = confinement_end
    transition2_end   = transition2_begin + transition2.value

    transition1_box.left  = transition1_begin
    transition1_box.right = transition1_end
    confinement_box.left  = confinement_begin
    confinement_box.right = confinement_end
    transition2_box.left  = transition2_begin
    transition2_box.right = transition2_end

    stats_str     = 'Transmissions: ' + str(ar_stats[0]) + '<br/>Recoveries: ' + str(ar_stats[1]) + '<br/>Deaths: ' + str(ar_stats[2])
    stats.text = stats_str

def reset_data():
    population.value   = POP_START
    iinfections.value  = IIF_START
    period.value       = T_START
    period_stdev.value = T_STDEV_START
    incubation.value   = I_START
    duration1.value    = DUR1_START
    duration2.value    = DUR2_START
    transition1.value  = TRA1_START
    transition2.value  = TRA2_START
    beta1.value        = BETA1_START
    beta2.value        = BETA2_START
    beta3.value        = BETA3_START
    drate.value        = DRATE_START

    # we seem to need to pass something here because the slider callback needs to have a declaration of 3 parameters
    update_data('xxxx',0,0)

### Main

# Set up widgets
population  = Slider(title=POP_LABEL, value=POP_START, start=POP_MIN, end=POP_MAX, step=POP_STEP)
iinfections = Slider(title=IIF_LABEL, value=IIF_START, start=IIF_MIN, end=IIF_MAX, step=1)

period       = Slider(title=T_LABEL,       value=T_START,       start=T_MIN,       end=T_MAX,       step=1)
period_stdev = Slider(title=T_STDEV_LABEL, value=T_STDEV_START, start=T_STDEV_MIN, end=T_STDEV_MAX, step=1)

incubation = Slider(title=I_LABEL, value=I_START, start=I_MIN, end=I_MAX, step=1)

duration1 = Slider(title=DUR1_LABEL, value=DUR1_START, start=DUR_MIN, end=DUR1_MAX, step=1)
duration2 = Slider(title=DUR2_LABEL, value=DUR2_START, start=DUR_MIN, end=DUR2_MAX, step=1)

transition1 = Slider(title=TRA1_LABEL, value=TRA1_START, start=TRA_MIN, end=TRA1_MAX, step=1)
transition2 = Slider(title=TRA2_LABEL, value=TRA2_START, start=TRA_MIN, end=TRA2_MAX, step=1)

beta1 = Slider(title=BETA1_LABEL, value=BETA1_START, start=BETA_MIN, end=BETA1_MAX, step=BETA1_STEP)
beta2 = Slider(title=BETA2_LABEL, value=BETA2_START, start=BETA_MIN, end=BETA2_MAX, step=BETA2_STEP)
beta3 = Slider(title=BETA3_LABEL, value=BETA3_START, start=BETA_MIN, end=BETA3_MAX, step=BETA3_STEP)

drate = Slider(title=DRATE_LABEL, value=DRATE_START, start=DRATE_MIN, end=DRATE_MAX, step=DRATE_STEP)

button = Button(label="Reset", button_type="default")

# text widgets
intro   = Div(text='', width=TEXT_WIDTH)
summary = Div(text='', width=TEXT_WIDTH)
stats   = Div(text='', width=TEXT_WIDTH)
notes   = Div(text='', width=TEXT_WIDTH)

# Assign widgets to the call back function
# updates are on value_throtled because this is too slow for realtime updates
for w in [population, iinfections, period, period_stdev, incubation, duration1, duration2, transition1, transition2, beta1, beta2, beta3, drate ]:
    w.on_change('value_throttled', update_data)

# reset button call back
button.on_click(reset_data)

# initial plot
x = np.linspace(0, DAYS, DAYS)
y1, y2, y3, y4, y5, y6, y7, y8, y9, ar_stats = get_data(x, population.value, iinfections.value, period.value, period_stdev.value, incubation.value, duration1.value, duration2.value, transition1.value, transition2.value, beta1.value, beta2.value, beta3.value, DAYS, drate.value, True )

# Active, New, Recovered, Dead, Rt, % Immunine
source_active = ColumnDataSource(data=dict(x=x, y=y1))
source_new    = ColumnDataSource(data=dict(x=x, y=y2))
source_rec    = ColumnDataSource(data=dict(x=x, y=y3))
source_dead   = ColumnDataSource(data=dict(x=x, y=y4))
source_rt     = ColumnDataSource(data=dict(x=x, y=y5))
source_im     = ColumnDataSource(data=dict(x=x, y=y6))
source_na     = ColumnDataSource(data=dict(x=x, y=y7))
source_ra     = ColumnDataSource(data=dict(x=x, y=y8))
source_da     = ColumnDataSource(data=dict(x=x, y=y9))

# plot 1

hover = HoverTool(tooltips=[ (PLOT_X_LABEL, "$index"), (PLOT_Y_LABEL, "@y{0}")], mode="vline" )
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
hover2 = HoverTool(tooltips=[ (PLOT_X_LABEL, "$index"), (PLOT_Y_LABEL, "@y{0}")], mode="mouse" )
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
hover3 = HoverTool(tooltips=[ (PLOT_X_LABEL, "$index"), (PLOT_Y_LABEL2, "@y{0.00}")], mode="vline" )
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
hover4 = HoverTool(tooltips=[ (PLOT_X_LABEL, "$index"), (PLOT_Y_LABEL2, "@y{0.00}")], mode="vline" )
hover4.point_policy='snap_to_data'
hover4.line_policy='nearest'

plot4 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT4_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )
plot4.xaxis.axis_label = PLOT_X_LABEL
plot4.yaxis.axis_label = PLOT_Y_LABEL2
plot4.add_tools(hover4)
plot4.toolbar.active_inspect = None

plot4.line('x', 'y', source=source_im, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_ACTIVE_COLOR, legend_label='% Immune' )
plot4.legend.location = 'bottom_right'

# plot 5

hover5 = HoverTool(tooltips=[ (PLOT_X_LABEL, "$index"), (PLOT_Y_LABEL, "@y{0}")], mode="vline" )
hover5.point_policy='snap_to_data'
hover5.line_policy='nearest'

plot5 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT5_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )
plot5.xaxis.axis_label = PLOT_X_LABEL
plot5.yaxis.axis_label = PLOT_Y_LABEL
plot5.add_tools(hover5)
plot5.toolbar.active_inspect = None

plot5.line('x', 'y', source=source_dead, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_DEAD_COLOR, legend_label='Dead' )

# plot 6

# using mode="mouse" because the vline mode produces overlapping tooltips when multiple lines are used
hover6 = HoverTool(tooltips=[ (PLOT_X_LABEL, "$index"), (PLOT_Y_LABEL, "@y{0}")], mode="mouse" )
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

hover7 = HoverTool(tooltips=[ (PLOT_X_LABEL, "$index"), (PLOT_Y_LABEL, "@y{0}")], mode="vline" )
hover7.point_policy='snap_to_data'
hover7.line_policy='nearest'

plot7 = figure(plot_height=PLOT_HEIGHT, plot_width=PLOT_WIDTH, title=PLOT7_TITLE, tools=PLOT_TOOLS, x_range=[0, DAYS], )
plot7.xaxis.axis_label = PLOT_X_LABEL
plot7.yaxis.axis_label = PLOT_Y_LABEL
plot7.add_tools(hover7)
plot7.toolbar.active_inspect = None

plot7.line('x', 'y', source=source_da, line_width=PLOT_LINE_WIDTH, line_alpha=PLOT_LINE_ALPHA, line_color=PLOT_LINE_DEAD_COLOR, legend_label='Deaths' )
plot7.legend.location = 'bottom_right'

# highlight phases with boxes
transition1_begin = duration1.value
transition1_end   = transition1_begin + transition1.value
confinement_begin = transition1_end
confinement_end   = confinement_begin + ( duration2.value - transition1.value)
transition2_begin = confinement_end
transition2_end   = transition2_begin + transition2.value

transition1_box = BoxAnnotation(left=transition1_begin, right=transition1_end, fill_alpha=0.1, fill_color='red')
confinement_box = BoxAnnotation(left=confinement_begin, right=confinement_end, fill_alpha=0.2, fill_color='red')
transition2_box = BoxAnnotation(left=transition2_begin, right=transition2_end, fill_alpha=0.1, fill_color='yellow')

plot.add_layout(transition1_box)
plot.add_layout(confinement_box)
plot.add_layout(transition2_box)

plot2.add_layout(transition1_box)
plot2.add_layout(confinement_box)
plot2.add_layout(transition2_box)

plot3.add_layout(transition1_box)
plot3.add_layout(confinement_box)
plot3.add_layout(transition2_box)

plot4.add_layout(transition1_box)
plot4.add_layout(confinement_box)
plot4.add_layout(transition2_box)

plot5.add_layout(transition1_box)
plot5.add_layout(confinement_box)
plot5.add_layout(transition2_box)

plot6.add_layout(transition1_box)
plot6.add_layout(confinement_box)
plot6.add_layout(transition2_box)

plot7.add_layout(transition1_box)
plot7.add_layout(confinement_box)
plot7.add_layout(transition2_box)

# misc text
intro.text    = TEXT_INTRO
summary.text  = TEXT_SUMMARY
summary.style = { 'font-weight' : 'bold' }
stats_str     = 'Transmissions: ' + str(ar_stats[0]) + '<br/>Recoveries: ' + str(ar_stats[1]) + '<br/>Deaths: ' + str(ar_stats[2])
stats.text    = stats_str
notes.text    = TEXT_NOTES

# Set up layouts and add to document
notespacer = Spacer(width=TEXT_WIDTH, height=10, width_policy='auto', height_policy='fixed')
inputs = column(intro, population, iinfections, period, period_stdev, incubation, duration1, transition1, duration2, transition2, beta1, beta2, beta3, drate, button, summary, stats, notespacer, notes)

curdoc().title = PAGE_TITLE

# useful for mobile scrolling on the left side
leftmargin = Spacer(width=LMARGIN_WIDTH, height=400, width_policy='fixed', height_policy='auto')
curdoc().add_root( row(leftmargin,inputs, column(plot, plot2, plot5), column(plot4, plot6, plot7), column(plot3)) )
