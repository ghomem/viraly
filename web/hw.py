''' 
Hello world

'''
import subprocess
import json
import numpy as np

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider, TextInput, BoxAnnotation, HoverTool, Button, Spacer, Div
from bokeh.plotting import figure

TEXT_WIDTH     = 300
TEXT_STR       = 'Under construction L2'
PAGE_TITLE     = 'L2 Website HW'

hw             = Div(text=TEXT_STR, width=TEXT_WIDTH)
curdoc().title = PAGE_TITLE

curdoc().add_root( row( hw )) 
