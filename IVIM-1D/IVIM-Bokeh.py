#Author: Gustavo Solcia
#E-mail: gustavo.solcia@usp.br

"""Simple python code to plot an interactive graph of IVIM multiexponential model using Bokeh. 

"""

import numpy as np
from bokeh.layouts import column, row
from bokeh.models import CustomJS, Slider
from bokeh.plotting import ColumnDataSource, figure, show

PF0 = 1.5e-1
DC0 = 4.0E-2
PDC0 = 4.0E-2

x = np.linspace(0, 900, 500)
y = (1-PF0)*np.exp(-DC0*x) + PF0*np.exp(-PDC0*x)

source = ColumnDataSource(data=dict(x=x, y=y))

plot = figure(x_range=(1e0,1e3),y_range=(1e-6, 1e0),
         y_axis_type='log',
        plot_width=400, plot_height=400)

plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)

PF_slider = Slider(start=0.1, end=1, value=PF0, step=.01, title="Perfusion fraction")
DC_slider = Slider(start=5e-4, end=5e-3, value=DC0, step=1e-5, title="Diffusion coefficient")
PDC_slider = Slider(start=6e-3, end=9e-2, value=PDC0, step=1e-4, title="Pseudodiffusion coefficient")

callback = CustomJS(args=dict(source=source, pf=PF_slider, dc=DC_slider, pdc=PDC_slider),
                    code="""
    const data = source.data;
    const A = pf.value;
    const k = dc.value;
    const phi = pdc.value;
    const x = data['x']
    const y = data['y']
    const y1 = data['y']
    for (var i = 0; i < x.length; i++) {
        y[i] = (1-A)*Math.exp(-k*x[i])+A*Math.exp(-phi*x[i]);
    }
    source.change.emit();
""")

PF_slider.js_on_change('value', callback)
DC_slider.js_on_change('value', callback)
PDC_slider.js_on_change('value', callback)

layout = row(
    plot,
    column(PF_slider, DC_slider, PDC_slider),
)

show(layout)
