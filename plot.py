#plot.py

import plotly.graph_objects as go

from plotly.subplots import make_subplots
from heat_transfer import *

rs = [ rc + min(i,Na)*dra + max(0,i-Na)*drp for i in range(N) ]

fig = go.Figure()

fig.add_trace(go.Scatter(
  x=rs,y=T,
  mode="lines",
  showlegend = False,
))

fig.update_layout(
  title="Distribuição de temperatura no forninho",
  xaxis = dict(
    title= "Posição radial [m]",
    # gridcolor = "lightgrey",
  ),
  yaxis = dict(
    title="Temperatura [°C]",
    # gridcolor = "lightgrey"
  ),
  plot_bgcolor="White",
  paper_bgcolor="White",
  template = "plotly_white",
  width=1618,height=1000,
)

fig.show()