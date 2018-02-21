from __future__ import division
from __future__ import print_function
from plotly.offline import plot
import plotly.graph_objs as go
from cea.plots.variable_naming import LOGO
from cea.plots.color_code import ColorCodeCEA
COLOR = ColorCodeCEA()


def primary_energy(data_frame, analysis_fields, title, output_path):

    #CALCULATE GRAPH
    traces_graph = calc_graph(analysis_fields, data_frame)

    #CALCULATE TABLE
    traces_table = calc_table(analysis_fields, data_frame)

    #PLOT GRAPH
    traces_graph.append(traces_table)
    layout = go.Layout(images=LOGO,title=title, barmode='stack',
                       yaxis=dict(title='Consumption of Fossil Fuels [GJ Oil-eq/yr]', domain=[0.35, 1]),
                       xaxis=dict(title='Scenario Name'))
    fig = go.Figure(data=traces_graph, layout=layout)
    plot(fig, auto_open=False, filename=output_path)


def calc_graph(analysis_fields, data_frame):
    # calculate graph
    graph = []
    data_frame['total'] = total = data_frame[analysis_fields].sum(axis=1)
    data_frame = data_frame.sort_values(by='total', ascending=False) # this will get the maximum value to the left
    for field in analysis_fields:
        y = data_frame[field]
        total_perc = (y/total*100).round(2).values
        total_perc_txt = ["("+str(x)+" %)" for x in total_perc]
        trace = go.Bar(x=data_frame.index, y=y, name=field, text = total_perc_txt,
                       marker=dict(color=COLOR.get_color_rgb(field.split('_GJ', 1)[0])))
        graph.append(trace)

    return graph

def calc_table(analysis_fields, data_frame):
    median = data_frame[analysis_fields].median().round(2).tolist()
    median = [str('{:20,.2f}'.format(x)) for x in median]
    # calculate graph
    anchors = []
    for field in analysis_fields+["total"]:
        anchors.append(calc_top_three_anchor_loads(data_frame, field))

    table = go.Table(domain=dict(x=[0, 1], y=[0.0, 0.2]),
                            header=dict(values=['Category', 'Median [GJ Oil-eq/yr]', 'Top 3 most consuming scenarios']),
                            cells=dict(values=[analysis_fields, median, anchors]))

    return table

def calc_top_three_anchor_loads(data_frame, field):
    data_frame = data_frame.sort_values(by=field, ascending=False)
    anchor_list = data_frame[:3].index.values
    return anchor_list