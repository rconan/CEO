import plotly.plotly as ply
import plotly.graph_objs as plyg

#from IPython.display import display, clear_output, Javascript
#display(Javascript('IPython.notebook.kernel.execute("theNotebook = " + "\'"+IPython.notebook.notebook_name+"\'");'))

def heatmap(z,x=None,y=None,filename="Ipython Notebooks"):

    if x is None:
        x = range(z.shape[1])
    if y is None:
        y = range(z.shape[0])

    r = float(z.shape[1])/float(z.shape[0])

    trace = plyg.Heatmap( z=z, x=x, y=y,
                          colorscale="Portland")
    data = plyg.Data([trace])
    wh = 300
    mg = 200
    height = wh + mg
    width = wh*r + mg
        
    fig = plyg.Figure(data=data,layout=plyg.Layout(autosize=False,height=height,width=width))

    return ply.iplot(fig,filename=filename)
