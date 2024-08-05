# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2022 Apostolos Chalkis, Vissarion Fisikopoulos, Elias Tsigaridas

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px
from dingo.utils import compute_copula

def plot_copula(data_flux1, data_flux2, n = 5, width = 900 , height = 600, export_format = "svg"):
    """A Python function to plot the copula between two fluxes

    Keyword arguments:
    data_flux1: A list that contains: (i) the vector of the measurements of the first reaction,
                                      (ii) the name of the first reaction
    data_flux2: A list that contains: (i) the vector of the measurements of the second reaction,
                                      (ii) the name of the second reaction
    n: The number of cells
    """

    flux1 = data_flux1[0]
    flux2 = data_flux2[0]
    copula = compute_copula(flux1, flux2, n)

    fig = go.Figure(
            data   = [go.Surface(z=copula)],
            layout = go.Layout(
                height = height, 
                width  = width,
            )
        )


    fig.update_layout(
            title = 'Copula between '+ data_flux1[1] + ' and ' + data_flux2[1],
            scene = dict(
                    xaxis_title= data_flux1[1],
                    yaxis_title= data_flux2[1],
                    zaxis_title="prob, mass"
                ),
            margin=dict(r=30, b=30, l=30, t=50))

    fig.layout.template = None
    
    fig.show()
    fig_name = data_flux1[1] + "_" + data_flux2[1] + "_copula." + export_format

    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=1.25, y=1.25, z=1.25)
    )

    fig.update_layout(scene_camera=camera)
    fig.to_image(format = export_format, engine="kaleido")
    pio.write_image(fig, fig_name, scale=2)


def plot_histogram(reaction_fluxes, reaction, n_bins=40):
    """A Python function to plot the histogram of a certain reaction flux.

    Keyword arguments:
    reaction_fluxes -- a vector that contains sampled fluxes of a reaction
    reaction -- a string with the name of the reacion
    n_bins -- the number of bins for the histogram
    """

    plt.figure(figsize=(7, 7))

    n, bins, patches = plt.hist(
        reaction_fluxes, bins=n_bins, density=False, facecolor="red", ec="black"
    )

    plt.xlabel("Flux (mmol/gDW/h)", fontsize=16)
    plt.ylabel("Frequency (#samples: " + str(reaction_fluxes.size) + ")", fontsize=14)
    plt.grid(True)
    plt.title("Reaction: " + reaction, fontweight="bold", fontsize=18)
    plt.axis([np.amin(reaction_fluxes), np.amax(reaction_fluxes), 0, np.amax(n) * 1.2])

    plt.show()



def plot_corr_matrix(corr_matrix, reactions, removed_reactions=[], format="svg"):
    """A Python function to plot the heatmap of a model's pearson correlation matrix.

    Keyword arguments:
    corr_matrix -- A matrix produced from the "correlated_reactions" function
    reactions -- A list with the model's reactions
    removed_reactions -- A list with the removed reactions in case of a preprocess.
                         If provided removed reactions are not plotted.
    """
    
    sns_colormap = [[0.0, '#3f7f93'],
                    [0.1, '#6397a7'],
                    [0.2, '#88b1bd'],
                    [0.3, '#acc9d2'],
                    [0.4, '#d1e2e7'],
                    [0.5, '#f2f2f2'],
                    [0.6, '#f6cdd0'],
                    [0.7, '#efa8ad'],
                    [0.8, '#e8848b'],
                    [0.9, '#e15e68'],
                    [1.0, '#da3b46']]
    
    if removed_reactions != 0:
        for reaction in reactions:
            index = reactions.index(reaction)
            if reaction in removed_reactions:
               reactions[index] = None
 
    fig = px.imshow(corr_matrix, 
                    color_continuous_scale = sns_colormap,
                    x = reactions, y = reactions, origin="upper")
    
    fig.update_layout(
    xaxis=dict(tickfont=dict(size=5)),
    yaxis=dict(tickfont=dict(size=5)),
    width=900, height=900, plot_bgcolor="rgba(0,0,0,0)")
    
    fig.update_traces(xgap=1, ygap=1,   hoverongaps=False)
    
    fig.show()
    
    fig_name = "CorrelationMatrix." + format
    pio.write_image(fig, fig_name, scale=2)
