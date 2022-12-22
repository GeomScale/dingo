# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2022 Apostolos Chalkis, Vissarion Fisikopoulos, Elias Tsigaridas

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from dingo.utils import compute_copula

def plot_copula(data_flux1, data_flux2, n = 5):

    flux1 = data_flux1[0]
    flux2 = data_flux2[0]
    copula = compute_copula(flux1, flux2, n)

    fig = go.Figure(data=[go.Surface(z=copula)])

    fig.update_layout(title='Copula between '+data_flux1[1]+' and '+data_flux2[1], autosize=False,
                      width=500, height=500, scene = dict(xaxis_title=data_flux1[1],
                      yaxis_title=data_flux2[1], zaxis_title='prob. mass'),
                      margin=dict(l=65, r=50, b=65, t=90))
    
    fig.update_layout()
    fig.show()


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
