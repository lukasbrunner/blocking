import numpy as np
import matplotlib as mpl
import cartopy.crs as ccrs

boundaries = [0, .5, 1.1]

pcolormesh = dict(
    colors = ['white', 'red'],
    levels = boundaries,
    cbar_kwargs=dict(
        label= '',
        # extend = 'max',
        ticks = [.25, .75]
    )
)

subplots = dict(
    figsize = (8, 4),
    subplot_kw = dict(
        projection =
        # ccrs.Robinson(central_longitude=-90)
        # ccrs.Mollweide()
        ccrs.PlateCarree(central_longitude=0)
    )
)

subplots_adjust = dict(
    hspace=0,
    wspace=0,
    # left=0.1,
    # right=.88,
    # bottom=0.12,
    # top=.92
    left=.1,
    right=.85,
    bottom=0.1,
    top=.95

)

polygons = [
    mpl.patches.Rectangle(
        xy=[-180, -90],
        width=360,
        height=15,
        facecolor='black',
        edgecolor='none',
        alpha=.5,
        transform=ccrs.PlateCarree()
    ),
    mpl.patches.Rectangle(
        xy=[-180, 75],
        width=360,
        height=15,
        facecolor='black',
        edgecolor='none',
        alpha=.5,
        transform=ccrs.PlateCarree()
    ),
        mpl.patches.Rectangle(
        xy=[-180, -15],
        width=360,
        height=30,
        facecolor='black',
        edgecolor='none',
        alpha=.5,
        transform=ccrs.PlateCarree()
    )
]
