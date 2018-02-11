import numpy as np
import matplotlib as mpl
import cartopy.crs as ccrs

boundaries = np.arange(0., .031, .001)
cmap = mpl.pyplot.cm.get_cmap('YlOrRd', len(boundaries))
colors = list(cmap(np.arange(len(boundaries))))
colors[0] = 'white'
cmap.set_under(colors[-1])
cmap = mpl.colors.ListedColormap(colors[:-1])

pcolormesh = dict(
    cmap = cmap,
    norm = mpl.colors.BoundaryNorm(boundaries, ncolors=cmap.N, clip=False),
    extend='max',
    cbar_kwargs=dict(
        label= '',
        fraction=0.0228,
        pad=0.02,
        # extend = 'max',
        ticks = boundaries[::5]
    )
)

subplots = dict(
    figsize = (8, 4),
    subplot_kw = dict(
        projection =
        # ccrs.Robinson(central_longitude=-90)
        # ccrs.Mollweide()
        ccrs.PlateCarree(central_longitude=-90)
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
    right=.92,
    bottom=0.1,
    top=.95

)

polygons = [
    mpl.patches.Rectangle(
        xy=[-30, 40],
        width=90,
        height=35,
        facecolor='none',
        edgecolor='black',
        lw=1,
        transform=ccrs.PlateCarree()
    )
]
