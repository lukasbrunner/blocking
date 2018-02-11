import numpy as np
import matplotlib as mpl
import cartopy.crs as ccrs

pcolormesh = dict(
    cbar_kwargs=dict(
        label= '',
        fraction=0.024,
        pad=0.02,
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
    left=.1,
    right=.92,
    bottom=0.1,
    top=.95

)

# polygons = [
#     mpl.patches.Rectangle(
#         xy=[-30, 40],
#         width=90,
#         height=35,
#         facecolor='none',
#         edgecolor='black',
#         lw=1,
#         transform=ccrs.PlateCarree()
#     )
# ]
