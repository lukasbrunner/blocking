import numpy as np
import matplotlib as mpl

levels = np.arange(0., .031, .001)
cmap = mpl.pyplot.cm.get_cmap('YlOrRd')
cmap.set_under('white')

# TODO: set this to the last color of cmap
# NOTE: without this the color of extend will be the same as the one before
cmap.set_over('maroon')

pcolormesh = dict(
    cmap = cmap,
    vmin = levels[1],
    vmax = levels[-1],
    norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=False),
    extend='max',
    cbar_kwargs=dict(
        label= '',
        # extend = 'max',
        ticks = levels[::5]
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
    right=.87,
    bottom=0.1,
    top=.95

)


def kwargs_blocking(args):
    if args.range is None:
        if args.varn == 'Blocking':
            levels = np.arange(0., .03+.0001, .001)
            ticks = np.arange(0, .03+.001, .005)
        elif args.varn == 'IB':
            levels = np.arange(0., .08+.001, .002)
            ticks = np.arange(0, .08+.001, .01)
    else:
        levels = np.arange(*args.range)
        if args.ticks is None:
            ticks = None
        else:
            ticks = np.arange(*args.ticks)

    print levels, ticks

    kwargs = {}
    #kwargs['levels'] = levels
    cmap = plt.cm.get_cmap('YlOrRd')
    #colors = cmap(np.linspace(0, len(levels)/(len(levels)+1.), cmap.N))
    # https://stackoverflow.com/questions/40982050/matplotlib-how-to-cut-the-unwanted-part-of-a-colorbar#40983666
    #cmap_new = mpl.colors.LinearSegmentedColormap.from_list('temp', colors)
    # cmap.set_under('w')
    # cmap_new.set_over(cmap(np.linspace(0, 1, len(levels)+10)[-1]))
    # kwargs['extend'] = 'max'
    kwargs['cmap'] = cmap
    # kwargs['vmin'] = levels[1]
    import ipdb; ipdb.set_trace()
    levels = MaxNLocator(nbins=15).tick_values(levels[0], levels[-1])


    # http://matplotlib.org/examples/images_contours_and_fields/pcolormesh_levels.html
    kwargs['norm'] = mpl.colors.BoundaryNorm(
         levels, ncolors=cmap.N, clip=False)

    return kwargs
