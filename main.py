import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Rectangle
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
# from scipy.ndimage import interpolation as ipt
from scipy.ndimage import map_coordinates
import skimage.transform as transform
from shapely.geometry import box
from shapely.affinity import rotate as shapely_rotate
from descartes import PolygonPatch

df = pd.DataFrame.from_dict({
    'file name': [
        'agc4109natural_kmh20kms.edit4_mom0.fits',
        '9244_CLN_V2_0.edit3_mom0.fits',
        'agc6871all_natural_kmh20kms.edit3smo_mom0.fits',
        '9362_CLN_V2_4.edit1_mom0.fits',
        '8408_CLN_V2_0.edit1_mom0.fits',
        'B10_ALL_NAT.APCLN.edit6_mom0.fits',
    ],
    'rotation': [82.83, 39.01, 33.74, 64.5, -37.45, 40.42],
    'pixel size': [1.7, 2, 1.4, 2, 2, 1.7],
    'bar length': [15.0, 15.8, 13.4, 18.2, 10.4, 9.6],
    'beam major': [8.93745, 7.7544, 5.89444, 7.28172, 9.61848, 10.8678],
    'beam minor': [5.09657, 5.27256, 4.7826, 4.8024, 4.7592, 5.08985],
    'beam pa': [-56.1455, 69.09, -32.4369, 62.98, 80.73, -56.165],
    'ra': [119.0693, 216.53495, 178.44321, 218.32296, 200.7516, 185.00551],
    'dec': [11.66217, 5.23793, 10.40315, 3.9034, 13.95065, 8.60778],
})

matplotlib.rcParams.update({'font.size': 12})


def convert_to_solar_mass_per_pc(v):
    # converts from JyHz/beam to solar mass/pc^2 (for the colour bar
    return (v * 1.6737236 * 10**(-24) * 9.5214087 * 10**(36)) / (1.989 * 10**(33))
    return (v * 1.6737236 * 1E-24     * 9.5214087 * 1E36    ) / (1.989 * 1E33    )


def get_galaxy(i):
    # Create a WCS object
    w = WCS(df['file name'][i], naxis=2)
    # grab the image data
    hIFits = fits.open(df['file name'][i])
    data = hIFits[0].data
    print(np.min(data), np.max(data))

    # find the bar centre in pixels
    center = w.all_world2pix([[df['ra'][i], df['dec'][i]]], 0)[0]

    # Find the offset of the center point of the galaxy from that of the image
    offset = center - np.array(data.shape) / 2
    centred_image = transform.warp(
        data,
        transform.AffineTransform(translation=(offset)),
    )
    # create the centred image extent in arcseconds (for plotting)
    extent = np.repeat(
        np.array(centred_image.shape) * df['pixel size'][i] / 2,
        2
    ) * [-1, 1, -1, 1]
    # rotate so bar is in x axis (rotating about galaxy centre)
    rotated_image = transform.rotate(centred_image, df['rotation'][i])
    # rotate so bar in y-axis (for perpendicular cross section)
    rotated_image_perp = transform.rotate(centred_image, df['rotation'][i] - 90)

    y_density = 6

    X, Y = np.mgrid[
        0:rotated_image.shape[0]:rotated_image.shape[0]*1j,
        (rotated_image.shape[1]/2-1):(rotated_image.shape[1]/2+1):(y_density*1j)
    ]
    positions = np.vstack([X.ravel(), Y.ravel()])

    # average interpolated values along y-axis
    cross_section = map_coordinates(
        rotated_image,
        positions[::-1],
    ).reshape(-1, y_density).mean(axis=1)

    cross_section_perp = map_coordinates(
        rotated_image_perp,
        positions[::-1],
    ).reshape(-1, y_density).mean(axis=1)

    return {
        'file name': df['file name'][i],
        'pixel size': df['pixel size'][i],
        'rotation': df['rotation'][i],
        'data': centred_image,
        'extent': extent,
        'center': np.array((0, 0)),
        'cross_section': cross_section,
        'cross_section_perp': cross_section_perp,
        'bar length': df['bar length'][i] / df['pixel size'][i],
        'beam major': df['beam major'][i] / df['pixel size'][i],
        'beam minor': df['beam minor'][i] / df['pixel size'][i],
        'beam pa': df['beam pa'][i],
    }


for zoom in [False, True]:
    crop_size = 30 if zoom else 100
    fig, ax = plt.subplots(figsize=(18, 5*len(df)), ncols=3, nrows=len(df), dpi=100)
    for i in range(len(ax)):
        image_axis = ax[i][0]
        cross_section_axis = ax[i][1]
        cross_section_perp_axis = ax[i][2]

        # calculate required data for this galaxy
        g = get_galaxy(i)
        im = image_axis.imshow(
            g['data'], origin='lower',
            cmap='inferno', extent=g['extent']
        )
        bar_axis = shapely_rotate(
            box(
                g['center'][0] - 10000, g['center'][1] - 1*g['pixel size'],
                g['center'][0] + 10000, g['center'][1] + 1*g['pixel size'],
            ),
            g['rotation'],
        )
        bar = shapely_rotate(
            box(
                g['center'][0]-g['bar length'] / 2, g['center'][1]-1*g['pixel size'],
                g['center'][0]+g['bar length'] / 2, g['center'][1]+1*g['pixel size'],
            ),
            g['rotation'],
        )
        image_axis.add_patch(PolygonPatch(bar_axis, fc='C0', alpha=0.5))
        image_axis.add_patch(PolygonPatch(bar, fc='w', alpha=0.7))
        image_axis.add_patch(Ellipse(
            xy=(-60, -60),
            width=g['beam major'],
            height=g['beam minor'],
            angle=90 + g['beam pa'],
            edgecolor='w',
            facecolor='none',
        ))
        # crop to a useful size
        image_axis.set_xlim(-crop_size, crop_size)
        image_axis.set_ylim(-crop_size, crop_size)
        image_axis.set_xlabel('Arcseconds from centre')
        image_axis.set_ylabel('Arcseconds from centre')
        # add a colorbar
        c = plt.colorbar(im, ax=image_axis, fraction=0.046, pad=0.04)
        c.set_label('Intensity [JY.Hz]')

        # plot the cross section
        cross_section_axis.plot(
            np.linspace(g['extent'][0], g['extent'][1], g['data'].shape[0]),
            g['cross_section'],
        )
        cross_section_axis.add_line(Line2D(
            [g['center'][1]] * 2, [0, 100000],
            linestyle='--',
            c='k',
        ))

        # add a line denoting bar size
        barPosition = np.array([g['bar length']]*2) * [-0.5, 0.5] + g['center'][1]
        bar = Rectangle(
            [barPosition[0], 0],
            g['bar length'],
            10000,
            facecolor='k',
            alpha=0.1,
        )
        cross_section_axis.add_artist(bar)

        cross_section_axis.set_xlim(-crop_size, crop_size)
        cross_section_axis.set_ylim(bottom=0)
        cross_section_axis.set_xlabel('Arcseconds from centre')
        cross_section_perp_axis.set_ylabel('Intensity [JY.Hz]')

        # plot the perpendicular cross section
        cross_section_perp_axis.plot(
            np.linspace(g['extent'][0], g['extent'][1], g['data'].shape[0]),
            g['cross_section_perp'],
        )
        cross_section_perp_axis.set_xlim(-crop_size, crop_size)
        cross_section_perp_axis.set_ylim(bottom=0)
        cross_section_perp_axis.set_xlabel('Arcseconds from centre')
        cross_section_perp_axis.set_ylabel('Intensity [JY.Hz]')


    plt.tight_layout()
    plt.savefig(
        'output_images/cross-sections_bar{}'.format('_zoomed' if zoom else ''),
        bbox_inches='tight'
    )
