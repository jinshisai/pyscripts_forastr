### modules
import numpy as np
import sys
import subprocess
import matplotlib.pyplot as plt
import matplotlib.figure as figure
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import axes3d
from astropy.io import fits
import matplotlib.patches as patches
import mpl_toolkits.axes_grid1
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib.offsetbox import AnchoredText
import matplotlib as mpl
#plt.style.use('seaborn-dark')
#plt.style.use('ggplot')
#plt.style.use('seaborn-deep')
#plt.style.use('default')


### setting for figures
#mpl.use('Agg')
#mpl.rcParams['agg.path.chunksize'] = 100000
plt.rcParams['font.family'] ='Arial'    # font (Times New Roman, Helvetica, Arial)
plt.rcParams['xtick.direction'] = 'in'  # directions of x ticks ('in'), ('out') or ('inout')
plt.rcParams['ytick.direction'] = 'in'  # directions of y ticks ('in'), ('out') or ('inout')
#plt.rcParams['xtick.major.width'] = 1.0 # x ticks width
#plt.rcParams['ytick.major.width'] = 1.0 # y ticks width
plt.rcParams['font.size'] = 11           # fontsize
#plt.rcParams['axes.linewidth'] = 1.0    # edge linewidth



### parameters
formatlist = np.array(['eps','pdf','png','jpeg'])
clight     = 2.99792458e10 # light speed [cm s^-1]



### functions
def change_aspect_ratio(ax, ratio):
    '''
    This function change aspect ratio of figure.
    Parameters:
        ax: ax (matplotlit.pyplot.subplots())
            Axes object
        ratio: float or int
            relative x axis width compared to y axis width.
    '''
    aspect = (1/ratio) *(ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
    aspect = np.abs(aspect)
    aspect = float(aspect)
    ax.set_aspect(aspect)


def Idistmap(fitsdata, outname=None, imscale=None, outformat='eps', color=True, cmap='Greys',
             colorbar=False, cbaroptions=np.array(['right','5%','0%','Jy/beam']), vmin=None,vmax=None,
             contour=True, clevels=np.array([0.15, 0.3, 0.45, 0.6, 0.75, 0.9]), ccolor='k',
             xticks=np.empty, yticks=np.empty, relativecoords=True, csize=9, scalebar=np.empty(0),
             cstar=True, prop_star=np.array(['1','0.5','white']), locsym=0.1, bcolor='k',figsize=(11.69,8.27)):
    '''
    Make a figure from single image.

    Args
    fitsdata (fits file): Input fitsdata. It must be image data having 3 or 4 axes.
    outname (str): Output file name. Not including file extension.
    outformat (str): Extension of the output file. Default is eps.
    imscale (ndarray): Image scale [arcsec]. Input as np.array([xmin,xmax,ymin,ymax]).
    color (bool): If True, image will be described in color scale.
    cmap: Choose colortype of color scale
    colorbar (bool): If True, color bar will be showen. Default False.
    cbaroptions: Set colorbar options.
    vmin: Cut off of color scale. Abusolute value.
    contour (bool): If True, contour will be drawn.
    clevels (ndarray): Set contour levels. Abusolute value.
    ccolor: Set contour color.
    xticks, yticks: Optional setting. If input ndarray, set xticks and yticsk as well as input.
    relativecoords (bool): If True, the coordinate is shown in relativecoordinate. Default True.
    csize: Font size.
    scalebar: Optional setting. Input ndarray([barx, bary, barlength, textx, texty, text ]).
     barx and bary are the position where scalebar will be putted. [arcsec].
     barlength is the length of the scalebar in the figure, so in arcsec.
     textx and texty are the position where a label of scalebar will be putted. [arcsec].
     text is a text which represents the scale of the scalebar.
    cstar (bool): If True, central star position will be marked by cross.
    locsym: Set locations of beamsize.
    '''

    ### setting output file name & format
    if (outformat == formatlist).any():
        outname = outname + '.' + outformat
    else:
        print 'ERROR\tsingleim_to_fig: Outformat is wrong.'
        return


    ### reading fits files
    data, header = fits.getdata(fitsdata,header=True)

    # reading header info.
    xlabel   = header['CTYPE1']
    ylabel   = header['CTYPE2']
    try:
        restfreq = header['RESTFRQ'] # Hz
    except:
        restfreq = header['RESTFREQ'] # Hz
    refval_x = header['CRVAL1']*60.*60.
    refval_y = header['CRVAL2']*60.*60.
    refpix_x = int(header['CRPIX1'])
    refpix_y = int(header['CRPIX2'])
    del_x    = header['CDELT1']*60.*60. # deg --> arcsec
    del_y    = header['CDELT2']*60.*60.
    nx       = header['NAXIS1']
    ny       = header['NAXIS2']
    bmaj     = header['BMAJ']*60.*60.
    bmin     = header['BMIN']*60.*60.
    bpa      = header['BPA']  # [deg]
    unit     = header['BUNIT']
    print 'x, y axes are ', xlabel, ' and ', ylabel

    # setting axes in relative coordinate
    if relativecoords:
        refval_x, refval_y = [0,0]
        xlabel = 'RA offset (arcsec; J2000)'
        ylabel = 'Dec offset (arcsec; J2000)'
    else:
        pass
    xmin = refval_x + (1 - refpix_x)*del_x - 0.5*del_x # refpix_x is not 0 start
    xmax = refval_x + (nx - refpix_x)*del_x + 0.5*del_x
    ymin = refval_y + (1 - refpix_y)*del_y - 0.5*del_y
    ymax = refval_y + (ny - refpix_y)*del_y + 0.5*del_y

    # check data axes
    if len(data.shape) == 2:
        pass
    elif len(data.shape) == 3:
        data = data[0,:,:]
    elif len(data.shape) == 4:
        data = data[0,0,:,:]
    else:
        print 'Error\tsingleim_to_fig: Input fits size is not corrected.\
         It is allowed only to have 3 or 4 axes. Check the shape of the fits file.'
        return


    ### ploting
    # setting figure
    fig = plt.figure(figsize=figsize)
    ax  = fig.add_subplot(111)
    plt.rcParams['font.size'] = csize
    #print 'start plot'

    # showing in color scale
    if color:
        # color image
        imcolor = ax.imshow(data, cmap=cmap, origin='lower', extent=(xmin,xmax,ymin,ymax),vmin=vmin,vmax=vmax)
        # color bar
        print 'plot color'
        if colorbar:
            divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
            cax     = divider.append_axes('right', '5%', pad='0%')
            cbar    = fig.colorbar(imcolor, cax = cax)
            cbar.set_label('Jy/beam')

    if contour:
        imcont02 = ax.contour(data, colors=ccolor, origin='lower',extent=(xmin,xmax,ymin,ymax), levels=clevels,linewidths=1)
        print 'plot contour'

    # set axes
    figxmin, figxmax, figymin, figymax = imscale
    ax.set_xlim(figxmax,figxmin)
    ax.set_ylim(figymin,figymax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xticks != np.empty and yticks != np.empty:
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
    else:
        pass
    ax.set_aspect(1)
    ax.tick_params(which='both', direction='in',bottom=True, top=True, left=True, right=True)

    # plot beam size
    beam    = patches.Ellipse(xy=(figxmax-locsym*figxmax, figymin-locsym*figymin), width=bmin, height=bmaj, fc=bcolor, angle=-bpa)
    ax.add_patch(beam)

    # central star position
    if cstar:
        ll,lw, cl = prop_star
        ll = float(ll)
        lw = float(lw)

        cross01 = patches.Arc(xy=(refval_x,refval_y), width=ll, height=0.001, lw=lw, color=cl,zorder=11)
        cross02 = patches.Arc(xy=(refval_x,refval_y), width=0.001, height=ll, lw=lw, color=cl,zorder=12)
        ax.add_patch(cross01)
        ax.add_patch(cross02)

    # scale bar
    if len(scalebar) == 0:
        pass
    elif len(scalebar) == 7:
        barx, bary, barlength, textx, texty, text, colors = scalebar

        barx      = float(barx)
        bary      = float(bary)
        barlength = float(barlength)
        textx     = float(textx)
        texty     = float(texty)

        scale   = patches.Arc(xy=(barx,bary), width=barlength, height=0., lw=2, color=colors,zorder=10)
        ax.add_patch(scale)
        ax.text(textx,texty,text,color=colors)
    else:
        print 'scalebar must be 7 elements. Check scalebar.'

    fig.savefig(outname, transparent = True)

    return ax



### moment maps
def multiIdistmap(fitsdata, outname=None, imscale=None, outformat='eps', cmap='Greys',
             colorbar=False, cbaroptions=np.array(['right','5%','0%','Jy/beam']), vmin=None, vmax=None,
             clevels=np.array([0.15, 0.3, 0.45, 0.6, 0.75, 0.9]), ccolor='k',
             xticks=np.empty, yticks=np.empty, relativecoords=True, csize=9, scalebar=np.empty(0),
             cstar=True, prop_star=np.array(['1','0.5','white']), locsym=0.1, bcolor='k',logscale=False):
    '''
    Make a figure from multi images.

    Args
    fitsdata (fits file): Input fitsdata. It must be image data having 3 or 4 axes.
    outname (str): Output file name. Not including file extension.
    outformat (str): Extension of the output file. Default is eps.
    imscale (ndarray): Image scale [arcsec]. Input as np.array([xmin,xmax,ymin,ymax]).
    color (bool): If True, image will be described in color scale.
    cmap: Choose colortype of color scale
    colorbar (bool): If True, color bar will be showen. Default False.
    cbaroptions: Set colorbar options.
    vmin: Cut off of color scale. Abusolute value.
    contour (bool): If True, contour will be drawn.
    clevels (ndarray): Set contour levels. Abusolute value.
    ccolor: Set contour color.
    xticks, yticks: Optional setting. If input ndarray, set xticks and yticsk as well as input.
    relativecoords (bool): If True, the coordinate is shown in relativecoordinate. Default True.
    csize: Font size.
    scalebar: Optional setting. Input ndarray([barx, bary, barlength, textx, texty, text ]).
     barx and bary are the position where scalebar will be putted. [arcsec].
     barlength is the length of the scalebar in the figure, so in arcsec.
     textx and texty are the position where a label of scalebar will be putted. [arcsec].
     text is a text which represents the scale of the scalebar.
    cstar (bool): If True, central star position will be marked by cross.
    locsym: Set locations of beamsize.
    '''

    ### setting output file name & format
    if (outformat == formatlist).any():
        outname = outname + '.' + outformat
    else:
        print 'ERROR\tsingleim_to_fig: Outformat is wrong.'
        return


    ### setting initial parameters
    # for data
    data   = None
    data02 = None
    data03 = None
    data04 = None

    # for header
    header = None
    hd02   = None
    hd03   = None
    hd04   = None

    ### reading fits files
    if type(fitsdata) == tuple:
        nfits = len(fitsdata)
        if nfits == 1:
            # read fits
            data, header = fits.getdata(fitsdata[0],header=True)

            # data shape
            if len(data.shape) == 2:
                pass
            elif len(data.shape) == 3:
                data = data[0,:,:]
            elif len(data.shape) == 4:
                data = data[0,0,:,:]
            else:
                print 'Error\tmultiIdistmap: Input fits size is not corrected.\
                 It is allowed only to have 3 or 4 axes. Check the shape of the fits file.'
                return
        elif nfits == 2:
            # read fits
            data, header = fits.getdata(fitsdata[0],header=True)
            data02, hd02 = fits.getdata(fitsdata[1],header=True)
            #print len(data.shape)

            # data shape
            if len(data.shape) == 2:
                pass
            elif len(data.shape) == 3:
                data = data[0,:,:]
            elif len(data.shape) == 4:
                data = data[0,0,:,:]
            else:
                print 'Error\tmultiIdistmap: Input fits size is not corrected.\
                 It is allowed only to have 3 or 4 axes. Check the shape of the fits file.'
                return

            # data02 shape
            if len(data02.shape) == 2:
                pass
            elif len(data02.shape) == 3:
                data02 = data02[0,:,:]
            elif len(data02.shape) == 4:
                data02 = data02[0,0,:,:]
            else:
                #print 'here'
                print 'Error\tmultiIdistmap: Input fits size is not corrected.\
                 It is allowed only to have 3 or 4 axes. Check the shape of the fits file.'
                return
        else:
            print 'Cannot deal with images more than 3. Now under development.'
            print 'Only use the first and second one, and ignore rest.'
            data, header = fits.getdata(fitsdata[0],header=True)
            data02, hd02 = fits.getdata(fitsdata[1],header=True)

            # data shape
            if len(data.shape) == 2:
                pass
            elif len(data.shape) == 3:
                data = data[0,:,:]
            elif len(data.shape) == 4:
                data = data[0,0,:,:]
            else:
                print 'Error\tmultiIdistmap: Input fits size is not corrected.\
                 It is allowed only to have 3 or 4 axes. Check the shape of the fits file.'
                return

            # data02 shape
            if len(data02.shape) == 2:
                pass
            elif len(data02.shape) == 3:
                data02 = data02[0,:,:]
            elif len(data.shape) == 4:
                data02 = data02[0,0,:,:]
            else:
                print 'Error\tmultiIdistmap: Input fits size is not corrected.\
                 It is allowed only to have 3 or 4 axes. Check the shape of the fits file.'
                return
    else:
        print 'Input type is wrong. Input type must be tuple.'
        return

    # reading header info.
    xlabel   = header['CTYPE1']
    ylabel   = header['CTYPE2']
    try:
        restfreq = header['RESTFRQ'] # Hz
    except:
        restfreq = header['RESTFREQ'] # Hz
    refval_x = header['CRVAL1']*60.*60.
    refval_y = header['CRVAL2']*60.*60.
    refpix_x = int(header['CRPIX1'])
    refpix_y = int(header['CRPIX2'])
    del_x    = header['CDELT1']*60.*60. # deg --> arcsec
    del_y    = header['CDELT2']*60.*60.
    nx       = header['NAXIS1']
    ny       = header['NAXIS2']
    bmaj     = header['BMAJ']*60.*60.
    bmin     = header['BMIN']*60.*60.
    bpa      = header['BPA']  # [deg]
    unit     = header['BUNIT']
    print 'x, y axes are ', xlabel, ' and ', ylabel

    # setting axes in relative coordinate
    if relativecoords:
        refval_x, refval_y = [0,0]
        xlabel = 'RA offset (arcsec; J2000)'
        ylabel = 'Dec offset (arcsec; J2000)'
    else:
        pass
    xmin = refval_x + (1 - refpix_x)*del_x - 0.5*del_x
    xmax = refval_x + (nx - refpix_x)*del_x + 0.5*del_x
    ymin = refval_y + (1 - refpix_y)*del_y - 0.5*del_y
    ymax = refval_y + (ny - refpix_y)*del_y + 0.5*del_y

    # set colorscale
    if vmax:
        pass
    else:
        vmax = np.nanmax(data)

    # logscale
    if logscale:
        norm = mpl.colors.LogNorm(vmin=vmin,vmax=vmax)
    else:
        norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)

    ### ploting
    # setting figure
    fig = plt.figure(figsize=(11.69,8.27))
    ax  = fig.add_subplot(111)
    plt.rcParams['font.size'] = csize

    # showing in color scale
    if data is not None:
        # color image
        imcolor = ax.imshow(data, cmap=cmap, origin='lower', extent=(xmin,xmax,ymin,ymax),norm=norm)
        # color bar
        if colorbar:
            divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
            cax     = divider.append_axes('right', '5%', pad='0%')
            cbar    = fig.colorbar(imcolor, cax = cax)
            cbar.set_label('Jy/beam')

    if data02 is not None:
        imcont02 = ax.contour(data02, colors=ccolor, origin='lower',extent=(xmin,xmax,ymin,ymax), levels=clevels,linewidths=1)

    # set axes
    figxmin, figxmax, figymin, figymax = imscale
    ax.set_xlim(figxmax,figxmin)
    ax.set_ylim(figymin,figymax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xticks != np.empty and yticks != np.empty:
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
    else:
        pass
    ax.set_aspect(1)
    ax.tick_params(which='both', direction='in',bottom=True, top=True, left=True, right=True)

    # plot beam size
    beam    = patches.Ellipse(xy=(figxmax-locsym*figxmax, figymin-locsym*figymin), width=bmin, height=bmaj, fc=bcolor, angle=-bpa)
    ax.add_patch(beam)

    # central star position
    if cstar:
        ll,lw, cl = prop_star
        ll = float(ll)
        lw = float(lw)

        cross01 = patches.Arc(xy=(refval_x,refval_y), width=ll, height=0.001, lw=lw, color=cl,zorder=11)
        cross02 = patches.Arc(xy=(refval_x,refval_y), width=0.001, height=ll, lw=lw, color=cl,zorder=12)
        ax.add_patch(cross01)
        ax.add_patch(cross02)

    # scale bar
    if len(scalebar) == 0:
        pass
    elif len(scalebar) == 7:
        barx, bary, barlength, textx, texty, text, colors = scalebar

        barx      = float(barx)
        bary      = float(bary)
        barlength = float(barlength)
        textx     = float(textx)
        texty     = float(texty)

        scale   = patches.Arc(xy=(barx,bary), width=barlength, height=0., lw=2, color=colors,zorder=10)
        ax.add_patch(scale)
        ax.text(textx,texty,text,color=colors)
    else:
        print 'scalebar must consist of 7 elements. Check scalebar.'

    fig.savefig(outname, transparent = True)

    return ax


### channel map
def channelmap(fitsdata, outname=None, outformat='eps', imscale=None, color=False,cbaron=False,cmap='Greys', vmin=None, vmax=None,
                contour=True, clevels=np.array([0.15, 0.3, 0.45, 0.6, 0.75, 0.9]), ccolor='k',
                nrow=5, ncol=5,velmin=None, velmax=None, nskip=1,
                xticks=np.empty, yticks=np.empty, relativecoords=True, vsys=None, csize=9, scalebar=np.empty,
                cstar=True, prop_star=np.array(['1','0.5','red']), locsym=0.2, logscale=False, tickcolor='k',axiscolor='k',
                labelcolor='k',cbarlabel=None, txtcolor='k', bcolor='k', figsize=(8.27,11.69)):
    '''
    Make channel maps from single image.

    Args
    fitsdata: Input fitsdata. It must be an image cube having 3 or 4 axes.
    outname: Output file name. Not including file extension.
    outformat: Extension of the output file. Default is eps.
    imscale: scale to be shown (arcsec). It must be given as [xmin, xmax, ymin, ymax].
    color (bool): If True, images will be shown in colorscale. Default is False.
        cmap: color of the colorscale.
        vmin: Minimum value of colorscale. Default is None.
        vmax: Maximum value of colorscale. Default is the maximum value of the image cube.
        logscale (bool): If True, the color will be shown in logscale.
    contour (bool): If True, images will be shown with contour. Default is True.
        clevels (ndarray): Contour levels. Input will be treated as absolute values.
        ccolor: color of contour.
    nrow, ncol: the number of row and column of the channel map.
    relativecoords (bool): If True, the channel map will be produced in relative coordinate. Abusolute coordinate mode is (probably) coming soon.
    velmin, velmax: Minimum and maximum velocity to be shown.
    vsys: Systemic velicity [km/s]. If no input value, velocities will be described in LSRK.
    csize: Caracter size. Default is 9.
    cstar: If True, a central star or the center of an image will be shown as a cross.
    locsym: A factor to decide locations of symbols (beam and velocity label). It must be 0 - 1.
    tickcolor, axiscolor, labelcolor, txtcolor: Colors for each stuffs.
    scalebar (array): If it is given, scalebar will be drawn. It must be given as [barx, bary, bar length, textx, texty, text].
                       Barx, bary, textx, and texty are locations of a scalebar and a text in arcsec.
    nskip: the number of channel skipped
    '''

    ### setting output file name & format
    nmap = 1
    if (outformat == formatlist).any():
        outfile = outname + '_nmap{0:02d}'.format(nmap) + '.' + outformat
    else:
        print 'ERROR\tsingleim_to_fig: Outformat is wrong.'
        return


    ### reading fits files
    data, header = fits.getdata(fitsdata,header=True)

    # reading header info.
    xlabel   = header['CTYPE1']
    ylabel   = header['CTYPE2']
    try:
        restfreq = header['RESTFRQ'] # Hz
    except:
        restfreq = header['RESTFREQ'] # Hz
    refval_x = header['CRVAL1']*60.*60. # deg --> arcsec
    refval_y = header['CRVAL2']*60.*60.
    refval_v = header['CRVAL3']
    refpix_x = int(header['CRPIX1'])
    refpix_y = int(header['CRPIX2'])
    refpix_v = int(header['CRPIX3'])
    del_x    = header['CDELT1']*60.*60. # deg --> arcsec
    del_y    = header['CDELT2']*60.*60.
    del_v    = header['CDELT3']
    nx       = header['NAXIS1']
    ny       = header['NAXIS2']
    nchan    = header['NAXIS3']
    bmaj     = header['BMAJ']*60.*60.
    bmin     = header['BMIN']*60.*60.
    bpa      = header['BPA']  # [deg]
    unit     = header['BUNIT']
    print 'x, y axes are ', xlabel, ' and ', ylabel

    # frequency --> velocity
    print 'The third axis is [FREQUENCY]'
    print 'Convert frequency to velocity'
    del_v    = - del_v*clight/restfreq       # delf --> delv [cm/s]
    del_v    = del_v*1.e-5                   # cm/s --> km/s
    refval_v = clight*(1.-refval_v/restfreq) # radio velocity c*(1-f/f0) [cm/s]
    refval_v = refval_v*1.e-5                # cm/s --> km/s
    #print del_v

    # setting axes in relative coordinate
    if relativecoords:
        refval_x, refval_y = [0,0]
        xlabel = 'RA offset (arcsec; J2000)'
        ylabel = 'Dec offset (arcsec; J2000)'
    else:
        pass
    xmin = refval_x + (1 - refpix_x)*del_x - 0.5*del_x
    xmax = refval_x + (nx - refpix_x)*del_x + 0.5*del_x
    ymin = refval_y + (1 - refpix_y)*del_y - 0.5*del_y
    ymax = refval_y + (ny - refpix_y)*del_y + 0.5*del_y

    # setting velocity axis in relative velocity
    if vsys:
        refval_v = refval_v - vsys
        #print refval_v


    # check data axes
    if len(data.shape) == 4:
        pass
    else:
        print 'Error\tsingleim_to_fig: Input fits size is not corrected.\
         It is allowed only to have 3 or 4 axes. Check the shape of the fits file.'
        return

    # set colorscale
    if vmax:
        pass
    else:
        vmax = np.nanmax(data[0,:,:,:])

    # logscale
    if logscale:
        norm = mpl.colors.LogNorm(vmin=vmin,vmax=vmax)
    else:
        norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)

    ### ploting
    # setting figure
    plt.rcParams['font.size'] = csize
    fig = plt.figure(figsize=figsize)
    # setting colorbar
    if cbaron:
        cbar_mode = 'single'
    else:
        cbar_mode= None

    # setting grid
    grid = ImageGrid(fig, rect=111, nrows_ncols=(nrow,ncol),
            axes_pad=0,share_all=True,cbar_mode=cbar_mode)

    # setting parameters used to plot
    i, j, gridi = [0,0,0]
    gridimax    = nrow*ncol-1
    figxmin, figxmax, figymin, figymax = imscale
    xscale = np.abs((figxmax - figxmin)*0.5)
    yscale = np.abs((figymax - figymin)*0.5)

    nroop = nchan//nskip + 1
    for k in xrange(0,nroop):
        ichan = k*nskip
        if ichan >= nchan:
            continue
        # select channel
        dataim = data[0,ichan,:,:]

        # velocity at nchan
        vnchan = refval_v + (ichan + 1 - refpix_v)*del_v

        # check whether vnchan in setted velocity range
        if velmax:
            if vnchan < velmin or vnchan > velmax:
                continue
        elif velmin:
            if vnchan < velmin:
                continue
        else:
            pass

        # each plot
        #ax  = fig.add_subplot(gs[i,j])
        ax = grid[gridi]
        print 'channel ', '%s'%ichan, ', velocity: ', '%4.2f'%vnchan, ' km/s'

        # showing in color scale
        if color:
            imcolor = ax.imshow(dataim, cmap=cmap, origin='lower', extent=(xmin,xmax,ymin,ymax),norm=norm)

        if contour:
            imcont  = ax.contour(dataim, colors=ccolor, origin='lower',extent=(xmin,xmax,ymin,ymax), levels=clevels, linewidths=0.5)

        # set axes
        ax.set_xlim(figxmax,figxmin)
        ax.set_ylim(figymin,figymax)
        ax.spines["bottom"].set_color(axiscolor)
        ax.spines["top"].set_color(axiscolor)
        ax.spines["left"].set_color(axiscolor)
        ax.spines["right"].set_color(axiscolor)
        if xticks != np.empty and yticks != np.empty:
            ax.set_xticks(xticks)
            ax.set_yticks(yticks)
        else:
            pass
        ax.set_aspect(1)
        ax.tick_params(which='both', direction='in',bottom=True, top=True, left=True, right=True,colors=tickcolor)

        # velocity
        #vlabel = AnchoredText('%03.2f'%vnchan,loc=2,frameon=False)
        #ax.add_artist(vlabel)
        vlabel = '%03.2f'%vnchan
        ax.text(figxmax-locsym*3.*xscale, figymax-locsym*3.*yscale,vlabel,color=txtcolor,size=csize)

        # only on the bottom corner pannel
        if i == nrow-1 and j == 0:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.xaxis.label.set_color(labelcolor)
            ax.yaxis.label.set_color(labelcolor)

            # plot beam size
            beam = patches.Ellipse(xy=(figxmax-locsym*xscale, figymin+locsym*yscale), width=bmin, height=bmaj, fc=bcolor, angle=-bpa)
            ax.add_patch(beam)
            #beam = AnchoredEllipse(ax.transData,width=bmin, height=bmaj, angle=-bpa,loc=3, pad=0.5, borderpad=0.4, frameon=False)
            #ax.add_artist(beam)

            # scale bar
            if scalebar is np.empty:
                pass
            else:
                barx, bary, barlength = scalebar[0,1,2]
                textx, texty, text    = scalebar[3,4,5]

                barx      = float(barx)
                bary      = float(bary)
                barlength = float(barlength)
                textx     = float(textx)
                texty     = float(texty)

                scale   = patches.Arc(xy=(barx,bary), width=barlength, height=0., lw=2, color='k',zorder=10)
                ax.add_patch(scale)
                ax.text(textx,texty,text)
        else:
            ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)

        # central star position
        if cstar:
            ll,lw, cl = prop_star
            ll = float(ll)
            lw = float(lw)

            cross01 = patches.Arc(xy=(refval_x,refval_y), width=ll, height=0.001, lw=lw, color=cl, zorder=11)
            cross02 = patches.Arc(xy=(refval_x,refval_y), width=0.001, height=ll, lw=lw, color=cl, zorder=12)
            ax.add_patch(cross01)
            ax.add_patch(cross02)

        # counts
        j     = j + 1
        gridi = gridi+1

        if j == ncol:
            j = 0
            i = i + 1

        if i == nrow:
            #gs.tight_layout(fig,h_pad=0,w_pad=0)
            #plt.subplots_adjust(wspace=0., hspace=0.)
            if cbaron:
                # With cbar_mode="single", cax attribute of all axes are identical.
                cbar = ax.cax.colorbar(imcolor)
                ax.cax.toggle_label(True)
                cbar.ax.yaxis.set_tick_params(color=tickcolor) # tick color
                cbar.ax.spines["bottom"].set_color(axiscolor)  # axes color
                cbar.ax.spines["top"].set_color(axiscolor)
                cbar.ax.spines["left"].set_color(axiscolor)
                cbar.ax.spines["right"].set_color(axiscolor)
                if cbarlabel:
                    cbar.ax.set_ylabel(cbarlabel, color=labelcolor) # label

            fig.savefig(outfile, transparent = True)
            fig.clf()
            nmap      = nmap+1
            outfile   = outname + '_nmap{0:02d}'.format(nmap) + '.' + outformat
            fig       = plt.figure(figsize=figsize)
            grid      = ImageGrid(fig, rect=111, nrows_ncols=(nrow,ncol),
                axes_pad=0,share_all=True,cbar_mode=cbar_mode)
            i,j,gridi = [0,0,0]

    if color:
        # With cbar_mode="single", cax attribute of all axes are identical.
        cbar = ax.cax.colorbar(imcolor)
        ax.cax.toggle_label(True)
        cbar.ax.yaxis.set_tick_params(color=tickcolor) # tick color
        cbar.ax.spines["bottom"].set_color(axiscolor)  # axes color
        cbar.ax.spines["top"].set_color(axiscolor)
        cbar.ax.spines["left"].set_color(axiscolor)
        cbar.ax.spines["right"].set_color(axiscolor)
        if cbarlabel:
            cbar.ax.set_ylabel(cbarlabel,color=labelcolor) # label

    if gridi != gridimax+1 and gridi != 0:
        while gridi != gridimax+1:
            #print gridi
            ax = grid[gridi]
            ax.spines["right"].set_color("none")  # right
            ax.spines["left"].set_color("none")   # left
            ax.spines["top"].set_color("none")    # top
            ax.spines["bottom"].set_color("none") # bottom
            ax.axis('off')
            gridi = gridi+1

        #gs.tight_layout(fig,h_pad=0,w_pad=0)
        #plt.subplots_adjust(wspace=0, hspace=0)
        fig.savefig(outfile, transparent = True)

    # convert N maps into one map
    # ImageMagick is needed
    try:
        outfile = outname + '_nmap*' + '.' + outformat
        sumfile = outname + '.pdf'
        cmdop   = ' -density 600x600 '
        cmd = 'convert' + cmdop + outfile + ' ' + sumfile
        subprocess.call(cmd,shell=True)
    except:
        print 'Cannot convert N maps into one pdf file.'
        print 'Install ImageMagick if you want.'

    return ax




### channel map
def mltichannelmap(fits01, fits02, outname=None, outformat='eps', imscale=None, cmap='Greys', vmin=None, vmax=None,
                clevels=np.array([0.15, 0.3, 0.45, 0.6, 0.75, 0.9]), ccolor='k',
                nrow=5, ncol=5,velmin=None, velmax=None, nskip=1,
                xticks=np.empty, yticks=np.empty, relativecoords=True, vsys=None, csize=9, scalebar=np.empty,
                cstar=True, locsym=0.2, logscale=False, tickcolor='k',axiscolor='k',labelcolor='k',cbarlabel=None, txtcolor='k',
                bcolor='k', figsize=(8.27,11.69)):
    '''
    Make channel maps from images.

    Args
    fitsdata: Input fitsdata. It must be an image cube having 3 or 4 axes.
    imtype (str): Choose which image will be drawn with contour or color. 'cl' means color, 'cn' means contour.
    outname: Output file name. Not including file extension.
    outformat: Extension of the output file. Default is eps.
    imscale: scale to be shown (arcsec). It must be given as [xmin, xmax, ymin, ymax].
    color (bool): If True, images will be shown in colorscale. Default is False.
        cmap: color of the colorscale.
        vmin: Minimum value of colorscale. Default is None.
        vmax: Maximum value of colorscale. Default is the maximum value of the image cube.
        logscale (bool): If True, the color will be shown in logscale.
    contour (bool): If True, images will be shown with contour. Default is True.
        clevels (ndarray): Contour levels. Input will be treated as absolute values.
        ccolor: color of contour.
    nrow, ncol: the number of row and column of the channel map.
    relativecoords (bool): If True, the channel map will be produced in relative coordinate. Abusolute coordinate mode is (probably) coming soon.
    velmin, velmax: Minimum and maximum velocity to be shown.
    vsys: Systemic velicity [km/s]. If no input value, velocities will be described in LSRK.
    csize: Caracter size. Default is 9.
    cstar: If True, a central star or the center of an image will be shown as a cross.
    locsym: A factor to decide locations of symbols (beam and velocity label). It must be 0 - 1.
    tickcolor, axiscolor, labelcolor, txtcolor: Colors for each stuffs.
    scalebar (array): If it is given, scalebar will be drawn. It must be given as [barx, bary, bar length, textx, texty, text].
                       Barx, bary, textx, and texty are locations of a scalebar and a text in arcsec.
    '''

    ### setting output file name & format
    nmap = 1
    if (outformat == formatlist).any():
        outfile = outname + '_nmap{0:02d}'.format(nmap) + '.' + outformat
    else:
        print 'ERROR\tsingleim_to_fig: Outformat is wrong.'
        return


    ### reading fits files
    # color image
    data01, hd01 = fits.getdata(fits01,header=True)
    # contour image
    data02, hd02 = fits.getdata(fits02,header=True)

    # reading hd01 info.
    xlabel   = hd01['CTYPE1']
    ylabel   = hd01['CTYPE2']
    try:
        restfreq = hd01['RESTFRQ'] # Hz
    except:
        restfreq = hd01['RESTFREQ'] # Hz
    refval_x = hd01['CRVAL1']*60.*60. # deg --> arcsec
    refval_y = hd01['CRVAL2']*60.*60.
    refval_v = hd01['CRVAL3']
    refpix_x = int(hd01['CRPIX1'])
    refpix_y = int(hd01['CRPIX2'])
    refpix_v = int(hd01['CRPIX3'])
    del_x    = hd01['CDELT1']*60.*60. # deg --> arcsec
    del_y    = hd01['CDELT2']*60.*60.
    del_v    = hd01['CDELT3']
    nx       = hd01['NAXIS1']
    ny       = hd01['NAXIS2']
    nchan    = hd01['NAXIS3']
    bmaj     = hd01['BMAJ']*60.*60.
    bmin     = hd01['BMIN']*60.*60.
    bpa      = hd01['BPA']  # [deg]
    unit     = hd01['BUNIT']
    print 'x, y axes are ', xlabel, ' and ', ylabel

    # frequency --> velocity
    print 'The third axis is [FREQUENCY]'
    print 'Convert frequency to velocity'
    del_v    = - del_v*clight/restfreq       # delf --> delv [cm/s]
    del_v    = del_v*1.e-5                   # cm/s --> km/s
    refval_v = clight*(1.-refval_v/restfreq) # radio velocity c*(1-f/f0) [cm/s]
    refval_v = refval_v*1.e-5                # cm/s --> km/s
    #print del_v

    # setting axes in relative coordinate
    if relativecoords:
        refval_x, refval_y = [0,0]
        xlabel = 'RA offset (arcsec; J2000)'
        ylabel = 'Dec offset (arcsec; J2000)'
    else:
        pass
    xmin = refval_x + (1 - refpix_x)*del_x - 0.5*del_x
    xmax = refval_x + (nx - refpix_x)*del_x + 0.5*del_x
    ymin = refval_y + (1 - refpix_y)*del_y - 0.5*del_y
    ymax = refval_y + (ny - refpix_y)*del_y + 0.5*del_y

    # setting velocity axis in relative velocity
    if vsys:
        refval_v = refval_v - vsys
        #print refval_v


    # check data axes
    if len(data01.shape) == 4:
        pass
    else:
        print 'Error\tsingleim_to_fig: Input fits size is not corrected.\
         It is allowed only to have 4 axes. Check the shape of the fits file.'
        return

    # set colorscale
    if vmax:
        pass
    else:
        vmax = np.nanmax(data01[0,:,:,:])

    # logscale
    if logscale:
        norm = mpl.colors.LogNorm(vmin=vmin,vmax=vmax)
    else:
        norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)

    ### ploting
    # setting figure
    plt.rcParams['font.size'] = csize
    fig = plt.figure(figsize=figsize)
    # setting colorbar
    #if color:
    cbar_mode = 'single'
    #else:
        #cbar_mode= None

    # setting grid
    grid = ImageGrid(fig, rect=111, nrows_ncols=(nrow,ncol),
            axes_pad=0,share_all=True,cbar_mode=cbar_mode)

    # setting parameters used to plot
    i, j, gridi = [0,0,0]
    gridimax    = nrow*ncol-1
    figxmin, figxmax, figymin, figymax = imscale


    nroop = nchan//nskip + 1
    for k in xrange(0,nroop):
        # select channel
        ichan   = k*nskip
        if ichan >= nchan:
            break

        color   = data01[0,ichan,:,:]
        contour = data02[0,ichan,:,:]

        # velocity at nchan
        vnchan = refval_v + (ichan + 1 - refpix_v)*del_v

        # check whether vnchan in setted velocity range
        if velmax:
            if vnchan < velmin or vnchan > velmax:
                continue
        elif velmin:
            if vnchan < velmin:
                continue
        else:
            pass

        # each plot
        #ax  = fig.add_subplot(gs[i,j])
        ax = grid[gridi]
        print 'channel ', '%s'%ichan, ', velocity: ', '%4.2f'%vnchan, ' km/s'

        # showing in color scale
        imcolor = ax.imshow(color, cmap=cmap, origin='lower', extent=(xmin,xmax,ymin,ymax),norm=norm)

        imcont  = ax.contour(contour, colors=ccolor, origin='lower',extent=(xmin,xmax,ymin,ymax), levels=clevels, linewidths=0.5)

        # set axes
        ax.set_xlim(figxmax,figxmin)
        ax.set_ylim(figymin,figymax)
        ax.spines["bottom"].set_color(axiscolor)
        ax.spines["top"].set_color(axiscolor)
        ax.spines["left"].set_color(axiscolor)
        ax.spines["right"].set_color(axiscolor)
        if xticks != np.empty and yticks != np.empty:
            ax.set_xticks(xticks)
            ax.set_yticks(yticks)
        else:
            pass
        ax.set_aspect(1)
        ax.tick_params(which='both', direction='in',bottom=True, top=True, left=True, right=True,colors=tickcolor)

        # velocity
        #vlabel = AnchoredText('%03.2f'%vnchan,loc=2,frameon=False)
        #ax.add_artist(vlabel)
        vlabel = '%03.2f'%vnchan
        ax.text(figxmax-locsym*figxmax, figymax-locsym*figymax,vlabel,color=txtcolor,size=csize)

        # only on the bottom corner pannel
        if i == nrow-1 and j == 0:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.xaxis.label.set_color(labelcolor)
            ax.yaxis.label.set_color(labelcolor)

            # plot beam size
            beam = patches.Ellipse(xy=(figxmax-locsym*figxmax, figymin-locsym*figymin), width=bmin, height=bmaj, fc=bcolor, angle=-bpa)
            ax.add_patch(beam)
            #beam = AnchoredEllipse(ax.transData,width=bmin, height=bmaj, angle=-bpa,loc=3, pad=0.5, borderpad=0.4, frameon=False)
            #ax.add_artist(beam)

            # scale bar
            if scalebar is np.empty:
                pass
            else:
                barx, bary, barlength = scalebar[0,1,2]
                textx, texty, text    = scalebar[3,4,5]

                barx      = float(barx)
                bary      = float(bary)
                barlength = float(barlength)
                textx     = float(textx)
                texty     = float(texty)

                scale   = patches.Arc(xy=(barx,bary), width=barlength, height=0., lw=2, color='k',zorder=10)
                ax.add_patch(scale)
                ax.text(textx,texty,text)
        else:
            ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)

        # central star position
        if cstar:
            cross01 = patches.Arc(xy=(refval_x,refval_y), width=1, height=0.001, lw=0.5, color='red',zorder=11)
            cross02 = patches.Arc(xy=(refval_x,refval_y), width=0.001, height=1, lw=0.5, color='red',zorder=12)
            ax.add_patch(cross01)
            ax.add_patch(cross02)

        # counts
        j     = j + 1
        gridi = gridi+1

        if j == ncol:
            j = 0
            i = i + 1

        if i == nrow:
            #gs.tight_layout(fig,h_pad=0,w_pad=0)
            #plt.subplots_adjust(wspace=0., hspace=0.)

            # With cbar_mode="single", cax attribute of all axes are identical.
            cbar = ax.cax.colorbar(imcolor)
            ax.cax.toggle_label(True)
            cbar.ax.yaxis.set_tick_params(color=tickcolor) # tick color
            cbar.ax.spines["bottom"].set_color(axiscolor)  # axes color
            cbar.ax.spines["top"].set_color(axiscolor)
            cbar.ax.spines["left"].set_color(axiscolor)
            cbar.ax.spines["right"].set_color(axiscolor)
            if cbarlabel:
                cbar.ax.set_ylabel(cbarlabel, color=labelcolor) # label

            fig.savefig(outfile, transparent = True)
            fig.clf()
            nmap      = nmap+1
            outfile   = outname + '_nmap{0:02d}'.format(nmap) + '.' + outformat
            fig       = plt.figure(figsize=figsize)
            grid      = ImageGrid(fig, rect=111, nrows_ncols=(nrow,ncol),
                axes_pad=0,share_all=True,cbar_mode=cbar_mode)
            i,j,gridi = [0,0,0]


    # With cbar_mode="single", cax attribute of all axes are identical.
    cbar = ax.cax.colorbar(imcolor)
    ax.cax.toggle_label(True)
    cbar.ax.yaxis.set_tick_params(color=tickcolor) # tick color
    cbar.ax.spines["bottom"].set_color(axiscolor)  # axes color
    cbar.ax.spines["top"].set_color(axiscolor)
    cbar.ax.spines["left"].set_color(axiscolor)
    cbar.ax.spines["right"].set_color(axiscolor)
    if cbarlabel:
        cbar.ax.set_ylabel(cbarlabel,color=labelcolor) # label

    if gridi != gridimax+1 and gridi != 0:
        while gridi != gridimax+1:
            #print gridi
            ax = grid[gridi]
            ax.spines["right"].set_color("none")  # right
            ax.spines["left"].set_color("none")   # left
            ax.spines["top"].set_color("none")    # top
            ax.spines["bottom"].set_color("none") # bottom
            ax.axis('off')
            gridi = gridi+1

        #gs.tight_layout(fig,h_pad=0,w_pad=0)
        #plt.subplots_adjust(wspace=0, hspace=0)
        fig.savefig(outfile, transparent = True)

    # convert N maps into one map
    # ImageMagick is needed
    try:
        outfile = outname + '_nmap*' + '.' + outformat
        sumfile = outname + '.pdf'
        cmdop   = ' -density 600x600 '
        cmd = 'convert' + cmdop + outfile + ' ' + sumfile
        subprocess.call(cmd,shell=True)
    except:
        print 'Cannot convert N maps into one pdf file.'
        print 'Install ImageMagick if you want.'

    return ax



### PV diagram
def pvdiagram(fitsdata,outname,outformat='eps',color=True,cmap='Greys',
    vmin=None,vmax=None,vsys=None,contour=True,clevels=None,ccolor='k',
    vrel=False,logscale=False,x_offset=False,ratio=1.2, prop_vkep=None,fontsize=11,
    lw=1):

    # setting for figures
    plt.rcParams['font.size'] = fontsize           # fontsize

    # output file
    if (outformat == formatlist).any():
        outname = outname + '.' + outformat
    else:
        print 'ERROR\tsingleim_to_fig: Outformat is wrong.'
        return


    # reading fits data
    data, header   = fits.getdata(fitsdata,header=True)


    # figures
    fig = plt.figure(figsize=(11.69,8.27)) # figsize=(11.69,8.27)
    ax  = fig.add_subplot(111)


    # header info.
    naxis      = int(header['NAXIS'])
    noff       = int(header['NAXIS1'])
    nvel       = int(header['NAXIS2'])
    bmaj       = header['BMAJ']*60.*60. # in arcsec
    bmin       = header['BMIN']*60.*60. # in arcsec
    bpa        = header['BPA']          # [deg]
    offlabel   = header['CTYPE1']
    vellabel   = header['CTYPE2']
    thirdlabel = header['BUNIT']
    offunit    = header['CUNIT1']
    restfreq   = header['RESTFRQ'] # Hz
    refval_off = header['CRVAL1']  # in arcsec
    refval_vel = header['CRVAL2']
    refval_vel = clight*(restfreq - refval_vel)/restfreq # Hz --> radio velocity [cm/s]
    refval_vel = refval_vel*1.e-5         # cm/s --> km/s
    refpix_off = header['CRPIX1']
    refpix_vel = header['CRPIX2']
    del_off    = header['CDELT1']  # in arcsec
    del_vel    = header['CDELT2']
    del_vel    = - clight*del_vel/restfreq # Hz --> cm/s
    del_vel    = del_vel*1.e-5             # cm/s --> km/s
    #print refval_vel, del_vel
    #print refval_off, del_off


    # check unit
    if offunit == 'degree' or offunit == 'deg':
        refval_off = refval_off*60.*60.
        del_off    = del_off*60.*60.


    # relative velocity or LSRK
    offlabel = 'offset (arcsec)'
    if vrel:
        refval_vel = refval_vel - vsys
        vellabel   = 'relative velocity (km/s)'
        vcenter    = 0
    else:
        vellabel = 'LSRK velocity (km/s)'
        vcenter  = vsys


    # set extent of an image
    offmin = refval_off + (1 - refpix_off)*del_off - del_off*0.5
    offmax = refval_off + (noff - refpix_off)*del_off + del_off*0.5
    velmin = refval_vel + (1 - refpix_vel)*del_vel - del_vel*0.5
    velmax = refval_vel + (nvel - refpix_vel)*del_vel +del_vel*0.5



    # set axes
    if x_offset:
        data   = data[0,:,:]
        extent = (offmin,offmax,velmin,velmax)
        xlabel = offlabel
        ylabel = vellabel
        hline_params = [vsys,offmin,offmax]
        vline_params = [0.,velmin,velmax]
    else:
        data   = np.rot90(data[0,:,:])
        extent = (velmin,velmax,offmin,offmax)
        xlabel = vellabel
        ylabel = offlabel
        hline_params = [0.,velmin,velmax]
        vline_params = [vcenter,offmin,offmax]


    # set colorscale
    if vmax:
        pass
    else:
        vmax = np.nanmax(data)


    # logscale
    if logscale:
        norm = mpl.colors.LogNorm(vmin=vmin,vmax=vmax)
    else:
        norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)


    # plot images
    if color:
        imcolor = ax.imshow(data, cmap=cmap, origin='lower', extent=extent,norm=norm)

    if contour:
        imcont  = ax.contour(data, colors=ccolor, origin='lower',extent=extent, levels=clevels, linewidths=lw)


    # axis labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


    # lines showing offset 0 and relative velocity 0
    xline = plt.hlines(hline_params[0], hline_params[1], hline_params[2], ccolor, linestyles='dashed', linewidths = 0.5)
    yline = plt.vlines(vline_params[0], vline_params[1], vline_params[2], ccolor, linestyles='dashed', linewidths = 0.5)
    ax.tick_params(which='both', direction='in',bottom=True, top=True, left=True, right=True)


    # aspect ratio
    change_aspect_ratio(ax, ratio)


    # save figure
    fig.savefig(outname, transparent=True)

    return ax