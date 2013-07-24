"""
Script to make VLT/HAWK-I Service Mode finder charts.

The finder chart requirements are listed at:
http://www.eso.org/sci/observing/phase2/SMGuidelines/FindingCharts.HAWKI.html

Background images for 2MASS, SDSS, DSS can be obtained from the Montage Image Mosaic Service:
http://hachi.ipac.caltech.edu:8080/montage

The script is optimized for 15' (0.25 deg) Montage cutouts.

"""
import numpy as np
import threedhst
import unicorn
import pyfits

def make_regions(ra=177.39458, dec=22.401338, pa=0):
    """
    Make region files for the HAWK-I FOV
    """
    
    DX = 217./3600.
    sep = 15./3600.
    
    cosdec = np.cos(dec/360.*2*np.pi)
    DY = DX #*cosdec
    
    box_x = np.array([-0.5,0.5,0.5,-0.5,-0.5])*DX
    box_y = np.array([-0.5,-0.5,0.5,0.5,-0.5])*DY
    
    
    full_poly = ''
    
    hx = []
    hy = []
    for ix, iy in zip([1,-1,-1,1], [-1,-1,1,1]):
        poly = np.zeros(10)
        px = ra*0 + box_x + ix*(DX/2+sep/2.)
        py = dec*0 + box_y + iy*(DY/2+sep/2.)
        px, py = threedhst.utils.xyrot(px, py, pa, ccw=False)
        poly[::2] = px/cosdec + ra
        poly[1::2] = py + dec
        hx.append(px/cosdec + ra)
        hy.append(py + dec)
        polystr = '%f' %(poly[0])
        for p in poly[1:]:
            polystr += ',%f' %(p)
        #
        #fp.write('polygon(%s)\n' %(polystr))
        full_poly += polystr+', '
    #
    fp = open('hawki.reg','w')
    fp.write('fk5\n')
    fp.write('polygon(%s)\n' %(full_poly[:-2]))
    fp.close()
    
    x10 = np.array([-0.5,0.5,0.5,-0.5,-0.5])*(DX*2+sep)
    y10 = np.array([-0.5,-0.5,0.5,0.5,-0.5])*(DY*2+sep)
    x10, y10 = threedhst.utils.xyrot(x10, y10, pa, ccw=False)
    p2 = ','.join(['%.6f,%.6f' %(x10[i]/cosdec+ra, y10[i]+dec) for i in range(5)])
    
    fp = open('fov.reg','w')
    fp.write('fk5\n')
    fp.write('polygon(%s) # color=magenta\n' %(p2))
    fp.close()
    
    hx.append(x10/cosdec+ra)
    hy.append(y10)
    
    x10 = np.array([-0.5,0.5,0.5,-0.5,-0.5])*10./60
    y10 = np.array([-0.5,-0.5,0.5,0.5,-0.5])*10./60
    p2 = ','.join(['%.6f,%.6f' %(x10[i]/cosdec+ra, y10[i]+dec) for i in range(5)])
    
    fp = open('10arcmin.reg','w')
    fp.write('fk5\n')
    fp.write('polygon(%s) # color=black\n' %(p2))
    fp.close()
    
    hx.append(x10/cosdec+ra)
    hy.append(y10+dec)
    
    return hx, hy
    
def make_finder_charts():
    import finder
    
    px, py = finder.make_regions(ra=3.5313884, dec=-30.389218, pa=53.8-90)
    
    fig, ax = finder.make_finder(image='A2744_DSS_r.fits', ra=3.5313884, dec=-30.389218, pa=53.8-90, target='A2744 (Ks)', image_label='DSS-r',     progID = '092.A-0472(A)', PI = 'G. Brammer')
    unicorn.plotting.savefig(fig, 'A2744_DSS.pdf')
    
    fig, ax = finder.make_finder(image='A2744_2mass_Ks.fits', ra=3.5313884, dec=-30.389218, pa=53.8-90, target='A2744 (Ks)', image_label='2MASS-Ks',     progID = '092.A-0472(A)', PI = 'G. Brammer', vscale=(833.18, 850.34))
    unicorn.plotting.savefig(fig, 'A2744_2mass.pdf')
    
def make_finder(image='A2744_DSS_r.fits', ra=3.5313884, dec=-30.389218, pa=53.8-90, target='A2744 (Ks)', image_label='DSS-r',     progID = '092.A-0472(A)', PI = 'G. Brammer', vscale=None):
    
    import matplotlib.pyplot as plt
    import astropy.wcs 
    
    import finder
    
    # plt.gray()
    
    im = pyfits.open(image)[0].data
    wcs = astropy.wcs.WCS(pyfits.getheader(image))
    
    px, py = finder.make_regions(ra=ra, dec=dec, pa=pa)
    
    fig = unicorn.plotting.plot_init(xs=6, aspect=1.1, left=0.01, right=0.01, top=0.01, bottom=0.01)
    ax = fig.add_axes((0.01, 0.01, 0.98, 0.88)) #(111)
    if vscale is None:
        vscale = np.percentile(im, [0.1,99.8])
        print 'Vscale: %.2f, %.2f' %(vscale[0], vscale[1])
        
    ai = ax.imshow(0-im, vmin=-vscale[1], vmax=-vscale[0], aspect='auto', interpolation='nearest', label=image_label)
    
    colors = ['magenta', 'black','black','black']
    for i in range(4):
        bx, by = wcs.wcs_world2pix(px[i], py[i], 0)
        ax.plot(bx, by, color=colors[i], label='HAWKI, Chip 1'*(i==0))
    
    ### 10x10 box
    bx, by = wcs.wcs_world2pix(px[-1], py[-1], 0)
    ax.plot(bx, by, color='black', linestyle='--', label=r"$10\times10$ arcmin")
    
    ax.set_xlim(0, im.shape[1])
    ax.set_xticks([0, im.shape[1]])
    ax.set_ylim(0, im.shape[0])
    ax.set_yticks([0, im.shape[0]])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    
    ax.text(0.02, 0.98, 'Image: %s' %(image_label), ha='left', va='top', transform=ax.transAxes, size=14)
        
    axt = fig.add_axes((0.01, 0.89, 0.98, 0.1))
    axt.text(0.02, 0.9, '%s,  PI: %s' %(progID, PI), ha='left', va='top', size=17, transform=axt.transAxes)
    axt.text(0.02, 0.1, 'Target: %s' %(target), ha='left', va='bottom', size=17, transform=axt.transAxes)
    axt.text(0.98, 0.1, 'PA = %.1f' %(pa), ha='right', va='bottom', size=12, transform=axt.transAxes)

    axt.set_xticklabels([])
    axt.set_yticklabels([])
    axt.set_xlim(0,1); axt.set_xticks([])
    axt.set_ylim(0,1); axt.set_yticks([])
    
    ax.legend(loc='lower right')
    
    ax.arrow(0.15,0.05, -0.1, 0, transform=ax.transAxes, color='black', head_length=0.01, head_width=0.01)
    ax.arrow(0.15,0.05, 0, 0.1, transform=ax.transAxes, color='black', head_length=0.01, head_width=0.01)
    ax.text(0.05, 0.05+0.02, 'E', ha='center', va='bottom', transform=ax.transAxes)
    ax.text(0.15-0.02, 0.15, 'N', ha='right', va='center', transform=ax.transAxes)
    
    return fig, ax
        
         
    
    