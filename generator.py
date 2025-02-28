import argparse
import astropy
from astropy.io import fits
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('filename', help='name of the generated fits file')
parser.add_argument('size', type=int, help='size of the generated image')
parser.add_argument('--template', default='center', choices=['center', 'grid', 'array'])
parser.add_argument('--arrayfile', nargs=1, help='name of the optional input array used')
args = parser.parse_args()

def centerPoint(beam):
    center = (int)(np.size(beam, 1) / 2)
    beam[center, center] = 255
    return beam

def gridPoints(beam):
    division = (int)(np.size(beam, 1) / 10)
    for x in range(1,10):
        for y in range(1,10):
            beam[x * division, y * division] = 255
    return beam

test_beam = np.zeros((args.size, args.size))
if args.template == 'center':
    test_beam = centerPoint(test_beam)
elif args.template == 'grid':
    test_beam = gridPoints(test_beam)
elif args.template == 'array':
    test_beam = np.load(args.arrayfile[0])


hdu = fits.PrimaryHDU()
fits.writeto(args.filename, test_beam, header=hdu.header, output_verify='exception', overwrite=True, checksum=False)

