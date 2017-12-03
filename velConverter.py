

import numpy as np
from numpy import cos, sin

import astropy.coordinates as coord
import astropy.units as u
from math import *


# This is the default circular velocity and LSR peculiar velocity of the Sun
# TODO: make this a config item?
VCIRC = 220. # u.km/u.s
VLSR = [10., 5.25, 7.17] # *u.km/u.s

def vgsr_to_vhel(gl, gb, vgsr, vcirc=VCIRC, vlsr=VLSR):
    """ Convert a radial velocity in the Galactic standard of rest (GSR) to
        a barycentric radial velocity.

        Parameters
        ----------
        coordinate : l and b are galactic coordinates (gl, gb)
        vgsr : :class:`~astropy.units.Quantity`
            GSR line-of-sight velocity.
        vcirc : :class:`~astropy.units.Quantity`
            Circular velocity of the Sun.
        vlsr : :class:`~astropy.units.Quantity`
            Velocity of the Sun relative to the local standard
            of rest (LSR).

        Returns
        -------
        vhel : :class:`~astropy.units.Quantity`
            Radial velocity in a barycentric rest frame.

    """
    l = gl*pi/180.
    b = gb*pi/180.


    #if not isinstance(vgsr, u.Quantity):
        #raise TypeError("vgsr must be a Quantity subclass")

    # compute the velocity relative to the LSR
    lsr = vgsr - vcirc*sin(l)*cos(b)

    # velocity correction for Sun relative to LSR
    v_correct = vlsr[0]*cos(b)*cos(l) + \
        vlsr[1]*cos(b)*sin(l) + \
        vlsr[2]*sin(b)
    vhel = lsr - v_correct

    return vhel

def vhel_to_vgsr(gl, gb, vhel, vcirc=VCIRC, vlsr=VLSR):
    """ Convert a velocity from a heliocentric radial velocity to
        the Galactic standard of rest (GSR).

        Parameters
        ----------
        coordinate : :class:`~astropy.coordinates.SkyCoord`
            An Astropy SkyCoord object or anything object that can be passed
            to the SkyCoord initializer.
        vhel : :class:`~astropy.units.Quantity`
            Barycentric line-of-sight velocity.
        vcirc : :class:`~astropy.units.Quantity`
            Circular velocity of the Sun.
        vlsr : :class:`~astropy.units.Quantity`
            Velocity of the Sun relative to the local standard
            of rest (LSR).

        Returns
        -------
        vgsr : :class:`~astropy.units.Quantity`
            Radial velocity in a galactocentric rest frame.

    """
    l = gl*pi/180.
    b = gb*pi/180.

    if not isinstance(vhel, u.Quantity):
        raise TypeError("vhel must be a Quantity subclass")

    lsr = vhel + vcirc*sin(l)*cos(b)

    # velocity correction for Sun relative to LSR
    v_correct = vlsr[0]*cos(b)*cos(l) + \
        vlsr[1]*cos(b)*sin(l) + \
        vlsr[2]*sin(b)
    vgsr = lsr + v_correct

    return vgsr


