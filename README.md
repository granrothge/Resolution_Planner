# Resolution Planner

=======================

A set of python scripts for calculating Energy resolution and relative intensity for
Neutron Direct Geometry Spectrometers.

The main script to execute is res_ints_4.py

There are two major functions to use.
1. plot_flux(nu,Ei,Ef,Spec)

   nu is an array of Fermi chopper frequencies

   Ei is the incident Energy of the neutron

   Ef is the final Energy of the neutron

   Spec is an instance of the spectrometer class

2. plot_res_omega(nu,Ei,omega,Spec)

   nu is the Fermi chopper speed

   Ei is the incident energy

   omega is a numpy array of energy transfers

   Spec is an instance of the spectrometer class

The Following instances of the spectrometer class are defined at the end of the main file
   SEQUOIA, SEQUOIA_sloppy, SEQUOIA_700_superfine,
   SEQUOIA_1000, ARCS, ARCS_700_fine, ARCS_700_sloppy,
   ARCS_700_superfine

## UB tools

UB.py is a library of tools for generating and working with primarily B matrices
There are also rotation matrix generators to generate U matrices
