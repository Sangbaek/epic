#!/usr/bin/env python3

import uproot
import pyhepmc
import argparse
import awkward as ak
import pandas as pd
import numpy as np

He4_mass = 3.727504

def rotate_2d(x, y, theta):
  x_rot = x*np.cos(theta) - y*np.sin(theta)
  y_rot = x*np.sin(theta) + y*np.cos(theta)
  return x_rot, y_rot

def rotate_3d(x, y, z, theta):
  z, x = rotate_2d(z, x, theta)
  return x, y, z

def main(args):
  output_filename= args.output
  energies = np.random.uniform(args.emin, args.emax, args.numberOfEvents)
  momenta =   np.sqrt(energies **2 - He4_mass**2)

  costhetas  = np.random.uniform(np.cos(args.thetamax), np.cos(args.thetamin), args.numberOfEvents)
  thetas   = np.arccos(costhetas)
  print(thetas.max(), thetas.min())
  phis     = np.random.uniform(0, 2*np.pi, args.numberOfEvents)

  if args.direction:
    direction = args.direction
  else:
    direction = [-0.03203692472855905, 0, 0.9994866859813274]
  ref_theta = np.arctan(direction[0]/direction[2])
  # print(direction, ref_theta, thetas, phis)

  pxs = momenta*np.sin(thetas)*np.cos(phis)
  pys = momenta*np.sin(thetas)*np.sin(phis)
  pzs = momenta*np.cos(thetas)
  pxs, pys, pzs = rotate_3d(pxs, pys, pzs, ref_theta)


  with pyhepmc.open(output_filename, "w", precision=7) as output_file:
    for event_number in range(args.numberOfEvents):

      He4_outgoing_e = energies[event_number]
      He4_outgoing_px = pxs[event_number]
      He4_outgoing_py = pys[event_number]
      He4_outgoing_pz = pzs[event_number]
      He4_outgoing = pyhepmc.GenParticle((He4_outgoing_px
        , He4_outgoing_py
        , He4_outgoing_pz
        , He4_outgoing_e), pid = 1000020040, status = 1)
      He4_incoming = pyhepmc.GenParticle((He4_outgoing_px
        , He4_outgoing_py
        , He4_outgoing_pz
        , He4_outgoing_e), pid = 1000020040, status = 4)
      #He4_outgoing.generated_mass = He4_mass

      event = pyhepmc.GenEvent()
      event.event_number = event_number
  
      event.set_units(pyhepmc.Units.MomentumUnit(1), pyhepmc.Units.LengthUnit(0))
  
      event.add_particle(He4_incoming)
      event.add_particle(He4_outgoing)
      
      v1 = pyhepmc.GenVertex()
      v1.add_particle_in(He4_incoming)
      v1.add_particle_out(He4_outgoing)
      event.add_vertex(v1)

      output_file.write(event)

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-o', '--output')
  parser.add_argument('-t', '--target', default = 'He4')
  parser.add_argument('-emin', '--emin', default = 40, type = float)
  parser.add_argument('-emax', '--emax', default = 41, type = float)
  parser.add_argument('-thetamin', '--thetamin', default = 0, type = float)
  parser.add_argument('-thetamax', '--thetamax', default = 0.05, type = float)
  parser.add_argument('-dir', '--direction', nargs = 3, default = None, type = float)
  parser.add_argument('-n', '--numberOfEvents', default = 0, type = int)

  args = parser.parse_args()

  main(args)
