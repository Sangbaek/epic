#!/usr/bin/env python3

import uproot
import pyhepmc
import argparse
import awkward as ak
import pandas as pd
import numpy as np

def rotate_2d(x, y, theta):
  x_rot = x*np.cos(theta) - y*np.sin(theta)
  y_rot = x*np.sin(theta) + y*np.cos(theta)
  return x_rot, y_rot

def rotate_3d(x, y, z, theta):
  z, x = rotate_2d(z, x, theta)
  return x, y, z

def main(args):
  direction = args.direction 
  ref_theta = np.arctan(direction[0]/direction[2])

  input_filename = args.input
  output_filename= args.output

  file = uproot.open(input_filename)
  tree = file["TOPEG"]
  entries = tree.num_entries
  entry_start = min(entries-1, max(args.entry_start, 0))
  if args.entry_stop == -1:
    entry_stop = entries
  else:
    entry_stop = min(entries, max(args.entry_stop, entry_start+1))
  print(entry_start, entry_stop)


  data_incoming = {}
  keys_incoming_parts = {"ElBeam": "ElBeam", "EhBeam": "EhBeam"}
  for key in keys_incoming_parts:
  	data_incoming[key] = tree[key].array(library = "ak", entry_start = 0, entry_stop = 1)[0]

  electron_mass = 0.511e-3
  electron_incoming_e = data_incoming["ElBeam"]
  electron_incoming_px = 0
  electron_incoming_py = 0
  electron_incoming_pz = -np.sqrt(electron_incoming_e **2 - electron_mass**2)
  electron_incoming_px, electron_incoming_py, electron_incoming_pz = rotate_3d(electron_incoming_px, electron_incoming_py, electron_incoming_pz, ref_theta)
  electron_incoming = pyhepmc.GenParticle((0, 0, electron_incoming_pz, electron_incoming_e), pid = 11, status = 4)
  electron_incoming.generated_mass = electron_mass

  He4_mass = 3.727504
  He4_incoming_e = data_incoming["EhBeam"]
  He4_incoming_px = 0
  He4_incoming_py = 0
  He4_incoming_pz = np.sqrt(He4_incoming_e **2 - He4_mass**2)
  He4_incoming_px, He4_incoming_py, He4_incoming_pz = rotate_3d(He4_incoming_px, He4_incoming_py, He4_incoming_pz, ref_theta)
  He4_incoming = pyhepmc.GenParticle((0, 0, He4_incoming_pz, He4_incoming_e), pid = 1000020040, status = 4)
  He4_incoming.generated_mass = He4_mass

  with pyhepmc.open(output_filename, "w", precision=7) as output_file:
    for event_number in range(entry_start, entry_stop):
      event = pyhepmc.GenEvent()
      event.event_number = event_number
  
      event.set_units(pyhepmc.Units.MomentumUnit(1), pyhepmc.Units.LengthUnit(0))
  
      event.add_particle(electron_incoming)
      event.add_particle(He4_incoming)
    
      data_outgoing = {}
      keys_outgoing_parts = {"part_px": "px", "part_py": "py", "part_pz": "pz", "part_e": "e"}
      for key in keys_outgoing_parts:
        data_outgoing[keys_outgoing_parts[key]] = tree[key].array(library = "ak", entry_start = event_number, entry_stop = event_number+1)[0]
  
      #electron
      electron_outgoing_px  = data_outgoing["px"] [0]
      electron_outgoing_py  = data_outgoing["py"] [0]
      electron_outgoing_pz  = data_outgoing["pz"] [0]
      electron_outgoing_px, electron_outgoing_py, electron_outgoing_pz = rotate_3d(electron_outgoing_px, electron_outgoing_py, electron_outgoing_pz, ref_theta)
      electron_outgoing_e   = data_outgoing["e"]  [0]
      electron_outgoing     = pyhepmc.GenParticle((electron_outgoing_px, electron_outgoing_py, electron_outgoing_pz, electron_outgoing_e), pid = 11, status = 1)
      electron_outgoing.generated_mass = electron_mass
      event.add_particle(electron_outgoing)
  
      ##virtual photon
      #virtual_photon_internal_px  = electron_incoming_px - electron_outgoing_px
      #virtual_photon_internal_py  = electron_incoming_py - electron_outgoing_py
      #virtual_photon_internal_pz  = electron_incoming_pz - electron_outgoing_pz
      #virtual_photon_internal_px, virtual_photon_internal_py, virtual_photon_internal_pz = rotate_3d(virtual_photon_internal_px, virtual_photon_internal_py, virtual_photon_internal_pz, ref_theta)
      #virtual_photon_internal_e   = electron_incoming_e - electron_outgoing_e
      #virtual_photon_internal     = pyhepmc.GenParticle((virtual_photon_internal_px, virtual_photon_internal_py, virtual_photon_internal_pz, virtual_photon_internal_e), pid = 22, status = 3)
      #virtual_photon_internal.generated_mass = np.sqrt(-virtual_photon_internal_e**2 + virtual_photon_internal_px**2 + virtual_photon_internal_py**2 + virtual_photon_internal_pz**2)
      ## event.add_particle(virtual_photon_internal)
  
      #photon
      photon_outgoing_px  = data_outgoing["px"] [1]
      photon_outgoing_py  = data_outgoing["py"] [1]
      photon_outgoing_pz  = data_outgoing["pz"] [1]
      photon_outgoing_px, photon_outgoing_py, photon_outgoing_pz = rotate_3d(photon_outgoing_px, photon_outgoing_py, photon_outgoing_pz, ref_theta)
      photon_outgoing_e   = data_outgoing["e"]  [1]
      photon_outgoing     = pyhepmc.GenParticle((photon_outgoing_px, photon_outgoing_py, photon_outgoing_pz, photon_outgoing_e), pid = 22, status = 1)
      photon_outgoing.generated_mass = 0
      event.add_particle(photon_outgoing)
  
      #He4
      He4_outgoing_px  = data_outgoing["px"] [2]
      He4_outgoing_py  = data_outgoing["py"] [2]
      He4_outgoing_pz  = data_outgoing["pz"] [2]
      He4_outgoing_px, He4_outgoing_py, He4_outgoing_pz = rotate_3d(He4_outgoing_px, He4_outgoing_py, He4_outgoing_pz, ref_theta)
      He4_outgoing_e   = data_outgoing["e"]  [2]
      He4_outgoing     = pyhepmc.GenParticle((He4_outgoing_px, He4_outgoing_py, He4_outgoing_pz, He4_outgoing_e), pid = 1000020040, status = 1)
      He4_outgoing.generated_mass = He4_mass
      event.add_particle(He4_outgoing)
  
      v1 = pyhepmc.GenVertex()
      v1.add_particle_in(electron_incoming)
      # v1.add_particle_out(virtual_photon_internal)
      v1.add_particle_out(electron_outgoing)
      # v2.add_particle_in(virtual_photon_internal)
      v1.add_particle_in(He4_incoming)
      v1.add_particle_out(He4_outgoing)
      v1.add_particle_out(photon_outgoing)
      event.add_vertex(v1)

      output_file.write(event)

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input')
  parser.add_argument('-o', '--output')
  parser.add_argument('-t', '--target', default = 'He4')
  parser.add_argument('-s', '--entry_start', default = 0, type = int)
  parser.add_argument('-e', '--entry_stop', default = -1, type = int)
  parser.add_argument('-dir', '--direction', nargs = 3, default = [-0.03203692472855905, 0, 0.9994866859813274], type = float)

  args = parser.parse_args()
  main(args)
