#!/usr/bin/env python3

import uproot
import awkward as ak
import matplotlib.pyplot as plt
import numpy as np

file = uproot.open("output/Out-EIClow-M3-Nh.root")
tree = file["events"]

df_gen = ak.to_dataframe(tree["MCParticles"].array(library = "ak"))
columns = {string: string.replace("MCParticles.", "") for string in df_gen.columns}
df_gen = df_gen.rename(columns = columns)
columns = {string: string.replace("vertex.", "v") for string in df_gen.columns}
df_gen = df_gen.rename(columns = columns)
columns = {string: string.replace("momentum.", "p") for string in df_gen.columns}
df_gen = df_gen.rename(columns = columns)

df_primary = df_gen.loc[df_gen.index.get_level_values('subentry') == 4, :]

df = ak.to_dataframe(tree["ForwardRomanPotHits"].array(library = "ak"))
columns = {string: string.replace("ForwardRomanPotHits.", "") for string in df.columns}
df = df.rename(columns = columns)
columns = {string: string.replace("position.", "") for string in df.columns}
df = df.rename(columns = columns)
columns = {string: string.replace("momentum.", "p") for string in df.columns}
df = df.rename(columns = columns)

dir = [-0.03203692472855905, 0, 0.9994866859813274]
# dir = [0, 0, 1]
def inner(dir1, dir2):
  return dir1[0]*dir2[0] + dir1[1]*dir2[1] + dir1[2]*dir2[2]
def mag2(dir1):
  return inner(dir1, dir1)
def mag(dir1):
  return np.sqrt(mag2(dir1))
def angle(dir1, dir2):
  return np.arccos(inner(dir1, dir2)/mag(dir1)/mag(dir2))

generated_angles = angle(dir, [df_primary.px, df_primary.py, df_primary.pz])
print(np.histogram(generated_angles, bins = 100))

df = ak.to_dataframe(tree["ForwardRomanPotHits"].array(library = "ak"))
columns = {string: string.replace("ForwardRomanPotHits.", "") for string in df.columns}
df = df.rename(columns = columns)
columns = {string: string.replace("position.", "") for string in df.columns}
df = df.rename(columns = columns)
columns = {string: string.replace("momentum.", "p") for string in df.columns}
df = df.rename(columns = columns)

hit_angles = angle(dir, df.loc[:, ["x", "y", "z"]].to_numpy().T)
print(np.histogram(hit_angles, bins = 100))
