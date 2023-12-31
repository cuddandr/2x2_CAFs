import numpy as np
import uproot as ur
import matplotlib.pyplot as plt
import argparse

default_file = "./files/flat/outputCAF_noGENIEtruth_27023797_10.flat.root"

parser = argparse.ArgumentParser(description='Make plots from CAF(s)')
parser.add_argument('file', metavar='file', nargs='?', default=default_file, help='input CAF file')
args = parser.parse_args()

file_name = args.file
print("Opening {}".format(file_name))
caf_file = ur.open(file_name)
caf_tree = caf_file['cafTree']

num_ixn = caf_tree['rec.common.ixn.ndlp'].array(library="np")

reco_vtx_x = np.concatenate(caf_tree['rec.common.ixn.dlp.vtx.x'].array(library="np"))
reco_vtx_y = np.concatenate(caf_tree['rec.common.ixn.dlp.vtx.y'].array(library="np"))
reco_vtx_z = np.concatenate(caf_tree['rec.common.ixn.dlp.vtx.z'].array(library="np"))
true_vtx_x = np.concatenate(caf_tree['rec.mc.nu.vtx.x'].array(library="np"))
true_vtx_y = np.concatenate(caf_tree['rec.mc.nu.vtx.y'].array(library="np"))
true_vtx_z = np.concatenate(caf_tree['rec.mc.nu.vtx.z'].array(library="np"))

reco_vtx_x = reco_vtx_x[np.isfinite(reco_vtx_x)]
reco_vtx_y = reco_vtx_y[np.isfinite(reco_vtx_y)]
reco_vtx_z = reco_vtx_z[np.isfinite(reco_vtx_z)]
true_vtx_x = true_vtx_x[np.isfinite(true_vtx_x)]
true_vtx_y = true_vtx_x[np.isfinite(true_vtx_y)]
true_vtx_z = true_vtx_x[np.isfinite(true_vtx_z)]

num_reco_vtx = len(reco_vtx_x)
num_true_vtx = len(true_vtx_x)

part_cont = np.concatenate(caf_tree['rec.common.ixn.dlp.part.dlp.contained'].array(library="np"))
part_pdg  = np.concatenate(caf_tree['rec.common.ixn.dlp.part.dlp.pdg'].array(library="np"))
part_E    = np.concatenate(caf_tree['rec.common.ixn.dlp.part.dlp.E'].array(library="np"))
part_Em   = np.concatenate(caf_tree['rec.common.ixn.dlp.part.dlp.E_method'].array(library="np"))
part_px   = np.concatenate(caf_tree['rec.common.ixn.dlp.part.dlp.p.x'].array(library="np"))
part_py   = np.concatenate(caf_tree['rec.common.ixn.dlp.part.dlp.p.y'].array(library="np"))
part_pz   = np.concatenate(caf_tree['rec.common.ixn.dlp.part.dlp.p.z'].array(library="np"))

proton_mask = (part_pdg == 2212) #& (part_cont == True)
proton_kinE = part_E[proton_mask]
proton_pvec = np.stack((part_px[proton_mask], part_py[proton_mask], part_pz[proton_mask]), axis=1)

print("Making histogram(s)...")

fig_ixn = plt.figure()
plt.hist(num_ixn, range=[0,20], bins=20)
plt.ylabel("Num. occurances")
plt.xlabel("Num. interactions/spill")
plt.xticks(np.arange(0, 21))
plt.savefig("num_ixn.png")

nbins=50
style_hist = {"alpha" : 0.5}

E_calo  = part_E[part_Em == 3]
E_other = part_E[part_Em != 3]

fig_energy = plt.figure()
plt.hist(E_calo, range=[-0.1,0.9], bins=nbins, label='Calometric energy', **style_hist)
plt.hist(E_other, range=[-0.1,0.9], bins=nbins, label='MCS/range energy', **style_hist)
plt.ylabel("Num. events")
plt.xlabel("Kinetic energy [GeV]")
plt.legend(frameon=False, fontsize=14)
plt.savefig("kin_energy.png")

reco_label = "reco ({} entries)".format(num_reco_vtx)
true_label = "true ({} entries)".format(num_true_vtx)

fig_vtx, axes = plt.subplots(1, 3, figsize=(15, 5))
axes[0].hist(reco_vtx_x, bins=nbins, label=reco_label, **style_hist)
axes[0].hist(true_vtx_x, bins=nbins, label=true_label, **style_hist)
axes[0].set_ylabel("Num. events")
axes[0].set_xlabel("Vtx X [cm]")
axes[0].legend(frameon=False, fontsize=14)

axes[1].hist(reco_vtx_y, bins=nbins, label=reco_label, **style_hist)
axes[1].hist(true_vtx_y, bins=nbins, label=true_label, **style_hist)
axes[1].set_ylabel("Num. events")
axes[1].set_xlabel("Vtx Y [cm]")
axes[1].legend(frameon=False, fontsize=14)

axes[2].hist(reco_vtx_z, bins=nbins, label=reco_label, **style_hist)
axes[2].hist(true_vtx_z, bins=nbins, label=true_label, **style_hist)
axes[2].set_ylabel("Num. events")
axes[2].set_xlabel("Vtx Z [cm]")
axes[2].legend(frameon=False, fontsize=14)

fig_vtx.savefig("vtx_reco_true.png")

proton_p = np.linalg.norm(proton_pvec, axis=1)
contained = (part_cont[part_pdg == 2212] == True)

nbins=50
style_hist = {"alpha" : 0.5}

fig_proton, axes = plt.subplots(1, 2, figsize=(15, 5))
axes[0].hist(proton_p, bins=nbins, range=[0,1.0], label='All protons', **style_hist)
axes[0].hist(proton_p[contained], bins=nbins, range=[0,1.0], label='Only contained', **style_hist)
axes[0].set_ylabel("Num. tracks")
axes[0].set_xlabel("Momentum [GeV]")
axes[0].legend()

axes[1].hist(proton_kinE, bins=nbins, range=[0,0.5], label='All protons', **style_hist)
axes[1].hist(proton_kinE[contained], bins=nbins, range=[0,0.5], label='Only contained', **style_hist)
axes[1].set_ylabel("Num. tracks")
axes[1].set_xlabel("Kinetic Energy [GeV]")
axes[1].legend()

fig_proton.savefig("proton_kinematics.png")

print("Finished.")
