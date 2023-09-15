import numpy as np
import uproot as ur
import matplotlib.pyplot as plt

file_name = "./files/outputCAF_notruth_27023797_10.flat.root"

print("Opening {}".format(file_name))
caf_file = ur.open(file_name)
caf_tree = caf_file['cafTree']

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
part_px   = np.concatenate(caf_tree['rec.common.ixn.dlp.part.dlp.p.x'].array(library="np"))
part_py   = np.concatenate(caf_tree['rec.common.ixn.dlp.part.dlp.p.y'].array(library="np"))
part_pz   = np.concatenate(caf_tree['rec.common.ixn.dlp.part.dlp.p.z'].array(library="np"))

proton_mask = (part_pdg == 2212) #& (part_cont == True)
proton_kinE = part_E[proton_mask]
proton_pvec = np.stack((part_px[proton_mask], part_py[proton_mask], part_pz[proton_mask]), axis=1)

print("Making histogram(s)...")

nbins=50
style_hist = {"alpha" : 0.5}

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

fig = plt.figure()
plt.hist(proton_p, bins=nbins, range=[0,1000], label='All protons', **style_hist)
plt.hist(proton_p[contained], bins=nbins, range=[0,1000], label='Only contained', **style_hist)
plt.ylabel("Num. tracks")
plt.xlabel("Momentum [MeV]")
plt.legend()
plt.savefig("proton_momentum.png")

#pmass = 0.938
#proton_q = np.sqrt(np.square(proton_kinE + pmass) - pmass**2) * 1000

fig = plt.figure()
plt.hist(proton_kinE, bins=nbins, range=[0,0.5], label='All protons', **style_hist)
plt.hist(proton_kinE[contained], bins=nbins, range=[0,0.5], label='Only contained', **style_hist)
plt.ylabel("Num. tracks")
plt.xlabel("Kinetic Energy [GeV]")
plt.legend()
plt.savefig("proton_kinE.png")

print("Finished.")
