List of issues/bugs/strange things noticed while working with the CAFs.

### Common issues

Things noticed in both structured and flat CAFs:
+ Truth vertex position centered at zero --> missing TPC offset
    - Affected variables: `rec.mc.nu.vtx.{x,y,z}`
+ Number of reconstructed interactions in `rec.nd.lar.ndlp` always 0
+ Only charged leptons listed as primary particles?
    - Affected variables: `rec.mc.nu.prim.pdg`
+ Interactions with no primary particles?
+ Too few primary particles from interaction? Most common is only a single primary particle, which seems far too low given the beam energy.

Most of this may or not may be expected given the missing GENIE record
+ Indices of SRTrueInteractions(s) always 0
    - Affected variables: `rec.common.ixn.dlp.truth`
+ Indices of SRInteraction in the SRTruthBranch always -1
    - Affected variables: `rec.common.ixn.dlp.part.dlp.truth.ixn`
+ Variables assumed for truth matching only contain values of -1 or 0 (default values according to the SR code)
    - Affected variables: `rec.nd.lar.dlp.tracks.truth.{ixn,part,type}`
    - Affected variables: `rec.nd.lar.dlp.showers.truth.{ixn,part,type}`

### Flat CAF issues

+ Kinetic energy for non-calometric calculation always set to -1e-3
    - Affected variables: `rec.common.ixn.dlp.part.dlp.E`
+ Number of interactions per spill is always an even number
    - Affected variables: `rec.common.ixn.ndlp`

### Structured CAF issues
