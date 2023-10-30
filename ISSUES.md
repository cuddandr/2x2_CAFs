List of issues/bugs/strange things noticed while working with the CAFs.

Updated for PicoRun4.1 CAFs (2023/10/30)

### Common issues

Things noticed in both structured and flat CAFs:
+ Repeated G4IDs for primary particles
    - Affected variables: `rec.mc.nu.prim.G4ID`
+ Truth particle type for in common branch always 3 (kSecondary)
    - Affected variables: `rec.common.ixn.dlp.part.dlp.truth.type
+ Variables assumed for truth matching only contain values of -1 or 0 (default values according to the SR code)
    - Affected variables: `rec.nd.lar.dlp.tracks.truth.{ixn,part,type}`
    - Affected variables: `rec.nd.lar.dlp.showers.truth.{ixn,part,type}`

### Flat CAF issues

+ Kinetic energy for non-calometric calculation (E_method) always set to -1e-3
    - Affected variables: `rec.common.ixn.dlp.part.dlp.E`

### Structured CAF issues

+ Reading anything from the `rec.mc.nu` branches causes a crash after some number of entries
    - e.g. Reading `PicoRun4.1_1E17_RHC.larnd.00000.caf.root` will crash at entry 38 using `cafTree->Scan("mc.nu.id")`

### Not necessarily issues

+ Empty spills? What causes a spill to be empty?
