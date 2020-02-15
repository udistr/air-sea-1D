# air-sea-1D

A code for a simple 1D coupled ocean-atmosphere model.

Prognostic variables:
* Atmosphere: Heat, horizontal momentum (1 direction) and moisture.
* Ocean: Heat and momentum (salinity not yet implemented)

Grid:
* height coordinates

Processess included:
* vertical advaction, explicit integration
* vetical diffustion:
  * Large and Yeager (1994) for the ocean
  * Holtslag and Boville (1993) for the atmosphere
* Bulk formulea for the air-sea interface based on Large and Yeager (2004)
* Solar radiation including transmision to the deeper ocean
