# Dump Routine Module Alias Audit

Reference for Phase 3 (thread-safe conversion). Lists all module-level
aliases from POINT_TO_MESH used by each output function.

## GAS_PHASE_OUTPUT (dump.f90:7384, ~1432 lines, 69 aliases)

**Velocity**: U, V, W, US, VS, WS
**Thermodynamic**: RHO, TMP, H, HS, D, DS, KRES, MU, MU_DNS
**Species**: ZZ
**Radiation**: UII, QR, CHI_R, Q
**Turbulence**: CSD2, LES_FILTER_WIDTH, D_Z_MAX, MIX_TIME
**Pressure**: PBAR, PRESSURE_ZONE
**Work arrays**: CARTVELDIV
**Cell/geometry**: CELL, CELL_INDEX, IBAR, JBAR, KBAR, XC, YC, ZC, DX, DY, DZ, RDXN, RDYN, RDZN, RDX, RDY, RDZ
**Flux**: FVX, FVY, FVZ, FVX_B, FVY_B, FVZ_B, ADV_FX, ADV_FY, ADV_FZ, DIF_FX, DIF_FY, DIF_FZ
**Source**: REAC_SOURCE_TERM, CHEM_SUBIT
**Particle**: AVG_DROP_DEN, AVG_DROP_RAD, AVG_DROP_TMP, AVG_DROP_AREA
**CC_IBM**: CCVAR, FCVAR, RC_FACE, CUT_FACE
**Diagnostic**: PP_RESIDUAL, POIS_ERR, RESMAX, VN, CFL

## SOLID_PHASE_OUTPUT (dump.f90:8832, ~826 lines, 21 aliases)

**Velocity**: U, V, W
**Thermodynamic**: RHO, TMP, H, KRES, D
**Cell/geometry**: CELL, CELL_INDEX, DX, DY, DZ
**Flux**: ADV_FX, ADV_FY, ADV_FZ, DIF_FX, DIF_FY, DIF_FZ
**CC_IBM**: CUT_FACE, CUT_CELL

## DUMP_SLCF (dump.f90:5968, ~543 lines)

Uses: CELL, CELL_INDEX, CELL_COUNT, WORK1 (->B), WORK2 (->S),
WORK3 (->QUANTITY), QQ, IBP1, JBP1, KBP1, IBAR, JBAR, KBAR,
XC, YC, ZC, X, Y, Z.
Calls GAS_PHASE_OUTPUT (69 aliases above).

## DUMP_BNDF (dump.f90:10640, ~220 lines)

Uses: WALL, BOUNDARY_COORD, BOUNDARY_PROP1, CELL, PP, PPN, IBK, PATCH.
Calls SOLID_PHASE_OUTPUT (21 aliases above).

## DUMP_ISOF (dump.f90:4314, ~181 lines)

Uses: CELL, CELL_INDEX, WORK1 (->B), WORK2 (->S), WORK3, WORK4,
QQ, QQ2, IBP1, JBP1, KBP1, IBAR, JBAR, KBAR, IBLK.
Calls GAS_PHASE_OUTPUT.

## DUMP_PART (dump.f90:4192, ~114 lines)

Uses: LAGRANGIAN_PARTICLE, NLP, OMESH.
Does NOT call GAS_PHASE_OUTPUT or SOLID_PHASE_OUTPUT directly.

## DUMP_SMOKE3D (dump.f90:4503, ~65 lines)

Uses: WORK3 (->FF), QQ, IBP1, JBP1, KBP1, IBAR, JBAR, KBAR.
Calls GAS_PHASE_OUTPUT.

## Key Observations

1. GAS_PHASE_OUTPUT already takes NM as parameter and uses MESHES(NM)
   in CC_IBM paths — partial M% usage already exists.
2. SOLID_PHASE_OUTPUT already takes NM and uses optional index parameters
   to access WALL/BOUNDARY/CFACE data — relatively self-contained.
3. DUMP_PART does NOT use POINT_TO_MESH aliases for physics — it accesses
   LAGRANGIAN_PARTICLE arrays directly. Easiest to convert.
4. DUMP_SLCF and DUMP_ISOF are the hardest due to WORK array aliasing
   and GAS_PHASE_OUTPUT's 69 dependencies.
