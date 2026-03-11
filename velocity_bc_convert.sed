# sed script to convert remaining mesh-specific variables in VELOCITY_BC_KERNEL
# Run with: sed -i -f velocity_bc_convert.sed Source/velo.f90

# CELL references (within VELOCITY_BC_KERNEL only, lines 703-1512)
/^SUBROUTINE VELOCITY_BC_KERNEL/,/^END SUBROUTINE VELOCITY_BC_KERNEL/{
s/\bCELL(/M%CELL(/g
s/\bWALL(/M%WALL(/g
s/\bBOUNDARY_COORD(/M%BOUNDARY_COORD(/g
s/\bBOUNDARY_PROP1(/M%BOUNDARY_PROP1(/g
s/\bMU(/M%MU(/g
s/\bTMP(/M%TMP(/g
s/\bZC(/M%ZC(/g
s/\bOMESH(/M%OMESH(/g
s/\bIBAR\b/M%IBAR/g
s/\bJBAR\b/M%JBAR/g
s/\bKBAR\b/M%KBAR/g
s/\bIBP1\b/M%IBP1/g
s/\bJBP1\b/M%JBP1/g
s/\bKBP1\b/M%KBP1/g
s/\bDRAG_UVWMAX\b/M%DRAG_UVWMAX/g
s/\bT_USED(/M%T_USED(/g
}
