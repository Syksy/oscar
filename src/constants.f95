      MODULE constants
        USE, INTRINSIC :: iso_c_binding
        IMPLICIT NONE

        ! ** Double precision (i.e accuracy) **
        INTEGER(KIND=c_int), PARAMETER :: dp = c_double 
        INTEGER(KIND=c_int), PARAMETER :: ip = c_int
        !INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)  

        
  
      END MODULE constants
