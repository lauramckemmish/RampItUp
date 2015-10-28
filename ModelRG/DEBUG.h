! This is the file used to customise the operation and output of RampItUp. 
! RampItUp MUST be recompiled for these options to take effect; they are not input to the program for efficiency. 
! Turn on an option by removing the comment (the '!' symbol)
! Turn off an option by adding a comment (the '!' symbol)


! Turn the below all on if you want perform a full calculation. 
! If you want to time just the ramp-containing integral calculation, comment the KEEPGGGG and uncomment TIMINGRUN 
#define KEEPONEELEC
#define KEEPDIGEST
#define KEEPGGGG
#define KEEPLONG
#define KEEPSHORT
#define KEEPGGGG
!#define TIMINGRUN


! Turn the below on if you want to print the one- and two-electron integrals in the output file
#define ONE_ELEC_PRT 
#define INTPRINT


! Extra output options
!#define BASISPRT     ! prints details of basis set used and significant shell-pairs
!#define MDLPRT       ! prints out details of models for Rg shell-pairs
!#define INTCNTPRT    ! counts the number of primitive two-electron integrals approximately
!#define WRITESCFTIMINGS ! writes scf timings for every iteration; turn off for shorter output

! Turn the below on if you want to print the one- and two-electron integrals in the Integrals output file
!#define WRITETWOEE
!#define WRITEONE


