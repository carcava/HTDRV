!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE ht_pw_drv(lib_comm,nim,npt,npl,nta,nbn,ndg,retval,infile) BIND(C)
  !----------------------------------------------------------------------------
  !
  ! ... C wrapper for library interface to the Pwscf
  USE ISO_C_BINDING
  !
  IMPLICIT NONE
  !
  INTEGER (kind=C_INT), VALUE :: lib_comm, nim, npt, npl, nta, nbn, ndg
  INTEGER (kind=C_INT), INTENT(OUT) :: retval
  CHARACTER (kind=C_CHAR), INTENT(IN) :: infile(*)
  INTEGER  :: i, lib_comm_, nim_, npt_, npl_, nta_, nbn_, ndg_, retval_
  CHARACTER(LEN=80)  :: infile_
  !
  ! ... Copy C data types to Fortran data types
  lib_comm_ = lib_comm
  nim_ = nim
  npt_ = npt
  npl_ = npl
  nta_ = nta
  nbn_ = nbn
  ndg_ = ndg
  retval = 0
  infile_ = ' '
  !
  ! ... Copying a string from C to Fortran is a bit ugly.
  DO i=1,80
      IF (infile(i) == C_NULL_CHAR) EXIT
      infile_ = TRIM(infile_) // infile(i)
  END DO
  !
  CALL f_ht_pwscf_drv(lib_comm_,nim_,npt_,npl_,nta_,nbn_,ndg_,retval_,infile_)
  retval = retval_
  !
END SUBROUTINE ht_pw_drv
!
!----------------------------------------------------------------------------
SUBROUTINE f_ht_pwscf_drv(lib_comm,nim,npt,npl,nta,nbn,ndg,retval,infile)
  !----------------------------------------------------------------------------
  !
  ! ... Library interface to the Plane Wave Self-Consistent Field code
  !
  USE environment,       ONLY : environment_start
  USE mp_global,         ONLY : mp_startup
  USE read_input,        ONLY : read_input_file
  USE command_line_options, ONLY: set_command_line
  USE parallel_include
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: lib_comm, nim, npt, npl, nta, nbn, ndg
  INTEGER, INTENT(INOUT) :: retval
  CHARACTER(LEN=80)      :: infile
  !
  INTEGER :: me, num, ierr
  CALL MPI_COMM_SIZE(lib_comm,num,ierr)
  IF (ierr /= MPI_SUCCESS) THEN
      CALL MPI_ERROR_STRING(ierr, infile, 80, retval)
      PRINT*,'MPI Error: ', infile
      STOP 100
  END IF
  CALL MPI_COMM_RANK(lib_comm,me,ierr)
  IF (me == 0) THEN
      PRINT*, 'Calling PW library interface with these flags:'
      PRINT*, 'communicator index: ', lib_comm
      PRINT*, 'communicator size:  ', num
      PRINT*, 'nimage: ', nim
      PRINT*, 'npool:  ', npl
      PRINT*, 'ntaskg: ', nta
      PRINT*, 'nband:  ', nbn
      PRINT*, 'ndiag:  ', ndg
      PRINT*, 'input:  "',TRIM(infile),'"'
  END IF
  !
  CALL set_command_line( nimage=nim, npool=npl, ntg=nta, &
      nband=nbn, ndiag=ndg )
  CALL mp_startup ( my_world_comm=lib_comm )
  CALL ht_environment_start ( 'PWSCF' )
  !
  CALL read_input_file ('PW', infile )
  !
  ! ... Perform actual calculation
  !
  CALL run_pwscf  ( retval )
  !
  CALL stop_run( retval )
  !
END SUBROUTINE f_ht_pwscf_drv

SUBROUTINE ht_environment_start( code )

  USE kinds, ONLY: DP
  USE io_files, ONLY: crash_file, nd_nmbr
  USE io_global, ONLY: stdout, meta_ionode
  USE mp_world,  ONLY: nproc
  USE mp_images, ONLY: me_image, my_image_id, root_image, nimage, nproc_image
  USE mp_pools,  ONLY: npool
  USE mp_bands,  ONLY: ntask_groups, nproc_bgrp, nbgrp
  USE global_version, ONLY: version_number, svn_revision
  USE environment

    CHARACTER(LEN=*), INTENT(IN) :: code
    LOGICAL           :: exst, debug = .false.
    CHARACTER(LEN=80) :: code_version, uname
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    INTEGER :: ios, crashunit
    INTEGER, EXTERNAL :: find_free_unit

    ! ... The Intel compiler allocates a lot of stack space
    ! ... Stack limit is often small, thus causing SIGSEGV and crash
    ! ... One may use "ulimit -s unlimited" but it doesn't always work
    ! ... The following call does the same and always works
    !
#if defined(__INTEL_COMPILER)
    CALL remove_stack_limit ( )
#endif
    ! ... use ".FALSE." to disable all clocks except the total cpu time clock
    ! ... use ".TRUE."  to enable clocks

    CALL init_clocks( .TRUE. )
    CALL start_clock( TRIM(code) )

    code_version = TRIM (code) // " v." // TRIM (version_number)
    IF ( TRIM (svn_revision) /= "unknown" ) code_version = &
         TRIM (code_version) // " (svn rev. " // TRIM (svn_revision) // ")"

    ! ... for compatibility with PWSCF

#if defined(__MPI)
    nd_nmbr = TRIM ( int_to_char( me_image+1 ))
#else
    nd_nmbr = ' '
#endif

    IF( meta_ionode ) THEN

       ! ...  search for file CRASH and delete it

       INQUIRE( FILE=TRIM(crash_file), EXIST=exst )
       IF( exst ) THEN
          crashunit = find_free_unit()
          OPEN( UNIT=crashunit, FILE=TRIM(crash_file), STATUS='OLD',IOSTAT=ios )
          IF (ios==0) THEN
             CLOSE( UNIT=crashunit, STATUS='DELETE', IOSTAT=ios )
          ELSE
             WRITE(stdout,'(5x,"Remark: CRASH file could not be deleted")')
          END IF
       END IF
    END IF

    ! ... one processor per image (other than meta_ionode)
    ! ... or, for debugging purposes, all processors,
    ! ... open their own standard output file

    uname = 'out.' // trim(int_to_char( my_image_id )) // '_' // trim(int_to_char( me_image))
    OPEN ( unit = stdout, file = TRIM(uname),status='unknown')
    !
    ! Need to modify "Modules/environment.f90" to compile and print these info
    CALL opening_message( code_version )
    CALL compilation_info ( )
#if defined(__MPI)
    CALL parallel_info ( )
#else
    CALL serial_info()
#endif

    !
END SUBROUTINE
