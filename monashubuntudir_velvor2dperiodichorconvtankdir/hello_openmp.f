      PROGRAM HELLO
       !     USE omp_lib
       IMPLICIT NONE
       INTEGER OMP_GET_MAX_THREADS
       INTEGER OMP_GET_NUM_THREADS
       INTEGER OMP_GET_THREAD_NUM
       
       write(6,"(a, i3)") " OpenMP max threads: ", OMP_GET_MAX_THREADS()
!$OMP PARALLEL
       write(6,"(2(a,i3))") " OpenMP: N_threads = ",     
     &   OMP_GET_NUM_THREADS()," thread = ", OMP_GET_THREAD_NUM()
!$OMP END PARALLEL
      END PROGRAM

