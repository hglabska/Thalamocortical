               PROGRAM IB


        PARAMETER (numcomp_nRT     = 59)
        PARAMETER ( num_nRT = 1,num_suppyrRS =0)

c Assorted parameters
         double precision, parameter :: dt = 0.002d0
         double precision, parameter :: Mg = 1.50
! Castro-Alamancos J Physiol, disinhib. neocortex in vitro, uses
! Mg = 1.3
         double precision, parameter :: NMDA_saturation_fact
!    &                                   = 5.d0
     &                                   = 80.d0

        double precision:: dexptablesmall(0:5000)
        double precision::  dexptablebig(0:10000)
c NMDA conductance developed on one postsynaptic compartment,
c from one type of presynaptic cell, can be at most this
c factor x unitary conductance
c UNFORTUNATELY, with this scheme,if one NMDA cond. set to 0
c on a cell type, all NMDA conductances will be forced to 0
c on that cell type...

       double precision, parameter :: thal_cort_delay = 1.d0
       double precision, parameter :: cort_thal_delay = 5.d0
       integer, parameter :: how_often = 50
! how_often defines how many time steps between synaptic conductance
! updates, and between broadcastings of axonal voltages.
       double precision, parameter :: axon_refrac_time = 1.5d0

       double precision::
     &  V_nRT    (numcomp_nRT,   num_nRT)


       double precision::
     &  curr_nRT     (numcomp_nRT,   num_nRT)

! define membrane and Ca state variables that must be passed
! to subroutines

       real*8  chi_nRT(numcomp_nRT,num_nRT)
       real*8  mnaf_nRT(numcomp_nRT,num_nRT),
     & mnap_nRT(numcomp_nRT,num_nRT),
     x hnaf_nRT(numcomp_nRT,num_nRT),
     x mkdr_nRT(numcomp_nRT,num_nRT),
     x mka_nRT(numcomp_nRT,num_nRT),
     x hka_nRT(numcomp_nRT,num_nRT),
     x mk2_nRT(numcomp_nRT,num_nRT), 
     x hk2_nRT(numcomp_nRT,num_nRT),
     x mkm_nRT(numcomp_nRT,num_nRT),
     x mkc_nRT(numcomp_nRT,num_nRT),
     x mkahp_nRT(numcomp_nRT,num_nRT),
     x mcat_nRT(numcomp_nRT,num_nRT),
     x hcat_nRT(numcomp_nRT,num_nRT),
     x mcal_nRT(numcomp_nRT,num_nRT),
     x mar_nRT(numcomp_nRT,num_nRT)



       double precision
     &    ranvec_nRT    (num_nRT),
     &    seed /137.d0/
       integer, parameter :: totaxgj_nRT    = 0 
       integer, parameter :: totaxgj_mix      =  0 ! decr. antidr. bursts in IB 
       integer gjtable_nRT   (totaxgj_nRT,4)
       integergjtable_mix     (totaxgj_mix,4)
       double precision, parameter :: gapcon_nRT    = 4.d-3

c Define arrays for distal axon voltages which will be shared
c between nodes.
         double precision::
     &  distal_axon_nRT    (num_nRT),
     &  distal_axon_global    (14000) !14000 = 14 x num_nRT
! distal_axon_global will be concatenation of individual
! distal_axon vectors       
! positions 1      -  1000  suppyrRS axons
!           1001   -  2000  suppyrFRB axons
!           2001   -  3000  supbask
!           3001   -  4000  supaxax
!           4001   -  5000  nRT
!           5001   -  6000  spinstell
!           6001   -  7000  tuftIB
!           7001   -  8000  tuftRS
!           8001   -  9000  nontuftRS
!           9001   - 10000  deepbask
!          10001   - 11000  deepaxax
!          11001   - 12000  deepLTS
!          12001   - 13000  TCR
!          13001   - 14000  nRT


! define arrays for axonal voltges, needed for mixed gj
         double precision ::
     &    vax_nRT (num_nRT)

         double precision::
     &  outtime_nRT    (5000, num_nRT)
     

         INTEGER
     &  outctr_nRT    (num_nRT)


       REAL*8  time1, time2, time, timtot,end_stim,start_stim, tmp1, tmp2
       REAL*8 presyntime, delta, dexparg, dexparg1, dexparg2
       INTEGER thisno, display /0/, O
       REAL*8 z, z1, z2, outrcd(20)

       REAL dtime, t_time(2), e_time ! e_time is dtime output

       double precision::
     & gAMPA_nRT     (numcomp_nRT,   num_nRT),
     & gNMDA_nRT     (numcomp_nRT,   num_nRT),
     & gGABA_A_nRT     (numcomp_nRT,   num_nRT)

         time1 = 0.d0	!	 used to be =gettime()
c Setup tables for calculating exponentials
          call dexptablesmall_setup (dexptablesmall)
           do j = 1, numcomp_nRT
c       WRITE(6,919) j
c919     FORMAT(I3,X)
       mnaf_nRT(j,num_nRT) = 0
       mnap_nRT(j,num_nRT) = 0
       hnaf_nRT(j,num_nRT) = 0
       mkdr_nRT(j,num_nRT) = 0
       mka_nRT(j,num_nRT) = 0
       hka_nRT(j,num_nRT) =0 
       mk2_nRT(j,num_nRT) =0
       hk2_nRT(j,num_nRT) =0
       mkm_nRT(j,num_nRT)=0
       mkc_nRT(j,num_nRT) =0
       mkahp_nRT(j,num_nRT)=0
       mcat_nRT(j,num_nRT)=0
       hcat_nRT(j,num_nRT)=0
       mcal_nRT(j,num_nRT)=0
       mar_nRT(j,num_nRT)=0
       gAMPA_nRT     (j,   num_nRT) =0
       gNMDA_nRT     (j,   num_nRT) =0
       gGABA_A_nRT     (j,   num_nRT) =0
       curr_nRT     (j,   num_nRT) = 0
       V_nRT    (j,   num_nRT) =0
       chi_nRT(numcomp_nRT,num_nRT) =0
		end do

           do j = 1, num_nRT   
        outtime_nRT(1,j)               = -1.d5
	outctr_nRT(j)		  = 0
           end do ! j

          call dexptablebig_setup(dexptablebig)
          timtot = 200.0d0 ! debug time does test 0.05 interval 
	ranvec_nRT(1) =0

c begin define input current
	end_stim = 200
	start_stim = 50.0d0
	O = 0
1000    O = O + 1
        time = time + dt
        if (time.gt.timtot) goto 2000
	curr_nRT     (1,   num_nRT)  = -0.4
 	if (time.lt.1000.d0) then
		if ( time.gt.start_stim.and.time.le.end_stim) then
    	curr_nRT     (1,   num_nRT)  = 0.4
c		tmp1 = 1*(-dexp((start_stim-time)*2))
c		tmp2 = 1*(dexp((start_stim-time)/20))
c		tmp1 = (1-exp((start_stim-time)/10))
c		tmp2 = (exp((start_stim-time)/20))
c        curr_nRT(43,num_nRT)= 3*tmp1 *  tmp2
c        WRITE(6,927)time,curr_nRT(43,num_nRT),tmp1,tmp2,start_stim
c927     FORMAT(E14.7,X,E14.7,X,E14.7,X, E14.7,X,E14.7,X)
		else 
			if (time.gt.end_stim) then
			curr_nRT(1,num_nRT) = 0
			endif
		 endif
	 endif

c end define input current
       CALL INTEGRATE_nRT (O, time, num_nRT,
     &    V_nRT, curr_nRT,
     & gAMPA_nRT, gNMDA_nRT , gGABA_A_nRT ,
     & Mg, 
     & gapcon_nRT,totaxgj_nRT,gjtable_nRT, dt,
     &  chi_nRT,mnaf_nRT,mnap_nRT,
     &  hnaf_nRT,mkdr_nRT,mka_nRT,
     &  hka_nRT,mk2_nRT,hk2_nRT,
     &  mkm_nRT,mkc_nRT,mkahp_nRT,
     &  mcat_nRT,hcat_nRT,mcal_nRT,
     &  mar_nRT)

        outrcd( 1) = time
        outrcd( 2) = v_nRT   (1,1)
        outrcd( 3) = v_nRT   (59   ,1)
        outrcd( 4) = v_nRT   (43,1)
         z = 0.d0
          do i = 1, num_nRT   
           z = z - v_nRT(1,i)
          end do
        outrcd( 5) = z / dble(num_nRT   ) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_nRT   
           z = z + gAMPA_nRT   (i,1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_nRT   
           z = z + gNMDA_nRT   (i,1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_nRT   
           z = z + gGABA_A_nRT   (i,1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_nRT   (1,3)
        outrcd(10) = v_nRT   (1,2)
        outrcd(11) =  curr_nRT(1,   num_nRT)
c        outrcd(12) = field_2mm_nRT

      OPEN(17,FILE='GROUCHO110.nRT')
      WRITE (17,FMT='(19F10.4)') (OUTRCD(I),I=1,11)

        goto 1000

2000    CONTINUE
        e_time = dtime(t_time)
        time2 = t_time(1)	! used to be = gettime()
        print *,'elapsed:',e_time,', user:',t_time(1),', sys:',t_time(2)
        write(6,3434) time2 - time1
3434    format(' elapsed time = ',f8.0,' secs')

!        call mpi_finalize (info)
             END





