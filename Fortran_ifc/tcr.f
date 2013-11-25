               PROGRAM IB


        PARAMETER (numcomp_tcr     = 137)
        PARAMETER ( num_tcr = 1)

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
     &  V_tcr    (numcomp_tcr,   num_tcr)


       double precision::
     &  curr_tcr     (numcomp_tcr,   num_tcr)

! define membrane and Ca state variables that must be passed
! to subroutines

       real*8  chi_tcr(numcomp_tcr,num_tcr)
       real*8  mnaf_tcr(numcomp_tcr,num_tcr),
     & mnap_tcr(numcomp_tcr,num_tcr),
     x hnaf_tcr(numcomp_tcr,num_tcr),
     x mkdr_tcr(numcomp_tcr,num_tcr),
     x mka_tcr(numcomp_tcr,num_tcr),
     x hka_tcr(numcomp_tcr,num_tcr),
     x mk2_tcr(numcomp_tcr,num_tcr), 
     x hk2_tcr(numcomp_tcr,num_tcr),
     x mkm_tcr(numcomp_tcr,num_tcr),
     x mkc_tcr(numcomp_tcr,num_tcr),
     x mkahp_tcr(numcomp_tcr,num_tcr),
     x mcat_tcr(numcomp_tcr,num_tcr),
     x hcat_tcr(numcomp_tcr,num_tcr),
     x mcal_tcr(numcomp_tcr,num_tcr),
     x mar_tcr(numcomp_tcr,num_tcr)



       double precision
     &    ranvec_tcr    (num_tcr),
     &    seed /137.d0/
       integer, parameter :: totaxgj_tcr    = 0 
       integer, parameter :: totaxgj_mix      =  0 ! decr. antidr. bursts in IB 
       integer gjtable_tcr   (totaxgj_tcr,4)
       integergjtable_mix     (totaxgj_mix,4)
       double precision, parameter :: gapcon_tcr    = 4.d-3

c Define arrays for distal axon voltages which will be shared
c between nodes.
         double precision::
     &  distal_axon_tcr    (num_tcr),
     &  distal_axon_global    (14000) !14000 = 14 x num_tcr
! distal_axon_global will be concatenation of individual
! distal_axon vectors       
! positions 1      -  1000  suppyrRS axons
!           1001   -  2000  suppyrFRB axons
!           2001   -  3000  supbask
!           3001   -  4000  supaxax
!           4001   -  5000  tcr
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
     &    vax_tcr (num_tcr)

         double precision::
     &  outtime_tcr    (5000, num_tcr)
     

         INTEGER
     &  outctr_tcr    (num_tcr)


       REAL*8  time1, time2, time, timtot,end_stim,start_stim, tmp1, tmp2
       REAL*8 presyntime, delta, dexparg, dexparg1, dexparg2
       INTEGER thisno, display /0/, O
       REAL*8 z, z1, z2, outrcd(20)

       REAL dtime, t_time(2), e_time ! e_time is dtime output

       double precision::
     & gAMPA_tcr     (numcomp_tcr,   num_tcr),
     & gNMDA_tcr     (numcomp_tcr,   num_tcr),
     & gGABA_A_tcr     (numcomp_tcr,   num_tcr)

         time1 = 0.d0	!	 used to be =gettime()
c Setup tables for calculating exponentials
          call dexptablesmall_setup (dexptablesmall)
           do j = 1, numcomp_tcr
c       WRITE(6,919) j
c919     FORMAT(I3,X)
       mnaf_tcr(j,num_tcr) = 0
       mnap_tcr(j,num_tcr) = 0
       hnaf_tcr(j,num_tcr) = 0
       mkdr_tcr(j,num_tcr) = 0
       mka_tcr(j,num_tcr) = 0
       hka_tcr(j,num_tcr) =0 
       mk2_tcr(j,num_tcr) =0
       hk2_tcr(j,num_tcr) =0
       mkm_tcr(j,num_tcr)=0
       mkc_tcr(j,num_tcr) =0
       mkahp_tcr(j,num_tcr)=0
       mcat_tcr(j,num_tcr)=0
       hcat_tcr(j,num_tcr)=0
       mcal_tcr(j,num_tcr)=0
       mar_tcr(j,num_tcr)=0
       gAMPA_tcr     (j,   num_tcr) =0
       gNMDA_tcr     (j,   num_tcr) =0
       gGABA_A_tcr     (j,   num_tcr) =0
       curr_tcr     (j,   num_tcr) = 0
       V_tcr    (j,   num_tcr) =0
       chi_tcr(numcomp_tcr,num_tcr) =0
		end do

           do j = 1, num_tcr   
        outtime_tcr(1,j)               = -1.d5
	outctr_tcr(j)		  = 0
           end do ! j

          call dexptablebig_setup(dexptablebig)
          timtot = 350.0d0 ! debug time does test 0.05 interval 
	ranvec_tcr(1) =0

c begin define input current
	end_stim = 350
	start_stim = 300.0d0
	O = 0
1000    O = O + 1
        time = time + dt
        if (time.gt.timtot) goto 2000
	curr_tcr     (1,   num_tcr)  = -0.9
 	if (time.lt.1000.d0) then
		if ( time.gt.start_stim.and.time.le.end_stim) then
    	curr_tcr     (1,   num_tcr)  = -0.3
c		tmp1 = 1*(-dexp((start_stim-time)*2))
c		tmp2 = 1*(dexp((start_stim-time)/20))
c		tmp1 = (1-exp((start_stim-time)/10))
c		tmp2 = (exp((start_stim-time)/20))
c        curr_tcr(43,num_tcr)= 3*tmp1 *  tmp2
c        WRITE(6,927)time,curr_tcr(43,num_tcr),tmp1,tmp2,start_stim
c927     FORMAT(E14.7,X,E14.7,X,E14.7,X, E14.7,X,E14.7,X)
		else 
			if (time.gt.end_stim) then
			curr_tcr(1,num_tcr) = 0
			endif
		 endif
	 endif

c end define input current
       CALL INTEGRATE_tcr (O, time, num_tcr,
     &    V_tcr, curr_tcr,
     & gAMPA_tcr, gNMDA_tcr , gGABA_A_tcr ,
     & Mg, 
     & gapcon_tcr,totaxgj_tcr,gjtable_tcr, dt,
     &  chi_tcr,mnaf_tcr,mnap_tcr,
     &  hnaf_tcr,mkdr_tcr,mka_tcr,
     &  hka_tcr,mk2_tcr,hk2_tcr,
     &  mkm_tcr,mkc_tcr,mkahp_tcr,
     &  mcat_tcr,hcat_tcr,mcal_tcr,
     &  mar_tcr)

        outrcd( 1) = time
        outrcd( 2) = v_tcr   (1,1)
        outrcd( 3) = v_tcr   (137   ,1)
        outrcd( 4) = v_tcr   (43,1)
         z = 0.d0
          do i = 1, num_tcr   
           z = z - v_tcr(1,i)
          end do
        outrcd( 5) = z / dble(num_tcr   ) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_tcr   
           z = z + gAMPA_tcr   (i,1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_tcr   
           z = z + gNMDA_tcr   (i,1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_tcr   
           z = z + gGABA_A_tcr   (i,1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_tcr   (1,3)
        outrcd(10) = v_tcr   (1,2)
        outrcd(11) =  curr_tcr(1,   num_tcr)
c        outrcd(12) = field_2mm_tcr

      OPEN(17,FILE='GROUCHO110.tcr')
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





