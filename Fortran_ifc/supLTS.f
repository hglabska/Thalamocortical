               PROGRAM IB


        PARAMETER (numcomp_supLTS     = 59)
        PARAMETER ( num_supLTS = 1,num_suppyrRS =0)

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
     &  V_supLTS    (numcomp_supLTS,   num_supLTS)


       double precision::
     &  curr_supLTS     (numcomp_supLTS,   num_supLTS)

! define membrane and Ca state variables that must be passed
! to subroutines

       real*8  chi_supLTS(numcomp_supLTS,num_supLTS)
       real*8  mnaf_supLTS(numcomp_supLTS,num_supLTS),
     & mnap_supLTS(numcomp_supLTS,num_supLTS),
     x hnaf_supLTS(numcomp_supLTS,num_supLTS),
     x mkdr_supLTS(numcomp_supLTS,num_supLTS),
     x mka_supLTS(numcomp_supLTS,num_supLTS),
     x hka_supLTS(numcomp_supLTS,num_supLTS),
     x mk2_supLTS(numcomp_supLTS,num_supLTS), 
     x hk2_supLTS(numcomp_supLTS,num_supLTS),
     x mkm_supLTS(numcomp_supLTS,num_supLTS),
     x mkc_supLTS(numcomp_supLTS,num_supLTS),
     x mkahp_supLTS(numcomp_supLTS,num_supLTS),
     x mcat_supLTS(numcomp_supLTS,num_supLTS),
     x hcat_supLTS(numcomp_supLTS,num_supLTS),
     x mcal_supLTS(numcomp_supLTS,num_supLTS),
     x mar_supLTS(numcomp_supLTS,num_supLTS)



       double precision
     &    ranvec_supLTS    (num_supLTS),
     &    seed /137.d0/
       integer, parameter :: totaxgj_supLTS    = 0 
       integer, parameter :: totaxgj_mix      =  0 ! decr. antidr. bursts in IB 
       integer gjtable_supLTS   (totaxgj_supLTS,4)
       integergjtable_mix     (totaxgj_mix,4)
       double precision, parameter :: gapcon_supLTS    = 4.d-3

c Define arrays for distal axon voltages which will be shared
c between nodes.
         double precision::
     &  distal_axon_supLTS    (num_supLTS),
     &  distal_axon_global    (14000) !14000 = 14 x num_supLTS
! distal_axon_global will be concatenation of individual
! distal_axon vectors       
! positions 1      -  1000  suppyrRS axons
!           1001   -  2000  suppyrFRB axons
!           2001   -  3000  supbask
!           3001   -  4000  supaxax
!           4001   -  5000  supLTS
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
     &    vax_supLTS (num_supLTS), vax_suppyrRS(num_suppyrRS)

         double precision::
     &  outtime_supLTS    (5000, num_supLTS)
     

         INTEGER
     &  outctr_supLTS    (num_supLTS)


       REAL*8  time1, time2, time, timtot,end_stim,start_stim, tmp1, tmp2
       REAL*8 presyntime, delta, dexparg, dexparg1, dexparg2
       INTEGER thisno, display /0/, O
       REAL*8 z, z1, z2, outrcd(20)

       REAL dtime, t_time(2), e_time ! e_time is dtime output

       double precision::
     & gAMPA_supLTS     (numcomp_supLTS,   num_supLTS),
     & gNMDA_supLTS     (numcomp_supLTS,   num_supLTS),
     & gGABA_A_supLTS     (numcomp_supLTS,   num_supLTS)

         time1 = 0.d0	!	 used to be =gettime()
c Setup tables for calculating exponentials
          call dexptablesmall_setup (dexptablesmall)
           do j = 1, numcomp_supLTS
c       WRITE(6,919) j
c919     FORMAT(I3,X)
       mnaf_supLTS(j,num_supLTS) = 0
       mnap_supLTS(j,num_supLTS) = 0
       hnaf_supLTS(j,num_supLTS) = 0
       mkdr_supLTS(j,num_supLTS) = 0
       mka_supLTS(j,num_supLTS) = 0
       hka_supLTS(j,num_supLTS) =0 
       mk2_supLTS(j,num_supLTS) =0
       hk2_supLTS(j,num_supLTS) =0
       mkm_supLTS(j,num_supLTS)=0
       mkc_supLTS(j,num_supLTS) =0
       mkahp_supLTS(j,num_supLTS)=0
       mcat_supLTS(j,num_supLTS)=0
       hcat_supLTS(j,num_supLTS)=0
       mcal_supLTS(j,num_supLTS)=0
       mar_supLTS(j,num_supLTS)=0
       gAMPA_supLTS     (j,   num_supLTS) =0
       gNMDA_supLTS     (j,   num_supLTS) =0
       gGABA_A_supLTS     (j,   num_supLTS) =0
       curr_supLTS     (j,   num_supLTS) = 0
       V_supLTS    (j,   num_supLTS) =0
       chi_supLTS(numcomp_supLTS,num_supLTS) =0
		end do

           do j = 1, num_supLTS   
        outtime_supLTS(1,j)               = -1.d5
	outctr_supLTS(j)		  = 0
           end do ! j

          call dexptablebig_setup(dexptablebig)
          timtot = 300.0d0 ! debug time does test 0.05 interval 
	ranvec_supLTS(1) =0

c begin define input current
	end_stim = 200
	start_stim = 100.0d0
	O = 0
1000    O = O + 1
        time = time + dt
        if (time.gt.timtot) goto 2000
	curr_supLTS     (1,   num_supLTS)  = -0.15
 	if (time.lt.1000.d0) then
		if ( time.gt.start_stim.and.time.le.end_stim) then
    	curr_supLTS     (1,   num_supLTS)  = 0.4
c		tmp1 = 1*(-dexp((start_stim-time)*2))
c		tmp2 = 1*(dexp((start_stim-time)/20))
c		tmp1 = (1-exp((start_stim-time)/10))
c		tmp2 = (exp((start_stim-time)/20))
c        curr_supLTS(43,num_supLTS)= 3*tmp1 *  tmp2
c        WRITE(6,927)time,curr_supLTS(43,num_supLTS),tmp1,tmp2,start_stim
c927     FORMAT(E14.7,X,E14.7,X,E14.7,X, E14.7,X,E14.7,X)
		else 
			if (time.gt.end_stim) then
			curr_supLTS(1,num_supLTS) = 0
			endif
		 endif
	 endif

c end define input current

       CALL INTEGRATE_supLTS (O, time, num_supLTS,
     &    V_supLTS, curr_supLTS,
     & gAMPA_supLTS, gNMDA_supLTS , gGABA_A_supLTS ,
     & Mg, 
     & gapcon_supLTS,totaxgj_supLTS,gjtable_supLTS, dt,
     &  chi_supLTS,mnaf_supLTS,mnap_supLTS,
     &  hnaf_supLTS,mkdr_supLTS,mka_supLTS,
     &  hka_supLTS,mk2_supLTS,hk2_supLTS,
     &  mkm_supLTS,mkc_supLTS,mkahp_supLTS,
     &  mcat_supLTS,hcat_supLTS,mcal_supLTS,
     &  mar_supLTS)

        outrcd( 1) = time
        outrcd( 2) = v_supLTS   (1,1)
        outrcd( 3) = v_supLTS   (39   ,1)
        outrcd( 4) = v_supLTS   (43,1)
         z = 0.d0
          do i = 1, num_supLTS   
           z = z - v_supLTS(1,i)
          end do
        outrcd( 5) = z / dble(num_supLTS   ) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_supLTS   
           z = z + gAMPA_supLTS   (i,1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_supLTS   
           z = z + gNMDA_supLTS   (i,1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_supLTS   
           z = z + gGABA_A_supLTS   (i,1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_supLTS   (1,3)
        outrcd(10) = v_supLTS   (1,2)
        outrcd(11) =  curr_supLTS(1,   num_supLTS)
c        outrcd(12) = field_2mm_supLTS

      OPEN(17,FILE='GROUCHO110.supLTS')
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





