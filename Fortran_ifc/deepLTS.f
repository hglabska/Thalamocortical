               PROGRAM IB


        PARAMETER (numcomp_deepLTS     = 59)
        PARAMETER ( num_deepLTS = 1,num_suppyrRS =0)

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
     &  V_deepLTS    (numcomp_deepLTS,   num_deepLTS)


       double precision::
     &  curr_deepLTS     (numcomp_deepLTS,   num_deepLTS)

! define membrane and Ca state variables that must be passed
! to subroutines

       real*8  chi_deepLTS(numcomp_deepLTS,num_deepLTS)
       real*8  mnaf_deepLTS(numcomp_deepLTS,num_deepLTS),
     & mnap_deepLTS(numcomp_deepLTS,num_deepLTS),
     x hnaf_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mkdr_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mka_deepLTS(numcomp_deepLTS,num_deepLTS),
     x hka_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mk2_deepLTS(numcomp_deepLTS,num_deepLTS), 
     x hk2_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mkm_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mkc_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mkahp_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mcat_deepLTS(numcomp_deepLTS,num_deepLTS),
     x hcat_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mcal_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mar_deepLTS(numcomp_deepLTS,num_deepLTS)



       double precision
     &    ranvec_deepLTS    (num_deepLTS),
     &    seed /137.d0/
       integer, parameter :: totaxgj_deepLTS    = 0 
       integer, parameter :: totaxgj_mix      =  0 ! decr. antidr. bursts in IB 
       integer gjtable_deepLTS   (totaxgj_deepLTS,4)
       integergjtable_mix     (totaxgj_mix,4)
       double precision, parameter :: gapcon_deepLTS    = 4.d-3

c Define arrays for distal axon voltages which will be shared
c between nodes.
         double precision::
     &  distal_axon_deepLTS    (num_deepLTS),
     &  distal_axon_global    (14000) !14000 = 14 x num_deepLTS
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
     &    vax_deepLTS (num_deepLTS), vax_suppyrRS(num_suppyrRS)

         double precision::
     &  outtime_deepLTS    (5000, num_deepLTS)
     

         INTEGER
     &  outctr_deepLTS    (num_deepLTS)


       REAL*8  time1, time2, time, timtot,end_stim,start_stim, tmp1, tmp2
       REAL*8 presyntime, delta, dexparg, dexparg1, dexparg2
       INTEGER thisno, display /0/, O
       REAL*8 z, z1, z2, outrcd(20)

       REAL dtime, t_time(2), e_time ! e_time is dtime output

       double precision::
     & gAMPA_deepLTS     (numcomp_deepLTS,   num_deepLTS),
     & gNMDA_deepLTS     (numcomp_deepLTS,   num_deepLTS),
     & gGABA_A_deepLTS     (numcomp_deepLTS,   num_deepLTS)

         time1 = 0.d0	!	 used to be =gettime()
c Setup tables for calculating exponentials
          call dexptablesmall_setup (dexptablesmall)
           do j = 1, numcomp_deepLTS
c       WRITE(6,919) j
c919     FORMAT(I3,X)
       mnaf_deepLTS(j,num_deepLTS) = 0
       mnap_deepLTS(j,num_deepLTS) = 0
       hnaf_deepLTS(j,num_deepLTS) = 0
       mkdr_deepLTS(j,num_deepLTS) = 0
       mka_deepLTS(j,num_deepLTS) = 0
       hka_deepLTS(j,num_deepLTS) =0 
       mk2_deepLTS(j,num_deepLTS) =0
       hk2_deepLTS(j,num_deepLTS) =0
       mkm_deepLTS(j,num_deepLTS)=0
       mkc_deepLTS(j,num_deepLTS) =0
       mkahp_deepLTS(j,num_deepLTS)=0
       mcat_deepLTS(j,num_deepLTS)=0
       hcat_deepLTS(j,num_deepLTS)=0
       mcal_deepLTS(j,num_deepLTS)=0
       mar_deepLTS(j,num_deepLTS)=0
       gAMPA_deepLTS     (j,   num_deepLTS) =0
       gNMDA_deepLTS     (j,   num_deepLTS) =0
       gGABA_A_deepLTS     (j,   num_deepLTS) =0
       curr_deepLTS     (j,   num_deepLTS) = 0
       V_deepLTS    (j,   num_deepLTS) =0
       chi_deepLTS(numcomp_deepLTS,num_deepLTS) =0
		end do

           do j = 1, num_deepLTS   
        outtime_deepLTS(1,j)               = -1.d5
	outctr_deepLTS(j)		  = 0
           end do ! j

          call dexptablebig_setup(dexptablebig)
          timtot = 300.0d0 ! debug time does test 0.05 interval 
	ranvec_deepLTS(1) =0

c begin define input current
	end_stim = 200
	start_stim = 100.0d0
	O = 0
1000    O = O + 1
        time = time + dt
        if (time.gt.timtot) goto 2000
	curr_deepLTS     (1,   num_deepLTS)  = -0.15
 	if (time.lt.1000.d0) then
		if ( time.gt.start_stim.and.time.le.end_stim) then
    	curr_deepLTS     (1,   num_deepLTS)  = 0.4
c		tmp1 = 1*(-dexp((start_stim-time)*2))
c		tmp2 = 1*(dexp((start_stim-time)/20))
c		tmp1 = (1-exp((start_stim-time)/10))
c		tmp2 = (exp((start_stim-time)/20))
c        curr_deepLTS(43,num_deepLTS)= 3*tmp1 *  tmp2
c        WRITE(6,927)time,curr_deepLTS(43,num_deepLTS),tmp1,tmp2,start_stim
c927     FORMAT(E14.7,X,E14.7,X,E14.7,X, E14.7,X,E14.7,X)
		else 
			if (time.gt.end_stim) then
			curr_deepLTS(1,num_deepLTS) = 0
			endif
		 endif
	 endif

c end define input current

       CALL INTEGRATE_deepLTS (O, time, num_deepLTS,
     &    V_deepLTS, curr_deepLTS,
     & gAMPA_deepLTS, gNMDA_deepLTS , gGABA_A_deepLTS ,
     & Mg, 
     & gapcon_deepLTS,totaxgj_deepLTS,gjtable_deepLTS, dt,
     &  chi_deepLTS,mnaf_deepLTS,mnap_deepLTS,
     &  hnaf_deepLTS,mkdr_deepLTS,mka_deepLTS,
     &  hka_deepLTS,mk2_deepLTS,hk2_deepLTS,
     &  mkm_deepLTS,mkc_deepLTS,mkahp_deepLTS,
     &  mcat_deepLTS,hcat_deepLTS,mcal_deepLTS,
     &  mar_deepLTS)

        outrcd( 1) = time
        outrcd( 2) = v_deepLTS   (1,1)
        outrcd( 3) = v_deepLTS   (39   ,1)
        outrcd( 4) = v_deepLTS   (43,1)
         z = 0.d0
          do i = 1, num_deepLTS   
           z = z - v_deepLTS(1,i)
          end do
        outrcd( 5) = z / dble(num_deepLTS   ) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_deepLTS   
           z = z + gAMPA_deepLTS   (i,1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_deepLTS   
           z = z + gNMDA_deepLTS   (i,1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_deepLTS   
           z = z + gGABA_A_deepLTS   (i,1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_deepLTS   (1,3)
        outrcd(10) = v_deepLTS   (1,2)
        outrcd(11) =  curr_deepLTS(1,   num_deepLTS)
c        outrcd(12) = field_2mm_deepLTS

      OPEN(17,FILE='GROUCHO110.deepLTS')
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





