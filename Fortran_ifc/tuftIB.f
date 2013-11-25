               PROGRAM IB


        PARAMETER (numcomp_tuftIB     = 61)
        PARAMETER ( num_tuftIB = 1,num_tuftRS =0)

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
     &  V_tuftIB    (numcomp_tuftIB,   num_tuftIB)


       double precision::
     &  curr_tuftIB     (numcomp_tuftIB,   num_tuftIB)

! define membrane and Ca state variables that must be passed
! to subroutines

       real*8  chi_tuftIB(numcomp_tuftIB,num_tuftIB)
       real*8  mnaf_tuftIB(numcomp_tuftIB,num_tuftIB),
     & mnap_tuftIB(numcomp_tuftIB,num_tuftIB),
     x hnaf_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mkdr_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mka_tuftIB(numcomp_tuftIB,num_tuftIB),
     x hka_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mk2_tuftIB(numcomp_tuftIB,num_tuftIB), 
     x hk2_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mkm_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mkc_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mkahp_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mcat_tuftIB(numcomp_tuftIB,num_tuftIB),
     x hcat_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mcal_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mar_tuftIB(numcomp_tuftIB,num_tuftIB)



       double precision
     &    ranvec_tuftIB    (num_tuftIB),
     &    seed /137.d0/
       integer, parameter :: totaxgj_tuftIB    = 0 
       integer, parameter :: totaxgj_tuft      =  0 ! decr. antidr. bursts in IB 
       integer gjtable_tuftIB   (totaxgj_tuftIB,4)
       integergjtable_tuft     (totaxgj_tuft,4)
       double precision, parameter :: gapcon_tuftIB    = 4.d-3

c Define arrays for distal axon voltages which will be shared
c between nodes.
         double precision::
     &  distal_axon_tuftIB    (num_tuftIB),
     &  distal_axon_global    (14000) !14000 = 14 x num_suppyrRS
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
     &    vax_tuftIB (num_tuftIB), vax_tuftRS (num_tuftRS)

         double precision::
     &  outtime_tuftIB    (5000, num_tuftIB)
     

         INTEGER
     &  outctr_tuftIB    (num_tuftIB)


       REAL*8  time1, time2, time, timtot,end_stim,start_stim, tmp1, tmp2
       REAL*8 presyntime, delta, dexparg, dexparg1, dexparg2
       INTEGER thisno, display /0/, O
       REAL*8 z, z1, z2, outrcd(20)

       REAL dtime, t_time(2), e_time ! e_time is dtime output

       double precision::
     & gAMPA_tuftIB     (numcomp_tuftIB,   num_tuftIB),
     & gNMDA_tuftIB     (numcomp_tuftIB,   num_tuftIB),
     & gGABA_A_tuftIB     (numcomp_tuftIB,   num_tuftIB)

         time1 = 0.d0	!	 used to be =gettime()
c Setup tables for calculating exponentials
          call dexptablesmall_setup (dexptablesmall)
           do j = 1, numcomp_tuftIB
c       WRITE(6,919) j
c919     FORMAT(I3,X)
       mnaf_tuftIB(j,num_tuftIB) = 0
       mnap_tuftIB(j,num_tuftIB) = 0
       hnaf_tuftIB(j,num_tuftIB) = 0
       mkdr_tuftIB(j,num_tuftIB) = 0
       mka_tuftIB(j,num_tuftIB) = 0
       hka_tuftIB(j,num_tuftIB) =0 
       mk2_tuftIB(j,num_tuftIB) =0
       hk2_tuftIB(j,num_tuftIB) =0
       mkm_tuftIB(j,num_tuftIB)=0
       mkc_tuftIB(j,num_tuftIB) =0
       mkahp_tuftIB(j,num_tuftIB)=0
       mcat_tuftIB(j,num_tuftIB)=0
       hcat_tuftIB(j,num_tuftIB)=0
       mcal_tuftIB(j,num_tuftIB)=0
       mar_tuftIB(j,num_tuftIB)=0
       gAMPA_tuftIB     (j,   num_tuftIB) =0
       gNMDA_tuftIB     (j,   num_tuftIB) =0
       gGABA_A_tuftIB     (j,   num_tuftIB) =0
       curr_tuftIB     (j,   num_tuftIB) = 0
       V_tuftIB    (j,   num_tuftIB) =0
       chi_tuftIB(numcomp_tuftIB,num_tuftIB) =0
		end do

           do j = 1, num_tuftIB   
        outtime_tuftIB(1,j)               = -1.d5
	outctr_tuftIB(j)		  = 0
           end do ! j

          call dexptablebig_setup(dexptablebig)
          timtot = 1000.0d0 ! debug time does test 0.05 interval 
	ranvec_tuftIB(1) =0

c begin define input current
	end_stim = timtot
	start_stim = 505.0d0
	O = 0
1000    O = O + 1
        time = time + dt
        if (time.gt.timtot) goto 2000
	curr_tuftIB     (1,   num_tuftIB) = 0.5
	curr_tuftIB     (43,   num_tuftIB)  = 0
 	if (time.lt.1000.d0) then
		if ( time.gt.start_stim.and.time.le.end_stim) then
c		tmp1 = 1*(-dexp((start_stim-time)*2))
c		tmp2 = 1*(dexp((start_stim-time)/20))
		tmp1 = (1-exp((start_stim-time)/10))
		tmp2 = (exp((start_stim-time)/20))
        curr_tuftIB(43,num_tuftIB)= 3*3*tmp1 *  tmp2 
c        WRITE(6,927)time,curr_tuftIB(43,num_tuftIB),tmp1,tmp2,start_stim
c927     FORMAT(E14.7,X,E14.7,X,E14.7,X, E14.7,X,E14.7,X)
		else 
			curr_tuftIB(43,num_tuftIB) = 0
		 endif
	 endif

c end define input current

       CALL INTEGRATE_tuftIB (O, time, num_tuftIB,
     &    V_tuftIB, curr_tuftIB,
     & gAMPA_tuftIB, gNMDA_tuftIB , gGABA_A_tuftIB ,
     & Mg, 
     & gapcon_tuftIB,totaxgj_tuftIB,gjtable_tuftIB, dt,
     & totaxgj_tuft  , gjtable_tuft  , 0  ,
     & vax_tuftRS   ,
     &  chi_tuftIB,mnaf_tuftIB,mnap_tuftIB,
     &  hnaf_tuftIB,mkdr_tuftIB,mka_tuftIB,
     &  hka_tuftIB,mk2_tuftIB,hk2_tuftIB,
     &  mkm_tuftIB,mkc_tuftIB,mkahp_tuftIB,
     &  mcat_tuftIB,hcat_tuftIB,mcal_tuftIB,
     &  mar_tuftIB,field_1mm_tuftIB,field_2mm_tuftIB)

        outrcd( 1) = time
        outrcd( 2) = v_tuftIB   (1,1)
        outrcd( 3) = v_tuftIB   (39   ,1)
        outrcd( 4) = v_tuftIB   (43,1)
         z = 0.d0
          do i = 1, num_tuftIB   
           z = z - v_tuftIB(1,i)
          end do
        outrcd( 5) = z / dble(num_tuftIB   ) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_tuftIB   
           z = z + gAMPA_tuftIB   (i,1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_tuftIB   
           z = z + gNMDA_tuftIB   (i,1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_tuftIB   
           z = z + gGABA_A_tuftIB   (i,1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_tuftIB   (1,3)
        outrcd(10) = v_tuftIB   (1,2)
        outrcd(11) =  curr_tuftIB(43,   num_tuftIB)
c        outrcd(12) = field_2mm_tuftIB

      OPEN(17,FILE='GROUCHO110.tuftIB')
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





