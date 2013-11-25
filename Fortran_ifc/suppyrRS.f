               PROGRAM IB


        PARAMETER (numcomp_suppyrRS     = 74)
        PARAMETER ( num_suppyrRS = 1,num_suppyrFRB =0)

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
     &  V_suppyrRS    (numcomp_suppyrRS,   num_suppyrRS)


       double precision::
     &  curr_suppyrRS     (numcomp_suppyrRS,   num_suppyrRS)

! define membrane and Ca state variables that must be passed
! to subroutines

       real*8  chi_suppyrRS(numcomp_suppyrRS,num_suppyrRS)
       real*8  mnaf_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     & mnap_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x hnaf_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mkdr_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mka_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x hka_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mk2_suppyrRS(numcomp_suppyrRS,num_suppyrRS), 
     x hk2_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mkm_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mkc_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mkahp_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mcat_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x hcat_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mcal_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mar_suppyrRS(numcomp_suppyrRS,num_suppyrRS)



       double precision
     &    ranvec_suppyrRS    (num_suppyrRS),
     &    seed /137.d0/
       integer, parameter :: totaxgj_suppyrRS    = 0 
       integer, parameter :: totaxgj_mix      =  0 ! decr. antidr. bursts in IB 
       integer gjtable_suppyrRS   (totaxgj_suppyrRS,4)
       integergjtable_mix     (totaxgj_mix,4)
       double precision, parameter :: gapcon_suppyrRS    = 4.d-3

c Define arrays for distal axon voltages which will be shared
c between nodes.
         double precision::
     &  distal_axon_suppyrRS    (num_suppyrRS),
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
     &    vax_suppyrRS (num_suppyrRS), vax_suppyrFRB(num_suppyrFRB)

         double precision::
     &  outtime_suppyrRS    (5000, num_suppyrRS)
     

         INTEGER
     &  outctr_suppyrRS    (num_suppyrRS)


       REAL*8  time1, time2, time, timtot,end_stim,start_stim, tmp1, tmp2
       REAL*8 presyntime, delta, dexparg, dexparg1, dexparg2
       INTEGER thisno, display /0/, O
       REAL*8 z, z1, z2, outrcd(20)

       REAL dtime, t_time(2), e_time ! e_time is dtime output

       double precision::
     & gAMPA_suppyrRS     (numcomp_suppyrRS,   num_suppyrRS),
     & gNMDA_suppyrRS     (numcomp_suppyrRS,   num_suppyrRS),
     & gGABA_A_suppyrRS     (numcomp_suppyrRS,   num_suppyrRS)

         time1 = 0.d0	!	 used to be =gettime()
c Setup tables for calculating exponentials
          call dexptablesmall_setup (dexptablesmall)
           do j = 1, numcomp_suppyrRS
c       WRITE(6,919) j
c919     FORMAT(I3,X)
       mnaf_suppyrRS(j,num_suppyrRS) = 0
       mnap_suppyrRS(j,num_suppyrRS) = 0
       hnaf_suppyrRS(j,num_suppyrRS) = 0
       mkdr_suppyrRS(j,num_suppyrRS) = 0
       mka_suppyrRS(j,num_suppyrRS) = 0
       hka_suppyrRS(j,num_suppyrRS) =0 
       mk2_suppyrRS(j,num_suppyrRS) =0
       hk2_suppyrRS(j,num_suppyrRS) =0
       mkm_suppyrRS(j,num_suppyrRS)=0
       mkc_suppyrRS(j,num_suppyrRS) =0
       mkahp_suppyrRS(j,num_suppyrRS)=0
       mcat_suppyrRS(j,num_suppyrRS)=0
       hcat_suppyrRS(j,num_suppyrRS)=0
       mcal_suppyrRS(j,num_suppyrRS)=0
       mar_suppyrRS(j,num_suppyrRS)=0
       gAMPA_suppyrRS     (j,   num_suppyrRS) =0
       gNMDA_suppyrRS     (j,   num_suppyrRS) =0
       gGABA_A_suppyrRS     (j,   num_suppyrRS) =0
       curr_suppyrRS     (j,   num_suppyrRS) = 0
       V_suppyrRS    (j,   num_suppyrRS) =0
       chi_suppyrRS(numcomp_suppyrRS,num_suppyrRS) =0
		end do

           do j = 1, num_suppyrRS   
        outtime_suppyrRS(1,j)               = -1.d5
	outctr_suppyrRS(j)		  = 0
           end do ! j

          call dexptablebig_setup(dexptablebig)
          timtot = 1000.0d0 ! debug time does test 0.05 interval 
	ranvec_suppyrRS(1) =0

c begin define input current
	end_stim = 700
	start_stim = 400.0d0
	O = 0
1000    O = O + 1
        time = time + dt
        if (time.gt.timtot) goto 2000
	curr_suppyrRS     (1,   num_suppyrRS)  = -0.4
 	if (time.lt.1000.d0) then
		if ( time.gt.start_stim.and.time.le.end_stim) then
    	curr_suppyrRS     (1,   num_suppyrRS)  = 0.75
c		tmp1 = 1*(-dexp((start_stim-time)*2))
c		tmp2 = 1*(dexp((start_stim-time)/20))
c		tmp1 = (1-exp((start_stim-time)/10))
c		tmp2 = (exp((start_stim-time)/20))
c        curr_suppyrRS(43,num_suppyrRS)= 3*tmp1 *  tmp2
c        WRITE(6,927)time,curr_suppyrRS(43,num_suppyrRS),tmp1,tmp2,start_stim
c927     FORMAT(E14.7,X,E14.7,X,E14.7,X, E14.7,X,E14.7,X)
		else 
			if (time.gt.end_stim) then
			curr_suppyrRS(1,num_suppyrRS) = 0
			endif
		 endif
	 endif

c end define input current

       CALL INTEGRATE_suppyrRS (O, time, num_suppyrRS,
     &    V_suppyrRS, curr_suppyrRS,
     & gAMPA_suppyrRS, gNMDA_suppyrRS , gGABA_A_suppyrRS ,
     & Mg, 
     & gapcon_suppyrRS,totaxgj_suppyrRS,gjtable_suppyrRS, dt,
     & totaxgj_mix  , gjtable_mix  , 0  ,
     & vax_suppyrFRB  ,
     &  chi_suppyrRS,mnaf_suppyrRS,mnap_suppyrRS,
     &  hnaf_suppyrRS,mkdr_suppyrRS,mka_suppyrRS,
     &  hka_suppyrRS,mk2_suppyrRS,hk2_suppyrRS,
     &  mkm_suppyrRS,mkc_suppyrRS,mkahp_suppyrRS,
     &  mcat_suppyrRS,hcat_suppyrRS,mcal_suppyrRS,
     &  mar_suppyrRS,field_1mm_suppyrRS,field_2mm_suppyrRS)

        outrcd( 1) = time
        outrcd( 2) = v_suppyrRS   (1,1)
        outrcd( 3) = v_suppyrRS   (39   ,1)
        outrcd( 4) = v_suppyrRS   (43,1)
         z = 0.d0
          do i = 1, num_suppyrRS   
           z = z - v_suppyrRS(1,i)
          end do
        outrcd( 5) = z / dble(num_suppyrRS   ) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_suppyrRS   
           z = z + gAMPA_suppyrRS   (i,1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_suppyrRS   
           z = z + gNMDA_suppyrRS   (i,1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_suppyrRS   
           z = z + gGABA_A_suppyrRS   (i,1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_suppyrRS   (1,3)
        outrcd(10) = v_suppyrRS   (1,2)
        outrcd(11) =  curr_suppyrRS(1,   num_suppyrRS)
c        outrcd(12) = field_2mm_suppyrRS

      OPEN(17,FILE='GROUCHO110.suppyrRS')
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





