               PROGRAM IB


        PARAMETER (numcomp_nontuftRS     = 50)
        PARAMETER ( num_nontuftRS = 1)

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
     &  V_nontuftRS    (numcomp_nontuftRS,   num_nontuftRS)


       double precision::
     &  curr_nontuftRS     (numcomp_nontuftRS,   num_nontuftRS)

! define membrane and Ca state variables that must be passed
! to subroutines

       real*8  chi_nontuftRS(numcomp_nontuftRS,num_nontuftRS)
       real*8  mnaf_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     & mnap_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x hnaf_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mkdr_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mka_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x hka_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mk2_nontuftRS(numcomp_nontuftRS,num_nontuftRS), 
     x hk2_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mkm_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mkc_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mkahp_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mcat_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x hcat_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mcal_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mar_nontuftRS(numcomp_nontuftRS,num_nontuftRS)



       double precision
     &    ranvec_nontuftRS    (num_nontuftRS),
     &    seed /137.d0/
       integer, parameter :: totaxgj_nontuftRS    = 0 
       integer, parameter :: totaxgj_tuft      =  0 ! decr. antidr. bursts in IB 
       integer gjtable_nontuftRS   (totaxgj_nontuftRS,4)
       integergjtable_tuft     (totaxgj_tuft,4)
       double precision, parameter :: gapcon_nontuftRS    = 4.d-3

c Define arrays for distal axon voltages which will be shared
c between nodes.
         double precision::
     &  distal_axon_nontuftRS    (num_nontuftRS),
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


         double precision::
     &  outtime_nontuftRS    (5000, num_nontuftRS)
     

         INTEGER
     &  outctr_nontuftRS    (num_nontuftRS)


       REAL*8  time1, time2, time, timtot,end_stim,start_stim, tmp1, tmp2
       REAL*8 presyntime, delta, dexparg, dexparg1, dexparg2
       INTEGER thisno, display /0/, O
       REAL*8 z, z1, z2, outrcd(20)

       REAL dtime, t_time(2), e_time ! e_time is dtime output

       double precision::
     & gAMPA_nontuftRS     (numcomp_nontuftRS,   num_nontuftRS),
     & gNMDA_nontuftRS     (numcomp_nontuftRS,   num_nontuftRS),
     & gGABA_A_nontuftRS     (numcomp_nontuftRS,   num_nontuftRS)

         time1 = 0.d0	!	 used to be =gettime()
c Setup tables for calculating exponentials
          call dexptablesmall_setup (dexptablesmall)
           do j = 1, numcomp_nontuftRS
c       WRITE(6,919) j
c919     FORMAT(I3,X)
       mnaf_nontuftRS(j,num_nontuftRS) = 0
       mnap_nontuftRS(j,num_nontuftRS) = 0
       hnaf_nontuftRS(j,num_nontuftRS) = 0
       mkdr_nontuftRS(j,num_nontuftRS) = 0
       mka_nontuftRS(j,num_nontuftRS) = 0
       hka_nontuftRS(j,num_nontuftRS) =0 
       mk2_nontuftRS(j,num_nontuftRS) =0
       hk2_nontuftRS(j,num_nontuftRS) =0
       mkm_nontuftRS(j,num_nontuftRS)=0
       mkc_nontuftRS(j,num_nontuftRS) =0
       mkahp_nontuftRS(j,num_nontuftRS)=0
       mcat_nontuftRS(j,num_nontuftRS)=0
       hcat_nontuftRS(j,num_nontuftRS)=0
       mcal_nontuftRS(j,num_nontuftRS)=0
       mar_nontuftRS(j,num_nontuftRS)=0
       gAMPA_nontuftRS     (j,   num_nontuftRS) =0
       gNMDA_nontuftRS     (j,   num_nontuftRS) =0
       gGABA_A_nontuftRS     (j,   num_nontuftRS) =0
       curr_nontuftRS     (j,   num_nontuftRS) = 0
       V_nontuftRS    (j,   num_nontuftRS) =0
       chi_nontuftRS(numcomp_nontuftRS,num_nontuftRS) =0
		end do

           do j = 1, num_nontuftRS   
        outtime_nontuftRS(1,j)               = -1.d5
	outctr_nontuftRS(j)		  = 0
           end do ! j

          call dexptablebig_setup(dexptablebig)
          timtot = 1000.0d0 ! debug time does test 0.05 interval 
	ranvec_nontuftRS(1) =0

c begin define input current
	end_stim = 700.0d0
	start_stim = 300.0d0
	O = 0
1000    O = O + 1
        time = time + dt
        if (time.gt.timtot) goto 2000
	curr_nontuftRS     (1,   num_nontuftRS) = -0.25
	curr_nontuftRS     (43,   num_nontuftRS)  = 0
 	if (time.lt.1000.d0) then
		if ( time.gt.start_stim.and.time.le.end_stim) then
c		tmp1 = 1*(-dexp((start_stim-time)*2))
c		tmp2 = 1*(dexp((start_stim-time)/20))
c		tmp1 = (1-exp((start_stim-time)/10))
c		tmp2 = (exp((start_stim-time)/20))
c        curr_nontuftRS(43,num_nontuftRS)= 3*3*tmp1 *  tmp2 
c        WRITE(6,927)time,curr_nontuftRS(43,num_nontuftRS),tmp1,tmp2,start_stim
c927     FORMAT(E14.7,X,E14.7,X,E14.7,X, E14.7,X,E14.7,X)
		curr_nontuftRS(1,num_nontuftRS) = 0.5
		else 
			if (time.gt.end_stim) then
			curr_nontuftRS(1,num_nontuftRS) = 0
			endif
		 endif
	 endif

c end define input current

       CALL INTEGRATE_nontuftRS (O, time, num_nontuftRS,
     &    V_nontuftRS, curr_nontuftRS,
     & gAMPA_nontuftRS, gNMDA_nontuftRS , gGABA_A_nontuftRS ,
     & Mg, 
     & gapcon_nontuftRS,totaxgj_nontuftRS,gjtable_nontuftRS, dt,
     &  chi_nontuftRS,mnaf_nontuftRS,mnap_nontuftRS,
     &  hnaf_nontuftRS,mkdr_nontuftRS,mka_nontuftRS,
     &  hka_nontuftRS,mk2_nontuftRS,hk2_nontuftRS,
     &  mkm_nontuftRS,mkc_nontuftRS,mkahp_nontuftRS,
     &  mcat_nontuftRS,hcat_nontuftRS,mcal_nontuftRS,
     &  mar_nontuftRS,field_1mm_nontuftRS,field_2mm_nontuftRS)

        outrcd( 1) = time
        outrcd( 2) = v_nontuftRS   (1,1)
        outrcd( 3) = v_nontuftRS   (39   ,1)
        outrcd( 4) = v_nontuftRS   (43,1)
         z = 0.d0
          do i = 1, num_nontuftRS   
           z = z - v_nontuftRS(1,i)
          end do
        outrcd( 5) = z / dble(num_nontuftRS   ) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_nontuftRS   
           z = z + gAMPA_nontuftRS   (i,1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_nontuftRS   
           z = z + gNMDA_nontuftRS   (i,1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_nontuftRS   
           z = z + gGABA_A_nontuftRS   (i,1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_nontuftRS   (1,3)
        outrcd(10) = v_nontuftRS   (1,2)
        outrcd(11) =  curr_nontuftRS(43,   num_nontuftRS)
c        outrcd(12) = field_2mm_nontuftRS

      OPEN(17,FILE='GROUCHO110.nontuftRS')
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





