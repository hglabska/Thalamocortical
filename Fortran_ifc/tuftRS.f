               PROGRAM IB


        PARAMETER (numcomp_tuftRS     = 61)
        PARAMETER ( num_tuftRS = 1,num_tuftIB =0)

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
     &  V_tuftRS    (numcomp_tuftRS,   num_tuftRS)


       double precision::
     &  curr_tuftRS     (numcomp_tuftRS,   num_tuftRS)

! define membrane and Ca state variables that must be passed
! to subroutines

       real*8  chi_tuftRS(numcomp_tuftRS,num_tuftRS)
       real*8  mnaf_tuftRS(numcomp_tuftRS,num_tuftRS),
     & mnap_tuftRS(numcomp_tuftRS,num_tuftRS),
     x hnaf_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mkdr_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mka_tuftRS(numcomp_tuftRS,num_tuftRS),
     x hka_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mk2_tuftRS(numcomp_tuftRS,num_tuftRS), 
     x hk2_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mkm_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mkc_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mkahp_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mcat_tuftRS(numcomp_tuftRS,num_tuftRS),
     x hcat_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mcal_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mar_tuftRS(numcomp_tuftRS,num_tuftRS)



       double precision
     &    ranvec_tuftRS    (num_tuftRS),
     &    seed /137.d0/
       integer, parameter :: totaxgj_tuftRS    = 0 
       integer, parameter :: totaxgj_tuft      =  0 ! decr. antidr. bursts in IB 
       integer gjtable_tuftRS   (totaxgj_tuftRS,4)
       integergjtable_tuft     (totaxgj_tuft,4)
       double precision, parameter :: gapcon_tuftRS    = 4.d-3

c Define arrays for distal axon voltages which will be shared
c between nodes.
         double precision::
     &  distal_axon_tuftRS    (num_tuftRS),
     &  distal_axon_global    (14000) !14000 = 14 x num_suppyrRS
! distal_axon_global will be concatenation of individual
! distal_axon vectors       
! positions 1      -  1000  suppyrRS axons
!           1001   -  2000  suppyrFRB axons
!           2001   -  3000  supbask
!           3001   -  4000  supaxax
!           4001   -  5000  supLTS
!           5001   -  6000  spinstell
!           6001   -  7000  tuftRS
!           7001   -  8000  tuftRS
!           8001   -  9000  nontuftRS
!           9001   - 10000  deepbask
!          10001   - 11000  deepaxax
!          11001   - 12000  deepLTS
!          12001   - 13000  TCR
!          13001   - 14000  nRT

! define arrays for axonal voltges, needed for mixed gj
         double precision ::
     &    vax_tuftRS (num_tuftRS), vax_tuftIB (num_tuftIB)

         double precision::
     &  outtime_tuftRS    (5000, num_tuftRS)
     

         INTEGER
     &  outctr_tuftRS    (num_tuftRS)


       REAL*8  time1, time2, time, timtot,end_stim,start_stim, tmp1, tmp2
       REAL*8 presyntime, delta, dexparg, dexparg1, dexparg2
       INTEGER thisno, display /0/, O
       REAL*8 z, z1, z2, outrcd(20)

       REAL dtime, t_time(2), e_time ! e_time is dtime output

       double precision::
     & gAMPA_tuftRS     (numcomp_tuftRS,   num_tuftRS),
     & gNMDA_tuftRS     (numcomp_tuftRS,   num_tuftRS),
     & gGABA_A_tuftRS     (numcomp_tuftRS,   num_tuftRS)

         time1 = 0.d0	!	 used to be =gettime()
c Setup tables for calculating exponentials
          call dexptablesmall_setup (dexptablesmall)
           do j = 1, numcomp_tuftRS
c       WRITE(6,919) j
c919     FORMAT(I3,X)
       mnaf_tuftRS(j,num_tuftRS) = 0
       mnap_tuftRS(j,num_tuftRS) = 0
       hnaf_tuftRS(j,num_tuftRS) = 0
       mkdr_tuftRS(j,num_tuftRS) = 0
       mka_tuftRS(j,num_tuftRS) = 0
       hka_tuftRS(j,num_tuftRS) =0 
       mk2_tuftRS(j,num_tuftRS) =0
       hk2_tuftRS(j,num_tuftRS) =0
       mkm_tuftRS(j,num_tuftRS)=0
       mkc_tuftRS(j,num_tuftRS) =0
       mkahp_tuftRS(j,num_tuftRS)=0
       mcat_tuftRS(j,num_tuftRS)=0
       hcat_tuftRS(j,num_tuftRS)=0
       mcal_tuftRS(j,num_tuftRS)=0
       mar_tuftRS(j,num_tuftRS)=0
       gAMPA_tuftRS     (j,   num_tuftRS) =0
       gNMDA_tuftRS     (j,   num_tuftRS) =0
       gGABA_A_tuftRS     (j,   num_tuftRS) =0
       curr_tuftRS     (j,   num_tuftRS) = 0
       V_tuftRS    (j,   num_tuftRS) =0
       chi_tuftRS(numcomp_tuftRS,num_tuftRS) =0
		end do

           do j = 1, num_tuftRS   
        outtime_tuftRS(1,j)               = -1.d5
	outctr_tuftRS(j)		  = 0
           end do ! j

          call dexptablebig_setup(dexptablebig)
          timtot = 1000.0d0 ! debug time does test 0.05 interval 
	ranvec_tuftRS(1) =0

c begin define input current
	end_stim = 700.0d0
	start_stim = 300.0d0
	O = 0
1000    O = O + 1
        time = time + dt
        if (time.gt.timtot) goto 2000
	curr_tuftRS     (1,   num_tuftRS) = -0.4
	curr_tuftRS     (43,   num_tuftRS)  = 0
 	if (time.lt.1000.d0) then
		if ( time.gt.start_stim.and.time.le.end_stim) then
c		tmp1 = 1*(-dexp((start_stim-time)*2))
c		tmp2 = 1*(dexp((start_stim-time)/20))
c		tmp1 = (1-exp((start_stim-time)/10))
c		tmp2 = (exp((start_stim-time)/20))
c        curr_tuftRS(43,num_tuftRS)= 3*3*tmp1 *  tmp2 
c        WRITE(6,927)time,curr_tuftRS(43,num_tuftRS),tmp1,tmp2,start_stim
c927     FORMAT(E14.7,X,E14.7,X,E14.7,X, E14.7,X,E14.7,X)
		curr_tuftRS(1,num_tuftRS) = 0.8
		else 
			if (time.gt.end_stim) then
			curr_tuftRS(1,num_tuftRS) = 0
			endif
		 endif
	 endif

c end define input current

       CALL INTEGRATE_tuftRS (O, time, num_tuftRS,
     &    V_tuftRS, curr_tuftRS,
     & gAMPA_tuftRS, gNMDA_tuftRS , gGABA_A_tuftRS ,
     & Mg, 
     & gapcon_tuftRS,totaxgj_tuftRS,gjtable_tuftRS, dt,
     & totaxgj_tuft  , gjtable_tuft  , 0  ,
     & vax_tuftRS   ,
     &  chi_tuftRS,mnaf_tuftRS,mnap_tuftRS,
     &  hnaf_tuftRS,mkdr_tuftRS,mka_tuftRS,
     &  hka_tuftRS,mk2_tuftRS,hk2_tuftRS,
     &  mkm_tuftRS,mkc_tuftRS,mkahp_tuftRS,
     &  mcat_tuftRS,hcat_tuftRS,mcal_tuftRS,
     &  mar_tuftRS,field_1mm_tuftRS,field_2mm_tuftRS)

        outrcd( 1) = time
        outrcd( 2) = v_tuftRS   (1,1)
        outrcd( 3) = v_tuftRS   (39   ,1)
        outrcd( 4) = v_tuftRS   (43,1)
         z = 0.d0
          do i = 1, num_tuftRS   
           z = z - v_tuftRS(1,i)
          end do
        outrcd( 5) = z / dble(num_tuftRS   ) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_tuftRS   
           z = z + gAMPA_tuftRS   (i,1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_tuftRS   
           z = z + gNMDA_tuftRS   (i,1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_tuftRS   
           z = z + gGABA_A_tuftRS   (i,1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_tuftRS   (1,3)
        outrcd(10) = v_tuftRS   (1,2)
        outrcd(11) =  curr_tuftRS(43,   num_tuftRS)
c        outrcd(12) = field_2mm_tuftRS

      OPEN(17,FILE='GROUCHO110.tuftRS')
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





