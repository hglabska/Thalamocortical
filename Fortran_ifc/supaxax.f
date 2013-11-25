               PROGRAM IB


        PARAMETER (numcomp_supaxax     = 59)
        PARAMETER ( num_supaxax = 1,num_suppyrRS =0)

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
     &  V_supaxax    (numcomp_supaxax,   num_supaxax)


       double precision::
     &  curr_supaxax     (numcomp_supaxax,   num_supaxax)

! define membrane and Ca state variables that must be passed
! to subroutines

       real*8  chi_supaxax(numcomp_supaxax,num_supaxax)
       real*8  mnaf_supaxax(numcomp_supaxax,num_supaxax),
     & mnap_supaxax(numcomp_supaxax,num_supaxax),
     x hnaf_supaxax(numcomp_supaxax,num_supaxax),
     x mkdr_supaxax(numcomp_supaxax,num_supaxax),
     x mka_supaxax(numcomp_supaxax,num_supaxax),
     x hka_supaxax(numcomp_supaxax,num_supaxax),
     x mk2_supaxax(numcomp_supaxax,num_supaxax), 
     x hk2_supaxax(numcomp_supaxax,num_supaxax),
     x mkm_supaxax(numcomp_supaxax,num_supaxax),
     x mkc_supaxax(numcomp_supaxax,num_supaxax),
     x mkahp_supaxax(numcomp_supaxax,num_supaxax),
     x mcat_supaxax(numcomp_supaxax,num_supaxax),
     x hcat_supaxax(numcomp_supaxax,num_supaxax),
     x mcal_supaxax(numcomp_supaxax,num_supaxax),
     x mar_supaxax(numcomp_supaxax,num_supaxax)



       double precision
     &    ranvec_supaxax    (num_supaxax),
     &    seed /137.d0/
       integer, parameter :: totaxgj_supaxax    = 0 
       integer, parameter :: totaxgj_mix      =  0 ! decr. antidr. bursts in IB 
       integer gjtable_supaxax   (totaxgj_supaxax,4)
       integergjtable_mix     (totaxgj_mix,4)
       double precision, parameter :: gapcon_supaxax    = 4.d-3

c Define arrays for distal axon voltages which will be shared
c between nodes.
         double precision::
     &  distal_axon_supaxax    (num_supaxax),
     &  distal_axon_global    (14000) !14000 = 14 x num_supaxax
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
     &    vax_supaxax (num_supaxax), vax_suppyrRS(num_suppyrRS)

         double precision::
     &  outtime_supaxax    (5000, num_supaxax)
     

         INTEGER
     &  outctr_supaxax    (num_supaxax)


       REAL*8  time1, time2, time, timtot,end_stim,start_stim, tmp1, tmp2
       REAL*8 presyntime, delta, dexparg, dexparg1, dexparg2
       INTEGER thisno, display /0/, O
       REAL*8 z, z1, z2, outrcd(20)

       REAL dtime, t_time(2), e_time ! e_time is dtime output

       double precision::
     & gAMPA_supaxax     (numcomp_supaxax,   num_supaxax),
     & gNMDA_supaxax     (numcomp_supaxax,   num_supaxax),
     & gGABA_A_supaxax     (numcomp_supaxax,   num_supaxax)

         time1 = 0.d0	!	 used to be =gettime()
c Setup tables for calculating exponentials
          call dexptablesmall_setup (dexptablesmall)
           do j = 1, numcomp_supaxax
c       WRITE(6,919) j
c919     FORMAT(I3,X)
       mnaf_supaxax(j,num_supaxax) = 0
       mnap_supaxax(j,num_supaxax) = 0
       hnaf_supaxax(j,num_supaxax) = 0
       mkdr_supaxax(j,num_supaxax) = 0
       mka_supaxax(j,num_supaxax) = 0
       hka_supaxax(j,num_supaxax) =0 
       mk2_supaxax(j,num_supaxax) =0
       hk2_supaxax(j,num_supaxax) =0
       mkm_supaxax(j,num_supaxax)=0
       mkc_supaxax(j,num_supaxax) =0
       mkahp_supaxax(j,num_supaxax)=0
       mcat_supaxax(j,num_supaxax)=0
       hcat_supaxax(j,num_supaxax)=0
       mcal_supaxax(j,num_supaxax)=0
       mar_supaxax(j,num_supaxax)=0
       gAMPA_supaxax     (j,   num_supaxax) =0
       gNMDA_supaxax     (j,   num_supaxax) =0
       gGABA_A_supaxax     (j,   num_supaxax) =0
       curr_supaxax     (j,   num_supaxax) = 0
       V_supaxax    (j,   num_supaxax) =0
       chi_supaxax(numcomp_supaxax,num_supaxax) =0
		end do

           do j = 1, num_supaxax   
        outtime_supaxax(1,j)               = -1.d5
	outctr_supaxax(j)		  = 0
           end do ! j

          call dexptablebig_setup(dexptablebig)
          timtot = 300.0d0 ! debug time does test 0.05 interval 
	ranvec_supaxax(1) =0

c begin define input current
	end_stim = 200
	start_stim = 100.0d0
	O = 0
1000    O = O + 1
        time = time + dt
        if (time.gt.timtot) goto 2000
	curr_supaxax     (1,   num_supaxax)  = -0.1
 	if (time.lt.1000.d0) then
		if ( time.gt.start_stim.and.time.le.end_stim) then
    	curr_supaxax     (1,   num_supaxax)  = 0.4
c		tmp1 = 1*(-dexp((start_stim-time)*2))
c		tmp2 = 1*(dexp((start_stim-time)/20))
c		tmp1 = (1-exp((start_stim-time)/10))
c		tmp2 = (exp((start_stim-time)/20))
c        curr_supaxax(43,num_supaxax)= 3*tmp1 *  tmp2
c        WRITE(6,927)time,curr_supaxax(43,num_supaxax),tmp1,tmp2,start_stim
c927     FORMAT(E14.7,X,E14.7,X,E14.7,X, E14.7,X,E14.7,X)
		else 
			if (time.gt.end_stim) then
			curr_supaxax(1,num_supaxax) = 0
			endif
		 endif
	 endif

c end define input current

       CALL INTEGRATE_supaxax (O, time, num_supaxax,
     &    V_supaxax, curr_supaxax,
     & gAMPA_supaxax, gNMDA_supaxax , gGABA_A_supaxax ,
     & Mg, 
     & gapcon_supaxax,totaxgj_supaxax,gjtable_supaxax, dt,
     &  chi_supaxax,mnaf_supaxax,mnap_supaxax,
     &  hnaf_supaxax,mkdr_supaxax,mka_supaxax,
     &  hka_supaxax,mk2_supaxax,hk2_supaxax,
     &  mkm_supaxax,mkc_supaxax,mkahp_supaxax,
     &  mcat_supaxax,hcat_supaxax,mcal_supaxax,
     &  mar_supaxax)

        outrcd( 1) = time
        outrcd( 2) = v_supaxax   (1,1)
        outrcd( 3) = v_supaxax   (39   ,1)
        outrcd( 4) = v_supaxax   (43,1)
         z = 0.d0
          do i = 1, num_supaxax   
           z = z - v_supaxax(1,i)
          end do
        outrcd( 5) = z / dble(num_supaxax   ) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_supaxax   
           z = z + gAMPA_supaxax   (i,1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_supaxax   
           z = z + gNMDA_supaxax   (i,1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_supaxax   
           z = z + gGABA_A_supaxax   (i,1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_supaxax   (1,3)
        outrcd(10) = v_supaxax   (1,2)
        outrcd(11) =  curr_supaxax(1,   num_supaxax)
c        outrcd(12) = field_2mm_supaxax

      OPEN(17,FILE='GROUCHO110.supaxax')
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





