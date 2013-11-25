               PROGRAM IB


        PARAMETER (numcomp_spinstell     = 59)
        PARAMETER ( num_spinstell = 1)

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
     &  V_spinstell    (numcomp_spinstell,   num_spinstell)


       double precision::
     &  curr_spinstell     (numcomp_spinstell,   num_spinstell)

! define membrane and Ca state variables that must be passed
! to subroutines

       real*8  chi_spinstell(numcomp_spinstell,num_spinstell)
       real*8  mnaf_spinstell(numcomp_spinstell,num_spinstell),
     & mnap_spinstell(numcomp_spinstell,num_spinstell),
     x hnaf_spinstell(numcomp_spinstell,num_spinstell),
     x mkdr_spinstell(numcomp_spinstell,num_spinstell),
     x mka_spinstell(numcomp_spinstell,num_spinstell),
     x hka_spinstell(numcomp_spinstell,num_spinstell),
     x mk2_spinstell(numcomp_spinstell,num_spinstell), 
     x hk2_spinstell(numcomp_spinstell,num_spinstell),
     x mkm_spinstell(numcomp_spinstell,num_spinstell),
     x mkc_spinstell(numcomp_spinstell,num_spinstell),
     x mkahp_spinstell(numcomp_spinstell,num_spinstell),
     x mcat_spinstell(numcomp_spinstell,num_spinstell),
     x hcat_spinstell(numcomp_spinstell,num_spinstell),
     x mcal_spinstell(numcomp_spinstell,num_spinstell),
     x mar_spinstell(numcomp_spinstell,num_spinstell)



       double precision
     &    ranvec_spinstell    (num_spinstell),
     &    seed /137.d0/
       integer, parameter :: totaxgj_spinstell    = 0 
       integer, parameter :: totaxgj_mix      =  0 ! decr. antidr. bursts in IB 
       integer gjtable_spinstell   (totaxgj_spinstell,4)
       integergjtable_mix     (totaxgj_mix,4)
       double precision, parameter :: gapcon_spinstell    = 4.d-3

c Define arrays for distal axon voltages which will be shared
c between nodes.
         double precision::
     &  distal_axon_spinstell    (num_spinstell),
     &  distal_axon_global    (14000) !14000 = 14 x num_spinstell
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
!          10001   - 11000  spinstell
!          11001   - 12000  deepLTS
!          12001   - 13000  TCR
!          13001   - 14000  nRT


! define arrays for axonal voltges, needed for mixed gj
         double precision ::
     &    vax_spinstell (num_spinstell)

         double precision::
     &  outtime_spinstell    (5000, num_spinstell)
     

         INTEGER
     &  outctr_spinstell    (num_spinstell)


       REAL*8  time1, time2, time, timtot,end_stim,start_stim, tmp1, tmp2
       REAL*8 presyntime, delta, dexparg, dexparg1, dexparg2
       INTEGER thisno, display /0/, O
       REAL*8 z, z1, z2, outrcd(20)

       REAL dtime, t_time(2), e_time ! e_time is dtime output

       double precision::
     & gAMPA_spinstell     (numcomp_spinstell,   num_spinstell),
     & gNMDA_spinstell     (numcomp_spinstell,   num_spinstell),
     & gGABA_A_spinstell     (numcomp_spinstell,   num_spinstell)

         time1 = 0.d0	!	 used to be =gettime()
c Setup tables for calculating exponentials
          call dexptablesmall_setup (dexptablesmall)
           do j = 1, numcomp_spinstell
c       WRITE(6,919) j
c919     FORMAT(I3,X)
       mnaf_spinstell(j,num_spinstell) = 0
       mnap_spinstell(j,num_spinstell) = 0
       hnaf_spinstell(j,num_spinstell) = 0
       mkdr_spinstell(j,num_spinstell) = 0
       mka_spinstell(j,num_spinstell) = 0
       hka_spinstell(j,num_spinstell) =0 
       mk2_spinstell(j,num_spinstell) =0
       hk2_spinstell(j,num_spinstell) =0
       mkm_spinstell(j,num_spinstell)=0
       mkc_spinstell(j,num_spinstell) =0
       mkahp_spinstell(j,num_spinstell)=0
       mcat_spinstell(j,num_spinstell)=0
       hcat_spinstell(j,num_spinstell)=0
       mcal_spinstell(j,num_spinstell)=0
       mar_spinstell(j,num_spinstell)=0
       gAMPA_spinstell     (j,   num_spinstell) =0
       gNMDA_spinstell     (j,   num_spinstell) =0
       gGABA_A_spinstell     (j,   num_spinstell) =0
       curr_spinstell     (j,   num_spinstell) = 0
       V_spinstell    (j,   num_spinstell) =0
       chi_spinstell(numcomp_spinstell,num_spinstell) =0
		end do

           do j = 1, num_spinstell   
        outtime_spinstell(1,j)               = -1.d5
	outctr_spinstell(j)		  = 0
           end do ! j

          call dexptablebig_setup(dexptablebig)


          timtot = 700.0d0 ! debug time does test 0.05 interval 
	ranvec_spinstell(1) =0

c begin define input current
	end_stim = 600.0d0
	start_stim = 200.0d0
	O = 0
1000    O = O + 1
        time = time + dt
        if (time.gt.timtot) goto 2000
	curr_spinstell     (1,   num_spinstell)  = -0.05
 	if (time.lt.1000.d0) then
		if ( time.gt.start_stim.and.time.le.end_stim) then
    	curr_spinstell     (1,   num_spinstell)  = 0.167
c		tmp1 = 1*(-dexp((start_stim-time)*2))
c		tmp2 = 1*(dexp((start_stim-time)/20))
c		tmp1 = (1-exp((start_stim-time)/10))
c		tmp2 = (exp((start_stim-time)/20))
c        curr_spinstell(43,num_spinstell)= 3*tmp1 *  tmp2
c        WRITE(6,927)time,curr_spinstell(43,num_spinstell),tmp1,tmp2,start_stim
c927     FORMAT(E14.7,X,E14.7,X,E14.7,X, E14.7,X,E14.7,X)
		else 
			if (time.gt.end_stim) then
			curr_spinstell(1,num_spinstell) = 0
			endif
		 endif
	 endif

c end define input current

       CALL INTEGRATE_spinstell (O, time, num_spinstell,
     &    V_spinstell, curr_spinstell,
     & gAMPA_spinstell, gNMDA_spinstell , gGABA_A_spinstell ,
     & Mg, 
     & gapcon_spinstell,totaxgj_spinstell,gjtable_spinstell, dt,
     &  chi_spinstell,mnaf_spinstell,mnap_spinstell,
     &  hnaf_spinstell,mkdr_spinstell,mka_spinstell,
     &  hka_spinstell,mk2_spinstell,hk2_spinstell,
     &  mkm_spinstell,mkc_spinstell,mkahp_spinstell,
     &  mcat_spinstell,hcat_spinstell,mcal_spinstell,
     &  mar_spinstell)

        outrcd( 1) = time
        outrcd( 2) = v_spinstell   (1,1)
        outrcd( 3) = v_spinstell   (59   ,1)
        outrcd( 4) = v_spinstell   (43,1)
         z = 0.d0
          do i = 1, num_spinstell   
           z = z - v_spinstell(1,i)
          end do
        outrcd( 5) = z / dble(num_spinstell   ) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_spinstell   
           z = z + gAMPA_spinstell   (i,1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_spinstell   
           z = z + gNMDA_spinstell   (i,1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_spinstell   
           z = z + gGABA_A_spinstell   (i,1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_spinstell   (1,3)
        outrcd(10) = v_spinstell   (1,2)
        outrcd(11) =  curr_spinstell(1,   num_spinstell)
c        outrcd(12) = field_2mm_spinstell

      OPEN(17,FILE='GROUCHO110.spinstell')
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





