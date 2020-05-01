module chemkin
    implicit none
    
    ! unit nuber for input
    integer, parameter :: unit_tplink = 20
    integer, parameter :: unit_cklink = 21
    integer, parameter :: unit_stdout = 6

    ! lenth of work array obtained from cklen
    integer len_logi_cklen
    integer len_int_cklen
    integer len_real_cklen
    integer len_char_cklen
    
    ! length of work array
    integer len_logi_ckwk
    integer len_int_ckwk
    integer len_real_ckwk
    integer len_char_ckwk

    ! chemkin work array
    integer,      allocatable :: int_ckwk(:)
    real(8),      allocatable :: real_ckwk(:)
    character(6), allocatable :: char_ckwk(:)*16
    
    ! transport work array
    integer,      allocatable :: int_tpwk(:)
    real(8),      allocatable :: real_tpwk(:)

    ! index of chemistry
    integer mm ! number of elements
    integer kk ! number of species
    integer ii ! number of reactions
    integer nfit ! number of fitting coefficients
    
    ! array pointers
    integer nsys, neq, lsdas, lidas, lrdas
    integer iprd, iprck, ipwt, ipwdot, ipu, iptot
    integer nrpar, nrdas, nsdas, natol, nrtol, nxmol, nz, nzp
    integer nipar, nidas, ipcck, nksym
    integer, parameter :: ipick = 20

    contains

    subroutine initialize_chemkin_workarray()

        ! open input unit
        open(unit_cklink, form='unformatted', file='./data/cklink')

        !   ------- get chemkin index ---------

        call cklen(unit_cklink, unit_stdout, len_int_cklen, len_real_cklen, len_char_cklen)
        call get_ckindex()

        !   ------- calculate pointer index ---------

        call calc_pointer()

        !   ------- initialize chemkin work array ---------

        allocate(int_ckwk(len_int_ckwk))
        allocate(real_ckwk(len_real_ckwk))
        allocate(char_ckwk(len_char_ckwk))

        call ckinit(len_int_ckwk, len_real_ckwk, len_char_ckwk, unit_cklink, &
                    unit_stdout, int_ckwk(ipick), real_ckwk(iprck), char_ckwk)

    end subroutine initialize_chemkin_workarray

    
    subroutine get_ckindex()

        ! temporary work array
        integer      :: int_ckwk_temp(len_int_cklen)
        real(8)      :: real_ckwk_temp(len_real_cklen)
        character(6) :: char_ckwk_temp(len_char_cklen)*16

        call ckinit(len_int_cklen, len_real_cklen, len_char_cklen, unit_cklink, &
                    unit_stdout, int_ckwk_temp, real_ckwk_temp, char_ckwk_temp)
        call ckindx(int_ckwk_temp, real_ckwk_temp, mm, kk, ii, nfit)

    end subroutine get_ckindex

    subroutine calc_pointer()

        NSYS  = KK + 1
        NEQ = NSYS
        LSDAS  = 1
        LIDAS  = 20 + NEQ
        LRDAS  = 40 + NEQ*9 + NSYS*NSYS
! C  
! C           SET RPAR POINTERS
! C       NOTE: MUST HAVE IPRD 1ST IN RPAR
! C  

        IPRD   = 1
        IPRCK  = IPRD  + II
        IPWT   = IPRCK + len_real_cklen
        IPWDOT = IPWT  + KK
        IPU    = IPWDOT+ KK
        IPTOT  = IPU   + KK
! C  
! C           APPORTION REAL WORK SPACE
! C  
        NRPAR  = 1
        NRDAS  = NRPAR + IPTOT
        NSDAS  = NRDAS + LRDAS
        NATOL  = NSDAS + LSDAS
        NRTOL  = NATOL + NEQ
        NXMOL  = NRTOL + NEQ
        NZ     = NXMOL + KK
        NZP    = NZ    + NEQ
        len_real_ckwk  = NZP   + NEQ
! C  
! C          APPORTION INTEGER WORK SPACE
! C  
        NIPAR  = 1
        NIDAS  = NIPAR + len_int_cklen + IPICK
        len_int_ckwk  = NIDAS + LIDAS
! C  
! C          APPORTION CHARACTER WORK SPACE
! C  
        IPCCK = 1
        NKSYM = IPCCK + len_char_cklen
        len_char_ckwk = NKSYM + KK

    end subroutine calc_pointer

    subroutine get_next_TY(p_cfd, t_cfd, y_cfd, delta_t_cfd, tols_cfd)

        real(8), intent(inout) :: t_cfd
        real(8), intent(in)    :: p_cfd
        real(8), intent(inout) :: y_cfd(kk)
        real(8), intent(in)    :: delta_t_cfd
        real(8), intent(in)    :: tols_cfd(4)

        integer :: icase = 1
        integer  lin, lout, linkck, lsave, lign, lrest
        logical :: lsens = .false.

        DATA LIN/5/, LOUT/6/, LINKCK/25/, LSAVE/7/, LIGN/9/, LREST/10/

        call begin(nsys, neq, icase, ii, kk, len_int_cklen, len_real_cklen,      &
            len_char_cklen, unit_cklink, lin, lout, lsave, lign, lrest, lsens,   &
            lidas, lrdas, lsdas, int_ckwk(nidas), real_ckwk(nrdas),              &
            real_ckwk(nsdas), real_ckwk(nrpar), int_ckwk(nipar), real_ckwk(nz),  &
            real_ckwk(nzp), real_ckwk(nrtol), real_ckwk(natol), real_ckwk(nxmol),&
            char_ckwk(nksym), char_ckwk(ipcck),                                  &
            p_cfd, t_cfd, y_cfd, delta_t_cfd, tols_cfd)

    end subroutine get_next_TY

end module chemkin