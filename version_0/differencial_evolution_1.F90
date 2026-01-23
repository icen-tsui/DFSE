!     
! File:   differencial_evolution.F90
! Author: zyy
!
! Created on 2013年4月23日, 下午9:47
!

module differencial_evolution
    contains
    logical function is_magnetic_atom(ptype)
        use kinds
        use parameters, only : altermagnetic, magnetic_types, num_magnetic_types
        implicit none
        integer(i4b), intent(in) :: ptype
        integer(i4b) :: idx

        is_magnetic_atom = .false.
        if(.not. altermagnetic) return

        do idx = 1, num_magnetic_types
            if(ptype == magnetic_types(idx)) then
                is_magnetic_atom = .true.
                exit
            end if
        end do
    end function is_magnetic_atom

    subroutine run_de_lammps(step)
        use kinds
        use lattice_mod, only : init_lat_nosym
        use init_struct, only : fill_atom
        use parameters, only : population, pstruct, pool, atom_dis, de_ratio
        use sort_atom_idx, only : kmatch
        use init_struct
        implicit none
        integer(i4b), intent(in) :: step
        integer(i4b) :: i,j,k
        real(dp) :: tmp, cr = 0.9, f = 0.5
        integer(i4b) :: a,b,c, a1, a2, n
        real(dp) :: p1(3), p2(3), p3(3), y(3), lat(3,3), dis
        logical :: f1
        real(dp) :: center(3)
        do i = 1, population
            if(i > floor(population * de_ratio)) then
                call init_lat_nosym(i)
                call fill_atom(i)
                cycle
            end if
            do
                call random_number(tmp)
                a = floor(tmp * population * de_ratio + 1.0)
                if(a /= i) exit
            end do
            do
                call random_number(tmp)
                b = floor(tmp * population * de_ratio + 1.0)
                if(b /= i .and. b /= a) exit
            end do
            do
                call random_number(tmp)
                c = floor(tmp * population * de_ratio + 1.0)
                if(c /= i .and. c /= a .and. c /= b) exit
            end do
            call kmatch(pstruct(i), pool(a))
            call kmatch(pstruct(i), pool(b))
            call kmatch(pstruct(i), pool(c))
            n = pstruct(i) % natom
            lat = pstruct(i) % lat
            do j = 1, n
                if(is_magnetic_atom(pstruct(i) % ptype(j))) cycle
                p1(:) = pool(a) % pos(:, j)
                p2(:) = pool(b) % pos(:, j)
                p3(:) = pool(c) % pos(:, j)
                y(:) = p1(:) + f * (p2(:) - p3(:))
                do k = 1, 3
                    y(k) = y(k) - lat(k, k) * real(floor(y(k) / lat(k, k)))
                end do
                f1 = .true.
                do k = 1, n
                    if(k == j) cycle
                    call find_dis(lat, pstruct(i) % pos(:, k), y(:), dis)
                    a1 = pstruct(i) % ptype(k)
                    a2 = pstruct(i) % ptype(j)
                    if(dis < atom_dis(a1, a2)) then
                        f1 = .false.
                    end if
                end do
                if(f1) then
                    pstruct(i) % pos(:, j) = y(:)
                end if
            end do
        end do
    end subroutine run_de_lammps
    
    subroutine run_de_vasp(step)
        use kinds
        use constants
        use parameters, only : pstruct, pool, population, de_ratio
        use parameters, only : Q2D, vacuum_layer, cluster
        use parameters, only : cluster_substrate, num_species, num_ele, SubstrateElements
        use parameters, only : fix_atom
        use parameters, only : polarization_2d,polarization_3d
        use sort, only : sort_results
        use init_struct_all, only : init_struct_tot
        use sort_atom_idx, only : kmatch
        implicit none
        integer(i4b), intent(in) :: step
        integer(i4b) :: i,j,k, ratio_cut, a,b, i1
        integer :: nat, jpair
        real(dp) :: tmp, y(3), p1(3), p2(3), p3(3), pbest(3)
!        real(dp) :: lambd = 0.2, F = 0.2
        real(dp) :: lambd = 0.1, F = 0.1
        real(dp) :: c1, c2, c3, ratio, max_z, min_z
        logical :: flag
        logical :: is_substrate(max_atom * max_type)

        write(*, *) "start de_operation"
        do i = 1, population
            write(*, *) "energy now : old", pstruct(i) % energy, pool(i) % energy
            if (pstruct(i) % energy < pool(i) % energy) then
                write(*, *) "renew: ", i
                pool(i) = pstruct(i)
            end if
        end do

        call sort_results()
        ratio_cut = floor(de_ratio * population)
        
        do i = 1, population
            write(*, *) "de_struct: ", i
            write(*, *) "atom_type: ", (pstruct(i) % ptype(k), k = 1, pstruct(i) % natom)
            if(i > ratio_cut) then
                call init_struct_tot(i)
                cycle
            end if
            do
                call random_number(tmp)
                a = floor(tmp * ratio_cut + 1)
                if(a /= i) exit
            end do
            do
                call random_number(tmp)
                b = floor(tmp * ratio_cut + 1)
                if(b /= i .and. b /= a) exit
            end do
            !call kmatch(pool(1), pstruct(i))
            !call kmatch(pool(1), pool(a))
            !call kmatch(pool(1), pool(b))
            !atoms of the substrate do not perform the de operation
            if(cluster_substrate) then
                i1 = 0
                do j = 1, num_species
                    do k = 1, num_ele(j)
                        i1 = i1 + 1
                        if(k <= SubstrateElements(j)) then
                            is_substrate(i1) = .true.
                        else
                            is_substrate(i1) = .false.
                        end if
                    end do
                end do
            end if
            do j = 1, pstruct(i) % natom
                if(is_magnetic_atom(pstruct(i) % ptype(j))) cycle
                pbest(:) = pool(1) % pos(:, j)
                p1(:) = pool(i) % pos(:, j)
                p2(:) = pool(a) % pos(:, j)
                p3(:) = pool(b) % pos(:, j)
                y(:) = lambd * pbest(:) + (1-lambd) * p1(:) + F * (p2(:) - p3(:))
                y(:) = y(:) - floor(y(:))
                if(cluster_substrate) then
                    if(.not. is_substrate(j)) then
                        pstruct(i) % pos(:, j) = y(:)
                    end if
                else
                    pstruct(i) % pos(:, j) = y(:)
                end if
            end do
 
            call check(i, flag)
            if (flag .eqv. .false.) then
                call init_struct_tot(i)    
                write(*,*) "===generate new structure==="            
            end if
            ! for quasi-2D case
            ! to keep the vacuum_layer after DE operations
            if(Q2D) then
                min_z = 1.0
                max_z = 0.0
                do j = 1, pstruct(i) % natom
                    if(min_z > pstruct(i) % pos(3, j)) min_z = pstruct(i) % pos(3, j)
                    if(max_z < pstruct(i) % pos(3, j)) max_z = pstruct(i) % pos(3, j)
                end do
                c1 = (1.0 - (max_z - min_z)) * pstruct(i) % lat(3, 3)
                c2 = (max_z - min_z) * pstruct(i) % lat(3, 3) + vacuum_layer
                ratio = c1 / c2
                do j = 1, pstruct(i) % natom
                    pstruct(i) % pos(3, j) = pstruct(i) % pos(3, j) * ratio
                end do
                pstruct(i) % lat(3, 3) = c2
            end if
        end do
    end subroutine run_de_vasp

    subroutine run_de_vasp_symmetric(step)
    use kinds
    use constants
    use parameters, only : pstruct, pool, population, de_ratio
    use parameters, only : Q2D, vacuum_layer, cluster
    use parameters, only : cluster_substrate, num_species, num_ele, SubstrateElements
    use parameters, only : fix_atom
    use parameters, only : polarization_2d, polarization_3d
    use parameters, only : altermagnetic
    use sort, only : sort_results
    use init_struct_all, only : init_struct_tot
    use sort_atom_idx, only : kmatch
    implicit none

    integer(i4b), intent(in) :: step
    integer(i4b) :: i, j, k, ratio_cut, a, b, i1
    integer :: nat, jpair, jj
    real(dp) :: tmp, y(3), p1(3), p2(3), p3(3), pbest(3)
    real(dp) :: lambd, F
    real(dp) :: c1, c2, ratio, max_z, min_z
    real(dp) :: center(3)
    real(dp) :: rj(3), rmirror(3)
    logical :: flag
    logical :: is_substrate(max_atom * max_type)
    logical :: is_self
    real(dp), parameter :: tol = 1.0e-6_dp

    lambd = 0.1_dp
    F     = 0.1_dp

    write(*,*) ">>> start symmetric DE operation <<<"

    ! update best pool
    do i = 1, population
        if (pstruct(i)%energy < pool(i)%energy) then
            pool(i) = pstruct(i)
        end if
    end do

    call sort_results()
    ratio_cut = floor(de_ratio * population)

    do i = 1, population
        write(*,*) "de_struct(symmetric): ", i
        write(*,*) "atom_type: ", (pstruct(i)%ptype(k), k = 1, pstruct(i)%natom)

        if (i > ratio_cut) then
            call init_struct_tot(i)
            cycle
        end if

        ! random pick a, b
        do
            call random_number(tmp)
            a = floor(tmp * ratio_cut + 1)
            if (a /= i) exit
        end do
        do
            call random_number(tmp)
            b = floor(tmp * ratio_cut + 1)
            if (b /= i .and. b /= a) exit
        end do

        ! mark substrate atoms if any
        if (cluster_substrate) then
            i1 = 0
            do j = 1, num_species
                do k = 1, num_ele(j)
                    i1 = i1 + 1
                    if (k <= SubstrateElements(j)) then
                        is_substrate(i1) = .true.
                    else
                        is_substrate(i1) = .false.
                    end if
                end do
            end do
        end if

        nat = pstruct(i)%natom

        ! fixed inversion center for ISYM=2
        center = (/ 0.5_dp, 0.5_dp, 0.5_dp /)

        ! symmetric pair-wise DE mutation
        if (polarization_2d .or. polarization_3d) then

            do j = 1, nat/2
                jpair = nat - j + 1

                if(altermagnetic) then
                    if(is_magnetic_atom(pstruct(i)%ptype(j)) .or. is_magnetic_atom(pstruct(i)%ptype(jpair))) cycle
                end if

                ! --- substrate logic ---
                if (cluster_substrate) then
                    if (is_substrate(j) .and. is_substrate(jpair)) cycle
                    if (is_substrate(j) .and. .not. is_substrate(jpair)) then
                        do k = 1, 3
                            pstruct(i)%pos(k, jpair) = 2.0_dp*center(k) - pstruct(i)%pos(k, j)
                            pstruct(i)%pos(k, jpair) = pstruct(i)%pos(k, jpair) - floor(pstruct(i)%pos(k, jpair))
                        end do
                        cycle
                    else if (.not. is_substrate(j) .and. is_substrate(jpair)) then
                        do k = 1, 3
                            pstruct(i)%pos(k, j) = 2.0_dp*center(k) - pstruct(i)%pos(k, jpair)
                            pstruct(i)%pos(k, j) = pstruct(i)%pos(k, j) - floor(pstruct(i)%pos(k, j))
                        end do
                        cycle
                    end if
                end if

                ! --- self-mirror check ---
                rj(:) = pstruct(i)%pos(:, j)
                do k = 1, 3
                    rmirror(k) = 2.0_dp*center(k) - rj(k)
                    rmirror(k) = rmirror(k) - floor(rmirror(k))
                end do
                is_self = all(abs(mod(2.0_dp*center(:) - rj(:), 1.0_dp) - rj(:)) < tol)
                if (is_self) then
                    ! keep atom fixed at its symmetric site
                    do k = 1, 3
                        pstruct(i)%pos(k, j) = rj(k)
                        pstruct(i)%pos(k, jpair) = rj(k)
                    end do
                    cycle
                end if

                ! --- perform DE mutation for non-self pairs ---
                p1(:) = pool(i)%pos(:, j)
                p2(:) = pool(a)%pos(:, j)
                p3(:) = pool(b)%pos(:, j)
                pbest(:) = pool(1)%pos(:, j)

                y(:) = lambd * pbest(:) + (1.0_dp - lambd) * p1(:) + F * (p2(:) - p3(:))

                ! wrap to [0,1)
                do k = 1, 3
                    y(k) = y(k) - floor(y(k))
                end do

                if (polarization_2d) y(3) = p1(3)

                pstruct(i)%pos(:, j) = y(:)

                do k = 1, 3
                    pstruct(i)%pos(k, jpair) = 2.0_dp*center(k) - y(k)
                    pstruct(i)%pos(k, jpair) = pstruct(i)%pos(k, jpair) - floor(pstruct(i)%pos(k, jpair))
                end do

            end do  ! end pair loop

            if (mod(nat, 2) == 1) then
                jj = (nat+1)/2
                if(.not. (altermagnetic .and. is_magnetic_atom(pstruct(i)%ptype(jj)))) then
                    pstruct(i)%pos(:, jj) = center(:)
                end if
            end if

        else
            ! ordinary DE mutation (no symmetry enforced)
            do j = 1, nat
                if(is_magnetic_atom(pstruct(i)%ptype(j))) cycle
                pbest(:) = pool(1)%pos(:, j)
                p1(:) = pool(i)%pos(:, j)
                p2(:) = pool(a)%pos(:, j)
                p3(:) = pool(b)%pos(:, j)
                y(:) = lambd * pbest(:) + (1.0_dp - lambd) * p1(:) + F * (p2(:) - p3(:))
                do k = 1, 3
                    y(k) = y(k) - floor(y(k))
                end do
                pstruct(i)%pos(:, j) = y(:)
            end do
        end if

        call check(i, flag)
        if (.not. flag) then
            call init_struct_tot(i)
            write(*,*) "=== regenerate invalid structure ==="
        end if

        if (Q2D) then
            min_z = 1.0_dp
            max_z = 0.0_dp
            do j = 1, nat
                if (min_z > pstruct(i)%pos(3, j)) min_z = pstruct(i)%pos(3, j)
                if (max_z < pstruct(i)%pos(3, j)) max_z = pstruct(i)%pos(3, j)
            end do
            c1 = (1.0_dp - (max_z - min_z)) * pstruct(i)%lat(3, 3)
            c2 = (max_z - min_z) * pstruct(i)%lat(3, 3) + vacuum_layer
            ratio = c1 / c2
            do j = 1, nat
                pstruct(i)%pos(3, j) = pstruct(i)%pos(3, j) * ratio
            end do
            pstruct(i)%lat(3,3) = c2
        end if

    end do  ! population loop

    write(*,*) ">>> symmetric DE operation done <<<"
end subroutine run_de_vasp_symmetric


    subroutine check(tag, flag)
        use kinds
        use parameters, only : pstruct, atom_dis
        use parameters, only : model_gb, bottom_height, gb_height, top_height
        use init_struct, only : dir2car, find_dis
        implicit none
        integer(i4b), intent(in) :: tag
        logical, intent(out) :: flag
        integer(i4b) :: i,j,k,a1,a2,n
        real(dp) :: dir0(3), dir1(3), dir2(3), car1(3), car2(3), lattice(3,3)
        real(dp) :: distance
        integer(i4b) :: dpos(3, 27), dnum
        n = pstruct(tag) % natom
        lattice = pstruct(tag) % lat
        
        dnum = 0
        do i = -1, 1
            do j = -1, 1
                do k = -1, 1
                    dnum = dnum + 1
                    dpos(1,dnum) = i
                    dpos(2,dnum) = j
                    dpos(3,dnum) = k
                end do
            end do
        end do
        
        
        flag = .true.
        do i = 1, n
            do j = 1, i - 1
                dir0 = pstruct(tag) % pos(:, i)
                dir2 = pstruct(tag) % pos(:, j)
                do k = 1, 27
                    dir1(:) = dpos(:, k) + dir0(:)
                    call dir2car(lattice, dir1, car1)
                    call dir2car(lattice, dir2, car2)
                    call find_dis(lattice, car1, car2, distance)
                    a1 = pstruct(tag) % ptype(i)
                    a2 = pstruct(tag) % ptype(j)
                    if(distance < atom_dis(a1, a2) * 0.4) then
                        flag = .false.
                        return
                    end if
                end do
            end do
        end do

        ! for grain boundary case
        ! if atom is not in the gb region, return false
        if(model_gb) then
            do i = 1, n
                dir0 = pstruct(tag) % pos(:, i)
                if(dir0(3) < (bottom_height - 0.2 * gb_height) / lattice(3, 3) .or. &
                    & dir0(3) > (1 - (top_height - 0.2 * gb_height)/ lattice(3, 3))) then
                    flag = .false.
                end if
            end do
        end if

    end subroutine check
    
    subroutine run_mode_vasp(step)
        use kinds
        use parameters, only : pstruct, pool, population, de_ratio
        use sort, only : sort_results_mode
        use init_struct_spg, only : init_struct_sym
        implicit none
        integer(i4b), intent(in) :: step
        integer(i4b) :: i,j,k,ratio_cut, a, b, c, nf
        real(dp) :: tmp, y(3), p1(3), p2(3), p3(3), pbest(3)
        real(dp) :: lambd = 0.1, F = 0.2
        logical :: flag
        do i = 1, population
            if(pstruct(i) % energy < pool(i) % energy .and. pstruct(i) % hardness < pool(i) % hardness) then
                pool(i) = pstruct(i)
            end if
            call random_number(tmp)
            if(pstruct(i) % energy < pool(i) % energy .or. pstruct(i) % hardness < pool(i) % hardness) then
                if(tmp < 0.5) then
                    pool(i) = pstruct(i)
                end if
            end if
        end do
        write(*, *) "start sort"
        call sort_results_mode()
        write(*, *) "end sort"
        ratio_cut = floor(de_ratio * population)
        nf = 0
        do i = 1, population
            if(pool(i) % frontier == 1) then
                nf = nf + 1
            end if
        end do
        write(*, *) "number of struct in 1st front: ", nf
        do i = 1, population
            write(*, *) "de on: ", i
            if(pool(i) % energy > 100 .or. pool(i) % hardness > 100) then
                write(*, *) "de new struct1: ", i
                call init_struct_sym(i)
                cycle
            end if
            if(i > ratio_cut) then
                write(*, *) "de new struct2: ", i
                call init_struct_sym(i)
                cycle
            end if
            write(*, *) "de operation: ", i
            do
                call random_number(tmp)
                a = floor(tmp * ratio_cut + 1)
                if(a /= i) exit
            end do
            do
                call random_number(tmp)
                b = floor(tmp * ratio_cut + 1)
                if(b /= i .and. b /= a) exit
            end do
            if(pool(i) % frontier /= 1) then
                call random_number(tmp)
                c = floor(tmp * nf + 1)
            end if
            if(pool(i) % frontier == 1) then
                do j = 1, pstruct(i) % natom
                    if(is_magnetic_atom(pstruct(i) % ptype(j))) cycle
                    p1(:) = pool(i) % pos(:, j)
                    p2(:) = pool(a) % pos(:, j)
                    p3(:) = pool(b) % pos(:, j)
                    y(:) = p1(:) + F * (p2(:) - p3(:))
                    y(:) = y(:) - floor(y(:))
                    pstruct(i) % pos(:, j) = y(:)
                end do
            else
                do j = 1, pstruct(i) % natom
                    if(is_magnetic_atom(pstruct(i) % ptype(j))) cycle
                    pbest(:) = pool(c) % pos(:, j)
                    p1(:) = pool(i) % pos(:, j)
                    p2(:) = pool(a) % pos(:, j)
                    p3(:) = pool(b) % pos(:, j)
                    y(:) = lambd * pbest(:) + (1 - lambd) * p1(:) + F * (p2(:) - p3(:))
                    y(:) = y(:) - floor(y(:))
                    pstruct(i) % pos(:, j) = y(:)
                end do
            end if
            call check(i,flag)
            if(flag .eqv. .false.) then
                call init_struct_sym(i)
            end if
        end do
    end subroutine run_mode_vasp
end module differencial_evolution
