program slp
implicit none

! =============== BASIC CONFIGURATION ==============

Integer, parameter :: nx = 30;
Integer, parameter :: ny = 20;
Real*8 :: p  = 3.0;
Real*8 :: volfrac = 0.4;

! Material Properties
Real*8, parameter :: E = 1.0; ! Young's Modulus
Real*8, parameter :: nu = 0.3; ! Poisson Ratio

! ============= END BASIC CONFIGURATION ============

Real*8, dimension(2 * (nx+1)*(ny+1)) :: F, U;
Real*8, dimension(8, 8) :: ke;
Logical, dimension(2, ny + 1, nx + 1) :: supports;
Integer :: support_count, i, j;
Real*8 :: change, obj, best_objective;

Real*8, dimension(ny, nx) :: x, old_x, best_x;
Real*8, dimension(2, ny*nx) :: table;
Real*8, dimension(size(table, 1)) :: rhs;
Integer :: since_last_best;

Character(len=128) :: material_output, compliance_output;
Character(len=10) :: pval, volume;
Character(len=60) :: outfmt = "(a4,i3,a4,i3,a4,f12.6,a10,f12.6,a4,f12.6)";

! =============== LOAD & SUPPORT CONFIGURATION =============

! Loads
F = 0.0;
F(2*(nx + 1) * (ny + 1)) = -1;

! Supports
supports = .false.;
supports(:, 1:ny+1, 1) = .true.;

! =============== END LOAD & SUPPORT CONFIGURATION =========

call getarg(1, pval);
call getarg(2, volume);

read(pval, '(f10.0)') p;
read(volume, '(f10.9)') volfrac;

write(*, *) "Running with p = ", p, " and volfrac = ", volfrac;

support_count = get_support_count();
x = volfrac;

call set_ke(ke);
call solve(U, x);
best_objective = objective(x, U);
best_x = x;
since_last_best = 0;

do i = 1, 40
    old_x = x;
    call solve(U, x);
    call init_table(x, table, rhs, U, ke, volfrac);
    call simplex(table, rhs, x);

    call solve(U, x);
    obj = objective(x, U);

    change = linf(reshape(x - old_x, [ny * nx]));
    write(*, *) "Iteration", i, "Compliance", obj, "Change", change, "Vol", sum(x);

    if (since_last_best > 15) then
        exit;
    end if

    if (obj <= best_objective) then
        best_x = x;
        best_objective = obj;
        since_last_best = 0;
    else
        since_last_best = since_last_best + 1;
    end if
end do

call getarg(3, material_output);
open(unit = 1, file = material_output, status="replace");

do i = 1, ny
    do j = 1, nx-1
        write(1, "(F8.6 A1)", advance="no") best_x(i, j), ",";
    end do
    write(1, "(F8.6)") best_x(i, nx);
end do
close(1)

call getarg(4, compliance_output);
open(unit = 2, file = compliance_output, status="unknown", position="append");
write(2, outfmt) "nx:", nx, "ny:", ny, "p:", p, "volfrac:", volfrac, "z:", best_objective;
close(2);



contains
    ! ================================================================
    ! ============== SEQUENTIAL LINEAR PROGRAMMING ===================
    ! ================================================================
    Real*8 function objective(xvec, U)
        Real*8, dimension(:,:), intent(in) :: xvec;
        Real*8, dimension(:), intent(in) :: U;

        Integer :: ex, ey;
        Integer :: n1, n2;
        Integer, dimension(8) :: edof;
        Real*8, dimension(8) :: ku;
        Real*8, dimension(1) :: uku;
        Real*8 :: total;

        total = 0.0;
        do ex = 1, size(xvec, 2)
            do ey = 1, size(xvec, 1)
                n1 = (ny + 1) * (ex - 1) + ey;
                n2 = (ny + 1) * ex + ey;

                edof = (/ 2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n2+1, 2*n2+2, 2*n1+1, 2*n1+2 /);

                ku = matmul(Ke, U(edof));
                uku = sum(U(edof) * ku);      
                total = total + xvec(ey, ex)**p * uku(1);
            end do
        end do

        objective = total;
    end function objective
        
    Real*8 function partial(xvec, ex, ey, U, Ke)
        Integer, intent(in) :: ex, ey;
        Real*8, dimension(:,:), intent(in) :: xvec;
        Real*8, dimension(:), intent(in) :: U;
        Real*8, dimension(8, 8), intent(in) :: Ke;

        Integer :: n1, n2;
        Integer, dimension(8) :: edof;
        Real*8, dimension(8) :: ku;
        Real*8, dimension(1) :: uku;

        n1 = (ny + 1) * (ex - 1) + ey;
        n2 = (ny + 1) * ex + ey;

        edof = (/ 2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n2+1, 2*n2+2, 2*n1+1, 2*n1+2 /);

        ku = matmul(Ke, U(edof));
        uku = sum(U(edof) * ku);

        partial = p*(xvec(ey,ex))**(p-1) * uku(1);
    end function partial

    subroutine init_table(xvec, table, rhs, U, Ke, volfrac)
        Real*8, dimension(:, :), intent(in) :: xvec;
        Real*8, dimension(:, :), intent(out) :: table;
        Real*8, dimension(:), intent(out) :: rhs;
        Real*8, dimension(:), intent(in) :: U;
        Real*8, dimension(8, 8), intent(in) :: Ke;
        Real*8, intent(in) :: volfrac;

        Integer :: i;
        Integer :: ex, ey;

        do ex = 1, size(xvec, 2)
            do ey = 1, size(xvec, 1)
                i = size(xvec, 1) * (ex - 1) + ey;
                table(1, i) = -partial(xvec, ex, ey, U, Ke);
                table(2, i) = 1.0;
            end do
        end do

        rhs(1) = 0.0;
        rhs(2) = 0.0;

    end subroutine init_table

    ! ===================================================
    ! ================= SIMPLEX METHOD ==================
    ! ===================================================

    subroutine simplex(table, rhs, xvec)
        Real*8, dimension(:, :), intent(inout) :: table;
        Real*8, dimension(:, :), intent(inout) :: xvec;
        Real*8, dimension(:), intent(inout) :: rhs;

        Real*8, dimension(:), allocatable :: temp_x, ub, lb;
        Logical, dimension(size(xvec, 1) * size(xvec, 2)) :: maxxed;
        Integer, dimension(size(table, 1) - 1) :: basic_vars;
        Real*8 :: density;
        Integer :: i;

        Integer :: entering, leaving_idx, leaving;
        Integer :: ex, ey;
        Logical :: opt = .false.;

        allocate(temp_x(size(xvec, 1) * size(xvec, 2)), ub(size(xvec, 1) * size(xvec, 2)));
        allocate(lb(size(xvec, 1) * size(xvec, 2)));
        temp_x = 0.0;
        ub = 1.0;
        maxxed = .false.;

        do ex = 1, size(xvec, 2)
            do ey = 1, size(xvec, 1)
                i = size(xvec, 1) * (ex - 1) + ey;
                density = xvec(ey, ex);
                lb(i) = max(0.001 - density, -0.2);
                ub(i) = min(1.0 - density, 0.2) - lb(i);
            end do
        end do

        call init_simplex(table, temp_x, rhs, lb, ub, maxxed, basic_vars);

        do
            call optimal(table, opt);
            if (opt) then
                exit;
            end if

            call entering_basic_var(table, entering);
            call leaving_basic_var(table, rhs, ub, temp_x, basic_vars, maxxed, leaving, leaving_idx);

            if (entering /= leaving) then
                basic_vars(leaving_idx - 1) = entering;
                call eliminate(table, rhs, entering, leaving_idx);
            end if
        end do

        do ex = 1, size(xvec, 2)
            do ey = 1, size(xvec, 1)
                i = size(xvec, 1) * (ex - 1) + ey;
                density = xvec(ey, ex);
                if (maxxed(i)) then
                    xvec(ey, ex) = density + (ub(i) - temp_x(i) + lb(i));
                else
                    xvec(ey, ex) = density + (temp_x(i) + lb(i));
                end if

                if (xvec(ey, ex) > 1.0 .or. xvec(ey, ex) < 0.001) then
                    write(*, *) "Element went OOB:", ex, ey, xvec(ey, ex)
                end if
            end do
        end do

        if (allocated(temp_x)) then
            deallocate(temp_x);
        end if
        if (allocated(ub)) then
            deallocate(ub)
        end if
        if (allocated(lb)) then
            deallocate(lb)
        end if
    end subroutine simplex

    subroutine optimal(table, cond)
        Real*8, Dimension(:,:), intent(in) :: table;
        Logical, intent(out) :: cond;
        Integer i;

        cond = .true.;
        do i=1, size(table, 2)
            if (table(1,i) < 0) then
                cond = .false.;
            end if
        end do
    end subroutine optimal

    subroutine init_simplex(table, xvec, rhs, lb, ub, maxxed, basic)
        Real*8, dimension(:, :), intent(inout) :: table;
        Real*8, dimension(:), intent(inout) :: xvec;
        Real*8, dimension(:), intent(inout) :: rhs;
        Real*8, dimension(:), intent(in) :: ub, lb;
        Logical, dimension(:), intent(inout) :: maxxed;
        Integer, dimension(:), intent(inout) :: basic;

        Integer :: i, j;

        rhs(2) = rhs(2) - sum(lb);
        do i = 1, size(xvec)
            xvec(i) = 0.0;
            maxxed(i) = .true.;

            do j = 1, size(table, 1)
                rhs(j) = rhs(j) - table(j, i) * ub(i);
                table(j, i) = -table(j, i);
            end do

            if (rhs(2) < ub(i+1)) then
                exit;
            end if
        end do

        i = i + 1;
        basic(1) = i;
        xvec(i) = rhs(2);
        rhs(2) = 0.0;

        call eliminate(table, rhs, i, 2);
    end subroutine init_simplex

    subroutine eliminate(table, rhs, entering, leaving)
        Real*8, Dimension(:, :), intent(inout) :: table;
        Real*8, Dimension(size(table, 1), size(table, 2)) :: output;
        Real*8, Dimension(:), intent(inout) :: rhs;
        Integer, intent(in) :: entering, leaving;
        Integer :: row, col;
        Real*8 :: pivot, in_pivot_c, in_pivot_r;

        output = table;
        pivot = table(leaving, entering);
        do col=1, size(table, 2)
            do row=1, size(table, 1)
                in_pivot_c = table(row, entering);
                in_pivot_r = table(leaving, col);

                if (row /= leaving) then
                    output(row, col) = table(row, col) - in_pivot_r * in_pivot_c / pivot;
                else
                    output(row, col) = table(row, col) / pivot;
                end if
            end do
        end do

        in_pivot_r = rhs(leaving);
        do row=1, size(table, 1)
            in_pivot_c = table(row, entering);
            if (row /= leaving) then
                rhs(row) = rhs(row) - in_pivot_r * in_pivot_c / pivot;
            else
                rhs(row) = rhs(row) / pivot;
            end if
        end do

        table = output;
    end subroutine eliminate

    subroutine entering_basic_var(table, idx)
        Real*8, Dimension(:,:), intent(in) :: table;
        Integer, intent(out) :: idx;
        Integer i;

        idx = 1;
        do i=1, size(table, 2)
            if (table(1, i) < table(1, idx)) then
                idx = i;
            end if
        end do
    end subroutine entering_basic_var

    subroutine leaving_basic_var(table, rhs, ub, xvec, basic, maxxed, leaving, idx)
        Real*8, dimension(:,:), intent(inout) :: table;
        Real*8, dimension(:), intent(inout) :: rhs;
        Real*8, dimension(:), intent(in) :: ub;
        Real*8, dimension(:), intent(inout) :: xvec;
        Integer, dimension(:), intent(inout) :: basic;
        Logical, dimension(:), intent(inout) :: maxxed;
        Integer, intent(out) :: idx;
        Integer, intent(out) :: leaving;
        Integer :: i;
        Integer :: entering;
        Real*8 :: ratio, ratio_idx;
        Logical :: idx_above_bound;

        call entering_basic_var(table, entering);
        leaving = entering;
        ratio_idx = ub(leaving);
        idx_above_bound = .true.;
        do i=2, size(table, 1)
            if (table(i, entering) > 0.0) then
                ratio = rhs(i) / table(i, entering);
                if (ratio < ratio_idx) then
                    idx = i;
                    leaving = basic(idx - 1);
                    ratio_idx = ratio;
                    idx_above_bound = .false.;
                end if
            end if
            if (table(i, entering) < 0.0 .and. basic(i-1) <= size(ub)) then
                ratio = (rhs(i) - ub(basic(i-1))) / table(i, entering);
                if (ratio < ratio_idx) then
                    idx = i;
                    leaving = basic(idx - 1);
                    ratio_idx = ratio;
                    idx_above_bound = .true.;
                end if
            end if
        end do

        if (entering <= size(xvec)) then
            xvec(entering) = ratio_idx;
        end if
        if (idx_above_bound) then
            do i = 1, size(table, 1)
                rhs(i) = rhs(i) - table(i, leaving) * ub(leaving);
                table(i, leaving) = -table(i, leaving);
            end do
            if (leaving <= size(xvec)) then
                xvec(leaving)  = 0.0;
                maxxed(leaving) = .not. maxxed(leaving);
            end if
        end if
    end subroutine leaving_basic_var

    ! ===============================================================
    ! ==================== DISPLACEMENT COMPUTATION =================
    ! ===============================================================

    subroutine solve(U, x)
        Real*8, dimension(:), intent(out) :: U;
        Real*8, dimension(:,:), intent(in) :: x;
        Real*8, dimension(2 * (nx+1)*(ny+1), 2*(nx+1)*(ny+1)) :: Kmat;

        Integer, dimension(size(U) - support_count) :: freedofs;
        Integer, dimension(support_count) :: fixeddofs;

        Real*8, dimension(size(freedofs)) :: Udof, S;
        Real*8, dimension(size(freedofs)) :: Fdof;
        Real*8, dimension(size(freedofs), size(freedofs)) :: Kdof, L;

        call set_kmat(Kmat, x, Ke);
        call get_freedofs(freedofs);
        call get_fixeddofs(fixeddofs);

        U = 0.0;
        Udof = 0.0;
        Fdof = F(freedofs);
        Kdof = Kmat(freedofs, freedofs);

        call chol_ldl(Kdof, L, S);
        call solve_ldl(L, S, Fdof, Udof);

        U(freedofs) = Udof;
        U(fixeddofs) = 0.0;
    end subroutine solve

    subroutine get_freedofs(D)
        Integer, dimension(:), intent(out) :: D;

        Integer supp, row, col, c;
        c = 1;

        do col = 1, size(supports, 3)
            do row = 1, size(supports, 2)
                do supp = 1, size(supports, 1)
                    if (.not. supports(supp, row, col)) then
                        D(c) = 2*size(supports, 2) * (col - 1) + 2*(row - 1) + supp;
                        c = c + 1;
                    end if
                end do
            end do
        end do
    end subroutine get_freedofs

    subroutine get_fixeddofs(D)
        Integer, dimension(:), intent(out) :: D;

        Integer supp, row, col, c;
        c = 1;

        do col = 1, size(supports, 3)
            do row = 1, size(supports, 2)
                do supp = 1, size(supports, 1)
                    if (supports(supp, row, col)) then
                        D(c) = 2*size(supports, 2) * (col - 1) + 2*(row - 1) + supp;
                        c = c + 1;
                    end if
                end do
            end do
        end do
    end subroutine get_fixeddofs

    Integer function get_support_count()
        Integer :: i, j, k, c;

        c = 0;
        do k=1, size(supports, 3)
            do j=1, size(supports, 2)
                do i=1, size(supports, 1)
                    if (supports(i, j, k)) then
                        c = c + 1
                    end if
                end do
            end do
        end do

        get_support_count = c;
    end function get_support_count

    subroutine set_kmat(kmat, x, ke)
        Real*8, dimension(:, :), intent(inout) :: kmat;
        Real*8, dimension(:, :), intent(in) :: x;
        Real*8, dimension(:, :), intent(in) :: ke;
        Integer, dimension(8) :: edof;
        Integer :: i, j;
        Integer :: n1, n2;

        kmat = 0.0;
        do j = 1, nx
            do i = 1, ny
                n1 = (ny+1) * (j-1) + i;
                n2 = (ny+1) * j + i;

                edof = (/ 2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n2+1, 2*n2+2, 2*n1+1, 2*n1+2 /);
                kmat(edof, edof) = kmat(edof, edof) + x(i, j)**p * ke;
            end do
        end do

    end subroutine set_kmat

    subroutine set_ke(ke)
        Real*8, dimension(8, 8), intent(inout) :: ke;
        Real*8, dimension(8), parameter :: k = (/ &
            0.5 - nu/6, &
            0.125 + nu/8, &
            -0.25 - nu/12, &
            -0.125 + 3 * nu / 8, &
            -0.25 + nu / 12, &
            -0.125 - nu /8, &
            nu / 6, &
            0.125 - 3 * nu / 8 &
        /);

        ke(1, :) = (/ k(1), k(2), k(3), k(4), k(5), k(6), k(7), k(8) /);
        ke(2, :) = (/ k(2), k(1), k(8), k(7), k(6), k(5), k(4), k(3) /);
        ke(3, :) = (/ k(3), k(8), k(1), k(6), k(7), k(4), k(5), k(2) /);
        ke(4, :) = (/ k(4), k(7), k(6), k(1), k(8), k(3), k(2), k(5) /);
        ke(5, :) = (/ k(5), k(6), k(7), k(8), k(1), k(2), k(3), k(4) /);
        ke(6, :) = (/ k(6), k(5), k(4), k(3), k(2), k(1), k(8), k(7) /);
        ke(7, :) = (/ k(7), k(4), k(5), k(2), k(3), k(8), k(1), k(6) /);
        ke(8, :) = (/ k(8), k(3), k(2), k(5), k(4), k(7), k(6), k(1) /);

        ke = E/(1 - nu**2) * ke;
    end subroutine set_ke

    subroutine chol_ldl(A, L, S)
        Real*8, dimension(:, :), intent(in) :: A;
        Real*8, dimension(:, :), intent(out) :: L;
        Real*8, dimension(:), intent(out) :: S;

        Real*8, dimension(size(A, 1), size(A, 2)) :: B;
        Integer i, j;

        B = A;
        L = 0.0;
        do i = 1, size(L, 1)
            L(i, i) = 1.0;
        end do

        do j = 1, size(L, 2)
            S(j) = B(j, j);
            L(j+1:, j) = B(j+1:, j) / B(j, j);
            do i = j + 1, size(L, 2)
                B(j+1:, i) = B(j+1:, i) - B(j, i) / B(j, j) * B(j+1:, j);
            end do
        end do
    end subroutine chol_ldl

    subroutine solve_ldl(L, S, b, x)
        Real*8, dimension(:, :), intent(in) :: L;
        Real*8, dimension(:), intent(in) :: S;
        Real*8, dimension(:), intent(in) :: b;
        Real*8, dimension(:), intent(inout) :: x;

        Real*8, dimension(size(x)) :: y;
        Integer :: i;

        y = 0.0;
        y(1) = b(1);
        do i = 2, size(y)
            y(i) = b(i) - dot_product(L(i, 1:i-1), y(1:i-1));
        end do

        y = y / S;

        x = 0.0;
        x(size(x)) = y(size(x));
        do i = size(x) - 1, 1, -1
            x(i) = y(i) - dot_product(L(i+1:size(x), i), x(i+1:size(x)));
        end do

    end subroutine solve_ldl
    
    Real*8 function linf(x)
        Real*8, dimension(:) :: x;
        Integer i;
        linf = abs(x(1));
        do i = 1, size(x)
            linf = max(linf, abs(x(i)))
        end do
    end function linf
end program slp
