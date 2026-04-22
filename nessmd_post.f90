program nessmd_periodic_flux_analysis
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: max_types = 100
  integer, parameter :: max_atoms_per_mol = 5000
  integer, parameter :: max_line = 1024
  real(dp), parameter :: hugev = 1.0d300

  ! Files
  character(len=256) :: field_file, config_file, history_file
  integer :: uf, uc, uh, ugeo, uout, usum, ios

  ! FIELD data
  integer :: nmoltypes, t, i, j, k
  character(len=64) :: molname(max_types)
  integer :: nmols(max_types), natoms_type(max_types)
  character(len=32), allocatable :: type_atom_labels(:,:)
  logical :: is_gas_type(max_types), is_carbon_type(max_types), is_wall_type(max_types)
  integer :: total_atoms_expected, total_molecules_expected
  integer :: atom_start_type(max_types), atom_end_type(max_types)
  integer :: mol_start_type(max_types), mol_end_type(max_types)

  ! CONFIG/HISTORY headers
  character(len=max_line) :: line, title
  integer :: levcfg_cfg, imcon_cfg, natms_cfg
  integer :: levcfg_hist, imcon_hist, natms_hist
  logical :: has_vel_cfg, has_force_cfg, has_vel_hist, has_force_hist
  real(dp) :: cell(3,3), lx, ly, lz, area

  ! Atom/molecule maps
  integer, allocatable :: atom_type(:), atom_mol_global(:), atom_local_in_mol(:)
  integer, allocatable :: mol_type_of_global(:), mol_first_atom(:), mol_natoms(:)
  logical, allocatable :: is_carbon_atom(:)

  ! Trajectory storage per frame
  real(dp), allocatable :: zatom(:), zcom(:), zcom_prev(:)
  logical, allocatable :: mol_has_prev(:)

  ! Gas species bookkeeping
  integer :: ngas_types, gas_type_list(max_types)
  integer :: out_cols_per_gas

  ! Geometry per frame
  real(dp) :: zc, zmin_rel, zmax_rel, slab_thickness
  real(dp) :: left_len, right_len

  ! Per-frame counts/densities
  integer :: stepnum, frame_count
  real(dp) :: dt_ps, time_ps
  integer :: nleft(max_types), nslab(max_types), nright(max_types)
  real(dp) :: rho_left(max_types), rho_right(max_types), delta_rho(max_types)
  integer :: nleft_tot, nslab_tot, nright_tot
  real(dp) :: rho_left_tot, rho_right_tot, delta_rho_tot

  ! Flux/deff accumulators per gas type and total
  integer :: nfwd(max_types), nbwd(max_types)
  integer :: nfwd_tot, nbwd_tot
  real(dp) :: flux(max_types), deff_signed(max_types), deff_abs(max_types)
  real(dp) :: flux_tot, deff_signed_tot, deff_abs_tot
  real(dp) :: plane_rel, plane_abs

  ! Optional analysis control
  logical :: have_input
  character(len=256) :: input_file
  integer :: uin
  real(dp) :: analysis_start_ps

  field_file = 'FIELD'
  config_file = 'CONFIG'
  history_file = 'HISTORY'
  input_file = 'nessmd_analysis.in'

  analysis_start_ps = 0.0_dp
  inquire(file=trim(input_file), exist=have_input)
  if (have_input) then
     open(newunit=uin, file=trim(input_file), status='old', action='read', iostat=ios)
     if (ios == 0) then
        read(uin,*,iostat=ios) dt_ps
        if (ios /= 0) stop 'ERROR reading timestep from nessmd_analysis.in'
        read(uin,*,iostat=ios) analysis_start_ps
        if (ios /= 0) analysis_start_ps = 0.0_dp
        close(uin)
     else
        stop 'ERROR opening nessmd_analysis.in'
     end if
  else
     dt_ps = -1.0_dp
  end if

  call read_field()
  call build_maps()
  call read_config_header_and_cell()
  call open_history_and_header()

  if (dt_ps <= 0.0_dp) then
     write(*,*) 'No nessmd_analysis.in found. Using timestep from HISTORY if present.'
  end if

  allocate(zatom(total_atoms_expected))
  allocate(zcom(total_molecules_expected))
  allocate(zcom_prev(total_molecules_expected))
  allocate(mol_has_prev(total_molecules_expected))
  mol_has_prev = .false.
  zcom_prev = 0.0_dp

  open(newunit=ugeo, file='nessmd_geometry.txt', status='replace', action='write')
  write(ugeo,'(A)') 'Detected molecule types from FIELD:'
  do t = 1, nmoltypes
     write(ugeo,'(I4,2X,A,2X,A,I8,2X,A,I8,2X,A,L1,2X,A,L1,2X,A,L1)') t, trim(molname(t)), &
          'NUMMOLS=', nmols(t), 'ATOMS=', natoms_type(t), 'gas=', is_gas_type(t), &
          'carbon=', is_carbon_type(t), 'wall=', is_wall_type(t)
  end do
  write(ugeo,'(A,F12.4)') 'Box Lx (A): ', lx
  write(ugeo,'(A,F12.4)') 'Box Ly (A): ', ly
  write(ugeo,'(A,F12.4)') 'Box Lz (A): ', lz
  write(ugeo,'(A,F12.4)') 'Cross-sectional area A (A^2): ', area
  write(ugeo,'(A,F12.6)') 'Analysis start time (ps): ', analysis_start_ps
  close(ugeo)

  open(newunit=uout, file='nessmd_regions_flux.dat', status='replace', action='write')
  call write_output_header(uout)

  nfwd = 0; nbwd = 0; nfwd_tot = 0; nbwd_tot = 0; frame_count = 0

  do
     call read_next_frame(stepnum, time_ps, ios, plane_abs)
     if (ios /= 0) exit
     frame_count = frame_count + 1

     call compute_carbon_geometry(zc, zmin_rel, zmax_rel, slab_thickness, left_len, right_len)
     plane_rel = 0.5_dp*(zmin_rel + zmax_rel)
     plane_abs = wrap01(zc + plane_rel, lz)

     call compute_all_mol_com(zc)
     call classify_and_flux(time_ps, zc, zmin_rel, zmax_rel, left_len, right_len, plane_rel)
     call write_frame(uout, stepnum, time_ps, zc, slab_thickness, left_len, right_len, plane_abs)
  end do

  close(uout)
  close(uh)

  open(newunit=usum, file='nessmd_summary.txt', status='replace', action='write')
  write(usum,'(A,I8)') 'Frames processed: ', frame_count
  write(usum,'(A,F12.6)') 'Timestep used (ps): ', dt_ps
  write(usum,'(A,F12.6)') 'Analysis start time (ps): ', analysis_start_ps
  write(usum,'(A)') ''
  do k = 1, ngas_types
     t = gas_type_list(k)
     write(usum,'(A)') 'Species: '//trim(molname(t))
     write(usum,'(A,I10)') '  Nfwd: ', nfwd(t)
     write(usum,'(A,I10)') '  Nbwd: ', nbwd(t)
     write(usum,'(A,ES16.8)') '  Final flux (molecules/A^2/ps): ', flux(t)
     write(usum,'(A,ES16.8)') '  Final D_eff_signed (A^2/ps): ', deff_signed(t)
     write(usum,'(A,ES16.8)') '  Final D_eff_abs (A^2/ps): ', deff_abs(t)
     write(usum,'(A,ES16.8)') '  Final D_eff_abs (m^2/s): ', deff_abs(t)*1.0d-8
     write(usum,'(A)') ''
  end do
  write(usum,'(A)') 'Total gas mixture'
  write(usum,'(A,I10)') '  Nfwd: ', nfwd_tot
  write(usum,'(A,I10)') '  Nbwd: ', nbwd_tot
  write(usum,'(A,ES16.8)') '  Final flux (molecules/A^2/ps): ', flux_tot
  write(usum,'(A,ES16.8)') '  Final D_eff_signed (A^2/ps): ', deff_signed_tot
  write(usum,'(A,ES16.8)') '  Final D_eff_abs (A^2/ps): ', deff_abs_tot
  write(usum,'(A,ES16.8)') '  Final D_eff_abs (m^2/s): ', deff_abs_tot*1.0d-8
  close(usum)

  write(*,*) 'Done.'
  write(*,*) 'Outputs:'
  write(*,*) '  nessmd_regions_flux.dat'
  write(*,*) '  nessmd_geometry.txt'
  write(*,*) '  nessmd_summary.txt'

contains

  subroutine read_field()
    integer :: atoms_read
    character(len=max_line) :: l
    character(len=64) :: key
    character(len=32) :: field_atom_label
    integer :: n, pos

    allocate(type_atom_labels(max_types, max_atoms_per_mol))
    type_atom_labels = ''
    molname = ''
    nmols = 0
    natoms_type = 0
    is_gas_type = .false.
    is_carbon_type = .false.
    is_wall_type = .false.

    open(newunit=uf, file=trim(field_file), status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'ERROR opening FIELD'

    call next_data_line(uf, l, ios)
    if (ios /= 0) stop 'ERROR reading FIELD title'

    nmoltypes = -1
    do
      call next_data_line(uf, l, ios)
      if (ios /= 0) exit
      if (starts_with_word(l, 'MOLECULES')) then
         read(l,*) key, nmoltypes
         exit
      end if
    end do
    if (nmoltypes <= 0) stop 'ERROR parsing FIELD MOLECULES line'
    if (nmoltypes > max_types) stop 'Increase max_types'

    do t = 1, nmoltypes
      call next_data_line(uf, l, ios)
      if (ios /= 0) stop 'ERROR reading molecule name in FIELD'
      molname(t) = adjustl(trim(l))

      call next_data_line(uf, l, ios)
      if (ios /= 0) stop 'ERROR reading NUMMOLS line in FIELD'
      read(l,*) key, nmols(t)

      call next_data_line(uf, l, ios)
      if (ios /= 0) stop 'ERROR reading ATOMS line in FIELD'
      read(l,*) key, natoms_type(t)
      if (natoms_type(t) > max_atoms_per_mol) stop 'Increase max_atoms_per_mol'

      ! Read atom definitions. DL_POLY FIELD may use one line per atom site,
      ! or a compressed form with a repeat count, e.g.:
      !   C   12.0110   0.0000   624   1
      ! meaning 624 repeated sites of label C.
      j = 0
      do while (j < natoms_type(t))
        call next_data_line(uf, l, ios)
        if (ios /= 0) stop 'ERROR reading atom definition in FIELD'
        call parse_field_atom_line(l, field_atom_label, k)
        if (k <= 0) k = 1
        if (j + k > natoms_type(t)) then
           stop 'ERROR: repeat count in FIELD atom line exceeds ATOMS count'
        end if
        type_atom_labels(t,j+1:j+k) = field_atom_label
        j = j + k
      end do

      do
        call next_data_line(uf, l, ios)
        if (ios /= 0) stop 'ERROR searching for FINISH in FIELD'
        if (trim(adjustl(l)) == 'FINISH') exit
      end do

      ! classify type
      is_carbon_type(t) = .false.
      is_wall_type(t) = .false.
      if (natoms_type(t) >= 1) then
         if (trim_array_eq(type_atom_labels(t,1:natoms_type(t)), 'C')) is_carbon_type(t) = .true.
         if (trim_array_eq(type_atom_labels(t,1:natoms_type(t)), 'OW')) is_wall_type(t) = .true.
      end if
      is_gas_type(t) = .not. is_carbon_type(t) .and. .not. is_wall_type(t)
    end do
    close(uf)

    ngas_types = 0
    do t = 1, nmoltypes
      if (is_gas_type(t)) then
        ngas_types = ngas_types + 1
        gas_type_list(ngas_types) = t
      end if
    end do
  end subroutine read_field

  subroutine build_maps()
    integer :: gmol
    total_atoms_expected = 0
    total_molecules_expected = 0
    do t = 1, nmoltypes
      if (nmols(t) < 0 .or. natoms_type(t) <= 0) stop 'ERROR invalid FIELD counts'
      total_atoms_expected = total_atoms_expected + nmols(t)*natoms_type(t)
      total_molecules_expected = total_molecules_expected + nmols(t)
    end do

    allocate(atom_type(total_atoms_expected), atom_mol_global(total_atoms_expected), atom_local_in_mol(total_atoms_expected))
    allocate(mol_type_of_global(total_molecules_expected), &
             mol_first_atom(total_molecules_expected), &
             mol_natoms(total_molecules_expected))
    allocate(is_carbon_atom(total_atoms_expected))

    atom_start_type = 0; atom_end_type = -1; mol_start_type = 0; mol_end_type = -1

    i = 0
    gmol = 0
    do t = 1, nmoltypes
      mol_start_type(t) = gmol + 1
      atom_start_type(t) = i + 1
      do j = 1, nmols(t)
        gmol = gmol + 1
        mol_type_of_global(gmol) = t
        mol_first_atom(gmol) = i + 1
        mol_natoms(gmol) = natoms_type(t)
        do k = 1, natoms_type(t)
          i = i + 1
          atom_type(i) = t
          atom_mol_global(i) = gmol
          atom_local_in_mol(i) = k
          is_carbon_atom(i) = is_carbon_type(t)
        end do
      end do
      mol_end_type(t) = gmol
      atom_end_type(t) = i
    end do
  end subroutine build_maps

  subroutine read_config_header_and_cell()
    open(newunit=uc, file=trim(config_file), status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'ERROR opening CONFIG'
    read(uc,'(A)',iostat=ios) title
    if (ios /= 0) stop 'ERROR reading CONFIG title'
    read(uc,'(A)',iostat=ios) line
    if (ios /= 0) stop 'ERROR reading CONFIG header'
    call parse_header(line, levcfg_cfg, imcon_cfg, natms_cfg)
    if (natms_cfg /= total_atoms_expected) stop 'ERROR FIELD/CONFIG total atoms mismatch'
    do i = 1, 3
      read(uc,*,iostat=ios) cell(i,1), cell(i,2), cell(i,3)
      if (ios /= 0) stop 'ERROR reading CONFIG cell'
    end do
    lx = cell(1,1); ly = cell(2,2); lz = cell(3,3); area = lx*ly
    has_vel_cfg = (levcfg_cfg >= 1)
    has_force_cfg = (levcfg_cfg >= 2)
    close(uc)
  end subroutine read_config_header_and_cell

  subroutine open_history_and_header()
    open(newunit=uh, file=trim(history_file), status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'ERROR opening HISTORY'
    read(uh,'(A)',iostat=ios) title
    if (ios /= 0) stop 'ERROR reading HISTORY title'
    read(uh,'(A)',iostat=ios) line
    if (ios /= 0) stop 'ERROR reading HISTORY header'
    call parse_header(line, levcfg_hist, imcon_hist, natms_hist)
    if (natms_hist /= total_atoms_expected) stop 'ERROR FIELD/HISTORY total atoms mismatch'
    has_vel_hist = (levcfg_hist >= 1)
    has_force_hist = (levcfg_hist >= 2)
  end subroutine open_history_and_header

  subroutine read_next_frame(stepnum, time_ps, ios_out, plane_abs_dummy)
    integer, intent(out) :: stepnum, ios_out
    real(dp), intent(out) :: time_ps, plane_abs_dummy
    real(dp) :: dt_hist
    real(dp) :: xx, yy, zz
    integer :: natms_tmp, keytrj, imcon_tmp

    ios_out = 0
    plane_abs_dummy = 0.0_dp
    do
      read(uh,'(A)',iostat=ios_out) line
      if (ios_out /= 0) return
      if (index(adjustl(line), 'timestep') == 1) exit
    end do

    call parse_timestep_line(line, stepnum, natms_tmp, keytrj, imcon_tmp, dt_hist)
    if (dt_ps <= 0.0_dp) dt_ps = dt_hist
    time_ps = real(stepnum,dp) * dt_ps

    do i = 1, 3
      read(uh,*,iostat=ios_out) cell(i,1), cell(i,2), cell(i,3)
      if (ios_out /= 0) return
    end do
    lx = cell(1,1); ly = cell(2,2); lz = cell(3,3); area = lx*ly

    do i = 1, total_atoms_expected
      read(uh,'(A)',iostat=ios_out) line
      if (ios_out /= 0) return
      read(uh,*,iostat=ios_out) xx, yy, zz
      if (ios_out /= 0) return
      zatom(i) = zz
      if (has_vel_hist) then
        read(uh,'(A)',iostat=ios_out) line
        if (ios_out /= 0) return
      end if
      if (has_force_hist) then
        read(uh,'(A)',iostat=ios_out) line
        if (ios_out /= 0) return
      end if
    end do
  end subroutine read_next_frame

  subroutine compute_carbon_geometry(zc_out, zmin_rel_out, zmax_rel_out, slab_t_out, left_len_out, right_len_out)
    real(dp), intent(out) :: zc_out, zmin_rel_out, zmax_rel_out, slab_t_out, left_len_out, right_len_out
    real(dp) :: csum, ssum, theta, zr
    integer :: nc

    csum = 0.0_dp
    ssum = 0.0_dp
    nc = 0
    do i = 1, total_atoms_expected
      if (is_carbon_atom(i)) then
        theta = twopi() * zatom(i) / lz
        csum = csum + cos(theta)
        ssum = ssum + sin(theta)
        nc = nc + 1
      end if
    end do
    if (nc == 0) stop 'ERROR: no carbon atoms detected from FIELD'
    theta = atan2(ssum, csum)
    if (theta < 0.0_dp) theta = theta + twopi()
    zc_out = lz * theta / twopi()

    zmin_rel_out = hugev
    zmax_rel_out = -hugev
    do i = 1, total_atoms_expected
      if (is_carbon_atom(i)) then
        zr = wrap_rel(zatom(i) - zc_out, lz)
        if (zr < zmin_rel_out) zmin_rel_out = zr
        if (zr > zmax_rel_out) zmax_rel_out = zr
      end if
    end do
    slab_t_out = zmax_rel_out - zmin_rel_out
    left_len_out = zmin_rel_out + 0.5_dp*lz
    right_len_out = 0.5_dp*lz - zmax_rel_out
  end subroutine compute_carbon_geometry

  subroutine compute_all_mol_com(zc_local)
    real(dp), intent(in) :: zc_local
    real(dp) :: zref, zr, sumr

    do i = 1, total_molecules_expected
      zref = zatom(mol_first_atom(i))
      sumr = 0.0_dp
      do j = 0, mol_natoms(i)-1
        zr = wrap_rel(zatom(mol_first_atom(i)+j) - zref, lz)
        sumr = sumr + zr
      end do
      zcom(i) = wrap01(zref + sumr/real(mol_natoms(i),dp), lz)
    end do
  end subroutine compute_all_mol_com

  subroutine classify_and_flux(time_ps, zc_local, zmin_rel_local, zmax_rel_local, &
       left_len_local, right_len_local, plane_rel_local)
    real(dp), intent(in) :: time_ps, zc_local, zmin_rel_local, zmax_rel_local, left_len_local, right_len_local, plane_rel_local
    integer :: g, mt
    real(dp) :: zrel, zprev_rel, drel, elapsed_ps

    nleft = 0; nslab = 0; nright = 0

    do g = 1, total_molecules_expected
      mt = mol_type_of_global(g)
      if (.not. is_gas_type(mt)) cycle

      zrel = wrap_rel(zcom(g) - zc_local, lz)
      if (zrel < zmin_rel_local) then
        nleft(mt) = nleft(mt) + 1
      else if (zrel > zmax_rel_local) then
        nright(mt) = nright(mt) + 1
      else
        nslab(mt) = nslab(mt) + 1
      end if

      if (time_ps >= analysis_start_ps) then
        if (mol_has_prev(g)) then
          zprev_rel = wrap_rel(zcom_prev(g) - zc_local, lz)
          drel = wrap_rel(zrel - zprev_rel, lz)
          if ((zprev_rel < plane_rel_local) .and. (zprev_rel + drel >= plane_rel_local)) then
            nfwd(mt) = nfwd(mt) + 1
            nfwd_tot = nfwd_tot + 1
          else if ((zprev_rel > plane_rel_local) .and. (zprev_rel + drel <= plane_rel_local)) then
            nbwd(mt) = nbwd(mt) + 1
            nbwd_tot = nbwd_tot + 1
          end if
        end if
      end if

      zcom_prev(g) = zcom(g)
      mol_has_prev(g) = .true.
    end do

    nleft_tot = 0; nslab_tot = 0; nright_tot = 0
    do t = 1, nmoltypes
      if (.not. is_gas_type(t)) cycle
      nleft_tot = nleft_tot + nleft(t)
      nslab_tot = nslab_tot + nslab(t)
      nright_tot = nright_tot + nright(t)
      rho_left(t) = real(nleft(t),dp)/(area*left_len_local)
      rho_right(t) = real(nright(t),dp)/(area*right_len_local)
      delta_rho(t) = rho_left(t) - rho_right(t)
    end do

    rho_left_tot = real(nleft_tot,dp)/(area*left_len_local)
    rho_right_tot = real(nright_tot,dp)/(area*right_len_local)
    delta_rho_tot = rho_left_tot - rho_right_tot

    elapsed_ps = max(time_ps - analysis_start_ps, 0.0_dp)
    do t = 1, nmoltypes
      if (.not. is_gas_type(t)) cycle
      if (elapsed_ps > 0.0_dp) then
        flux(t) = real(nfwd(t)-nbwd(t),dp)/(area*elapsed_ps)
      else
        flux(t) = 0.0_dp
      end if
      if (abs(delta_rho(t)) > 1.0d-30) then
        deff_signed(t) = -(flux(t))* (zmax_rel_local-zmin_rel_local) / delta_rho(t)
        deff_abs(t) = abs(deff_signed(t))
      else
        deff_signed(t) = hugev
        deff_abs(t) = hugev
      end if
    end do

    if (elapsed_ps > 0.0_dp) then
      flux_tot = real(nfwd_tot-nbwd_tot,dp)/(area*elapsed_ps)
    else
      flux_tot = 0.0_dp
    end if
    if (abs(delta_rho_tot) > 1.0d-30) then
      deff_signed_tot = -(flux_tot)*(zmax_rel_local-zmin_rel_local)/delta_rho_tot
      deff_abs_tot = abs(deff_signed_tot)
    else
      deff_signed_tot = hugev
      deff_abs_tot = hugev
    end if
  end subroutine classify_and_flux

  subroutine write_output_header(u)
    integer, intent(in) :: u
    write(u,'(A)', advance='no') '# step time_ps slab_center_z slab_thickness left_len right_len plane_z '
    do k = 1, ngas_types
      t = gas_type_list(k)
      write(u,'(A)', advance='no') trim(molname(t))//'_Nleft '//trim(molname(t))//'_Nslab '//trim(molname(t))//'_Nright '
      write(u,'(A)', advance='no') trim(molname(t))//'_rhoL '//trim(molname(t))//'_rhoR '//trim(molname(t))//'_deltaRho '
      write(u,'(A)', advance='no') trim(molname(t))//'_Nfwd '//trim(molname(t))//'_Nbwd '//trim(molname(t))//'_flux '
      write(u,'(A)', advance='no') trim(molname(t))//'_Deff_signed '//trim(molname(t))//'_Deff_abs '
    end do
    write(u,'(A)') 'Nleft_tot Nslab_tot Nright_tot rhoL_tot rhoR_tot deltaRho_tot ' // &
         'Nfwd_tot Nbwd_tot flux_tot Deff_signed_tot Deff_abs_tot'
  end subroutine write_output_header

  subroutine write_frame(u, stepnum, time_ps, zc_local, slab_t_local, left_len_local, right_len_local, plane_abs_local)
    integer, intent(in) :: u, stepnum
    real(dp), intent(in) :: time_ps, zc_local, slab_t_local, left_len_local, right_len_local, plane_abs_local
    write(u,'(I10,1X,F12.6,1X,5(ES16.8,1X))', advance='no') &
         stepnum, time_ps, zc_local, slab_t_local, left_len_local, right_len_local, &
         plane_abs_local
    do k = 1, ngas_types
      t = gas_type_list(k)
      write(u,'(3(I8,1X),3(ES16.8,1X),2(I8,1X),3(ES16.8,1X))', advance='no') &
           nleft(t), nslab(t), nright(t), &
           rho_left(t), rho_right(t), delta_rho(t), &
           nfwd(t), nbwd(t), flux(t), deff_signed(t), deff_abs(t)
    end do
    write(u,'(3(I8,1X),3(ES16.8,1X),2(I8,1X),3(ES16.8,1X))') &
         nleft_tot, nslab_tot, nright_tot, &
         rho_left_tot, rho_right_tot, delta_rho_tot, &
         nfwd_tot, nbwd_tot, flux_tot, deff_signed_tot, deff_abs_tot
  end subroutine write_frame


  subroutine parse_field_atom_line(l, label_out, nrep)
    character(len=*), intent(in) :: l
    character(len=32), intent(out) :: label_out
    integer, intent(out) :: nrep
    character(len=max_line) :: tmp
    integer :: ios_local, n1, n2, n3
    real(dp) :: r1, r2

    label_out = ''
    nrep = 1
    tmp = adjustl(l)

    ! Try the common compressed FIELD form:
    ! label mass charge repeat freeze
    read(tmp,*,iostat=ios_local) label_out, r1, r2, n1, n2
    if (ios_local == 0) then
       nrep = n1
       return
    end if

    ! Try a shorter form:
    ! label mass charge repeat
    read(tmp,*,iostat=ios_local) label_out, r1, r2, n1
    if (ios_local == 0) then
       nrep = n1
       return
    end if

    ! Fallback to single atom-site definition:
    ! label ...
    read(tmp,*,iostat=ios_local) label_out
    if (ios_local /= 0) stop 'ERROR parsing FIELD atom definition line'
    nrep = 1
  end subroutine parse_field_atom_line

  subroutine next_data_line(u, l, ios_out)
    integer, intent(in) :: u
    character(len=*), intent(out) :: l
    integer, intent(out) :: ios_out
    do
      read(u,'(A)',iostat=ios_out) l
      if (ios_out /= 0) return
      l = adjustl(l)
      if (len_trim(l) == 0) cycle
      if (l(1:1) == '#') cycle
      return
    end do
  end subroutine next_data_line

  logical function starts_with_word(l, word)
    character(len=*), intent(in) :: l, word
    character(len=max_line) :: tmp
    tmp = adjustl(l)
    starts_with_word = index(tmp, trim(word)) == 1
  end function starts_with_word

  subroutine parse_header(l, levcfg, imcon, natms)
    character(len=*), intent(in) :: l
    integer, intent(out) :: levcfg, imcon, natms
    integer :: ios_loc, d1, d2
    levcfg = 0; imcon = 0; natms = 0
    read(l,*,iostat=ios_loc) levcfg, imcon, natms, d1, d2
    if (ios_loc /= 0) then
      read(l,*,iostat=ios_loc) levcfg, imcon, natms
      if (ios_loc /= 0) stop 'ERROR parsing header'
    end if
  end subroutine parse_header

  subroutine parse_timestep_line(l, stepnum, natms_tmp, keytrj, imcon_tmp, dt_hist)
    character(len=*), intent(in) :: l
    integer, intent(out) :: stepnum, natms_tmp, keytrj, imcon_tmp
    real(dp), intent(out) :: dt_hist
    character(len=32) :: tag
    integer :: ios_loc
    dt_hist = -1.0_dp
    natms_tmp = 0; keytrj = 0; imcon_tmp = 0; stepnum = 0
    read(l,*,iostat=ios_loc) tag, stepnum, natms_tmp, keytrj, imcon_tmp, dt_hist
    if (ios_loc /= 0) then
      read(l,*,iostat=ios_loc) tag, stepnum
      if (ios_loc /= 0) stop 'ERROR parsing timestep line'
    end if
  end subroutine parse_timestep_line

  subroutine parse_atom_label(l, alab)
    character(len=*), intent(in) :: l
    character(len=32), intent(out) :: alab
    integer :: idx, ios_loc
    real(dp) :: a, b, c
    alab = ''
    read(l,*,iostat=ios_loc) alab, idx, a, b, c
    if (ios_loc /= 0) then
      read(l,*,iostat=ios_loc) alab, idx
      if (ios_loc /= 0) stop 'ERROR parsing atom label'
    end if
  end subroutine parse_atom_label

  pure real(dp) function wrap_rel(dz, box)
    real(dp), intent(in) :: dz, box
    wrap_rel = dz - box*dnint(dz/box)
  end function wrap_rel

  pure real(dp) function wrap01(z, box)
    real(dp), intent(in) :: z, box
    wrap01 = modulo(z, box)
    if (wrap01 < 0.0_dp) wrap01 = wrap01 + box
  end function wrap01

  pure real(dp) function twopi()
    twopi = 6.28318530717958647692_dp
  end function twopi

  pure logical function trim_array_eq(arr, val)
    character(len=*), intent(in) :: arr(:)
    character(len=*), intent(in) :: val
    integer :: ii
    trim_array_eq = .true.
    do ii = 1, size(arr)
      if (trim(arr(ii)) /= trim(val)) then
        trim_array_eq = .false.
        exit
      end if
    end do
  end function trim_array_eq

end program nessmd_periodic_flux_analysis

