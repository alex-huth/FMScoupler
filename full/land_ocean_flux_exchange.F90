!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS) Coupler.
!*
!* FMS Coupler is free software: you can redistribute it and/or modify
!* it under the terms of the GNU Lesser General Public License as
!* published by the Free Software Foundation, either version 3 of the
!* License, or (at your option) any later version.
!*
!* FMS Coupler is distributed in the hope that it will be useful, but
!* WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!* General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS Coupler.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> \file
!> \brief Handles flux exchanges and exchange grids between land and ocean grids
module land_ocean_flux_exchange_mod

!! FMS
  use FMS
  use FMSconstants, only: RADIUS
!! Components
  use land_model_mod,      only: land_data_type
  use ocean_model_mod,     only: ocean_public_type, land_ocean_boundary_type

  implicit none
  private


  !---- exchange grid maps -----

  type(FmsXgridXmap_type), save :: xmap_IS_smb
  integer         :: n_xgrid_IS_smb=0

  ! Exchange grid indices
  integer :: X2_GRID_LND, X2_GRID_OCN

  public :: flux_land_to_ocean, land_ocean_flux_exchange_init

  integer :: cplClock, fluxLandOceanClock
  logical :: do_IS_smb
  real    :: Dt_cpl
contains

  subroutine land_ocean_flux_exchange_init(Land, Ocean, land_ocean_boundary, Dt_cpl_in, do_IS_smb_in, cplClock_in)
    type(land_data_type),         intent(in)    :: Land !< A derived data type to specify land boundary data
    type(ocean_public_type),          intent(inout) :: Ocean !< A derived data type to specify ocean boundary data
    type(land_ocean_boundary_type), intent(inout) :: land_ocean_boundary !< A derived data type to specify properties
                                                                         !! and fluxes passed from land to ocean
    real,                         intent(in)    :: Dt_cpl_in
    logical,                      intent(in)    :: do_IS_smb_in
    integer,                      intent(in)    :: cplClock_in

    integer :: is, ie, js, je

    do_IS_smb = do_IS_smb_in
    cplClock = cplClock_in
    Dt_cpl   = Dt_cpl_in
    fluxLandOceanClock = fms_mpp_clock_id( 'Flux land to ice sheet', flags=fms_clock_flag_default, grain=CLOCK_ROUTINE )

    if (do_IS_smb) then
       call fms_xgrid_setup_xmap(xmap_IS_smb, (/ 'LND', 'OCN' /),       &
            (/ Land%Domain, Ocean%Domain /),                    &
            "INPUT/grid_spec.nc"             )
       ! exchange grid indices
       X2_GRID_LND = 1; X2_GRID_OCN = 2;
       n_xgrid_IS_smb = max(fms_xgrid_count(xmap_IS_smb),1)
       if (n_xgrid_IS_smb.eq.1) write (*,'(a,i6,6x,a)') 'PE = ', fms_mpp_pe(), 'ice sheet SMB  exchange size equals one.'
    endif

    call fms_mpp_domains_get_compute_domain( Ocean%domain, is, ie, js, je )

    !allocate land_ocean_boundary
    allocate( land_ocean_boundary%shelf_sfc_mass_flux(is:ie,js:je) )
    allocate( land_ocean_boundary%shelf_sfc_mass_hflx(is:ie,js:je) )
    ! initialize values for override experiments (mjh)
    land_ocean_boundary%shelf_sfc_mass_flux=0.0
    land_ocean_boundary%shelf_sfc_mass_hflx=0.0


  end subroutine land_ocean_flux_exchange_init

  !#######################################################################
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! flux_land_to_ocean - translate shelf_sfc_mass_flux from land to ocean grids  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  !> \brief Conservative transfer of surface mass flux from the land model to the ocean model.
  !!
  !! The following elements are transferred from the Land to the land_ocean_boundary:
  !! <pre>
  !!        IS_smb --> shelf_sfc_mass_flux (kg/m2)
  !! </pre>
  subroutine flux_land_to_ocean( Time, Land, Ocean, Land_Ocean_Boundary )
    type(FmsTime_type),             intent(in) :: Time !< Current time
    type(land_data_type),           intent(in) :: Land !< A derived data type to specify land boundary data
    type(ocean_public_type),          intent(in) :: Ocean !< A derived data type to specify ocean boundary data
    type(land_ocean_boundary_type), intent(inout):: Land_Ocean_Boundary !< A derived data type to specify properties and fluxes passed
                                                                        !! from land to ocean

    integer                         :: ier
    real, dimension(n_xgrid_IS_smb) :: ex_shelf_sfc_mass_flux, ex_shelf_sfc_mass_hflx
    real, dimension(size(Land_Ocean_Boundary%shelf_sfc_mass_flux,1),size(Land_Ocean_Boundary%shelf_sfc_mass_flux,2),1) :: ocn_buf

    !Balaji
    call fms_mpp_clock_begin(cplClock)
    call fms_mpp_clock_begin(fluxLandOceanClock)

    ! ccc = conservation_check(Land%IS_smb, 'LND', xmap_IS_smb)
    ! if (fms_mpp_pe()==fms_mpp_root_pe()) print *,'SHELF_SFC_MASS_FLUX', ccc

    if (do_IS_smb) then
       Land%IS_smb = 0.1 !for testing
       call fms_xgrid_put_to_xgrid ( Land%IS_smb,      'LND', ex_shelf_sfc_mass_flux,  xmap_IS_smb)
       call fms_xgrid_put_to_xgrid ( Land%IS_smb_heat,      'LND', ex_shelf_sfc_mass_hflx,  xmap_IS_smb)
       call fms_xgrid_get_from_xgrid (ocn_buf, 'OCN', ex_shelf_sfc_mass_flux,  xmap_IS_smb)
       Land_Ocean_Boundary%shelf_sfc_mass_flux = ocn_buf(:,:,1);
       call fms_xgrid_get_from_xgrid (ocn_buf, 'OCN', ex_shelf_sfc_mass_hflx,  xmap_IS_smb)
       Land_Ocean_Boundary%shelf_sfc_mass_hflx = ocn_buf(:,:,1);
       !Balaji
       call fms_data_override('OCN', 'shelf_sfc_mass_flux' , Land_Ocean_Boundary%shelf_sfc_mass_flux , Time)
       call fms_data_override('OCN', 'shelf_sfc_mass_hflx' , Land_Ocean_Boundary%shelf_sfc_mass_hflx , Time)

       ! compute stock increment
       ocn_buf(:,:,1) = Land_Ocean_Boundary%shelf_sfc_mass_flux
       call fms_xgrid_stock_move(from=fms_stock_constants_lnd_stock(ISTOCK_WATER), to=fms_stock_constants_ocn_stock(ISTOCK_WATER), &
            & grid_index=X2_GRID_OCN, &
            & stock_data3d=ocn_buf, &
            & xmap=xmap_IS_smb, &
            & delta_t=Dt_cpl, &
            & from_side=ISTOCK_SIDE, to_side=ISTOCK_SIDE, &
            & radius=Radius, ier=ier, verbose='stock move SHELF_SFC_MASS_FLUX (Lnd->Ocn) ')
    else
       Land_Ocean_Boundary%shelf_sfc_mass_flux = 0.0
       Land_Ocean_Boundary%shelf_sfc_mass_hflx = 0.0
    endif

    call fms_mpp_clock_end(fluxLandOceanClock)
    call fms_mpp_clock_end(cplClock)

  end subroutine flux_land_to_ocean


!#######################################################################

end module land_ocean_flux_exchange_mod
