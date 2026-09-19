target_sources(
  ${main_lib}
  PRIVATE mqc_fraglist.f90
          mqc_capi_status.f90
          mqc_capi_fraglist.f90
          mqc_capi_system.f90
          mqc_capi_bond_orders.F90
          mqc_validate.f90
          mqc_session.F90
          mqc_capi_session.f90
          mqc_capi_run.f90
          mqc_calculation_interface.f90)

# Charges and Mayer bond orders both need a real density, so they need the
# integrals backend. Without libcint there is nothing to partition and the
# symbols are simply absent -- a Python caller gets an AttributeError on the
# binding rather than a link error here. The xTB bond orders are unaffected:
# they are in mqc_capi_bond_orders.F90 above, which needs no density.
if(MQC_ENABLE_CZT)
  target_sources(${main_lib} PRIVATE mqc_capi_charges.f90 mqc_capi_mayer.f90)
endif()
