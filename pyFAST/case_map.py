"""Provides the list of cases that can be run."""

ROF = "openfast_regression"
RBD = "beamdyn_regression"
AD = "aerodyn"
BD = "beamdyn"
LIN = "linear"

CASE_MAP = {
    "AWT_YFix_WSt": [ROF],
    "AWT_WSt_StartUp_HighSpShutDown": [ROF],
    "AWT_YFree_WSt": [ROF],
    "AWT_YFree_WTurb": [ROF],
    "AWT_WSt_StartUpShutDown": [ROF],
    "AOC_WSt": [ROF],
    "AOC_YFree_WTurb": [ROF],
    "AOC_YFix_WSt": [ROF],
    "UAE_Dnwind_YRamp_WSt": [ROF],
    "UAE_Upwind_Rigid_WRamp_PwrCurve": [ROF],
    "WP_VSP_WTurb_PitchFail": [ROF],
    "WP_VSP_ECD": [ROF],
    "WP_VSP_WTurb": [ROF],
    "WP_Stationary_Linear": [ROF, LIN],
    "SWRT_YFree_VS_EDG01": [ROF],
    "SWRT_YFree_VS_EDC01": [ROF],
    "SWRT_YFree_VS_WTurb": [ROF],
    "5MW_Land_DLL_WTurb": [ROF],
    "5MW_OC3Mnpl_DLL_WTurb_WavesIrr": [ROF],
    "5MW_OC3Trpd_DLL_WSt_WavesReg": [ROF],
    "5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth": [ROF],
    "5MW_ITIBarge_DLL_WTurb_WavesIrr": [ROF],
    "5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti": [ROF],
    "5MW_OC3Spar_DLL_WTurb_WavesIrr": [ROF],
    "5MW_OC4Semi_WSt_WavesWN": [ROF],
    "5MW_Land_BD_DLL_WTurb": [ROF],
    "5MW_Land_BD_Linear": [ROF, LIN],
    "Ideal_Beam_Fixed_Free_Linear": [ROF, LIN],
    "Ideal_Beam_Free_Free_Linear": [ROF, LIN],
    "bd_5MW_dynamic": [RBD],
    "bd_5MW_dynamic_gravity_Az00": [RBD],
    "bd_5MW_dynamic_gravity_Az90": [RBD],
    "bd_curved_beam": [RBD],
    "bd_isotropic_rollup": [RBD],
    "bd_static_cantilever_beam": [RBD],
    "bd_static_twisted_with_k1": [RBD]
}
