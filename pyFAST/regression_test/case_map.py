"""
Provides the list of cases that can be run.
"""

OPENFAST = {
    "AWT_YFix_WSt": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn14", "elastodyn", "servodyn"],
        "reference_output": "AWT_YFix_WSt.outb",
        "turbine_directory": "AWT27",
        "tags": []
    },
    "AWT_WSt_StartUp_HighSpShutDown": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn"],
        "reference_output": "AWT_WSt_StartUp_HighSpShutDown.outb",
        "turbine_directory": "AWT27",
        "tags": []
    },
    "AWT_YFree_WSt": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn"],
        "reference_output": "AWT_YFree_WSt.outb",
        "turbine_directory": "AWT27",
        "tags": []
    },
    "AWT_YFree_WTurb": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn14", "elastodyn", "servodyn"],
        "reference_output": "AWT_YFree_WTurb.outb",
        "turbine_directory": "AWT27",
        "tags": []
    },
    "AWT_WSt_StartUpShutDown": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn"],
        "reference_output": "AWT_WSt_StartUpShutDown.outb",
        "turbine_directory": "AWT27",
        "tags": []
    },
    "AOC_WSt": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn14", "elastodyn", "servodyn"],
        "reference_output": "AOC_WSt.outb",
        "turbine_directory": "AOC",
        "tags": []
    },
    "AOC_YFree_WTurb": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn"],
        "reference_output": "AOC_YFree_WTurb.outb",
        "turbine_directory": "AOC",
        "tags": []
    },
    "AOC_YFix_WSt": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn"],
        "reference_output": "AOC_YFix_WSt.outb",
        "turbine_directory": "AOC",
        "tags": []
    },
    "UAE_Dnwind_YRamp_WSt": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn14", "elastodyn", "servodyn"],
        "reference_output": "UAE_Dnwind_YRamp_WSt.outb",
        "turbine_directory": "UAE_VI",
        "tags": []
    },
    "UAE_Upwind_Rigid_WRamp_PwrCurve": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn"],
        "reference_output": "UAE_Upwind_Rigid_WRamp_PwrCurve.outb",
        "turbine_directory": "UAE_VI",
        "tags": []
    },
    "WP_VSP_WTurb_PitchFail": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn14", "elastodyn", "servodyn"],
        "reference_output": "WP_VSP_WTurb_PitchFail.outb",
        "turbine_directory": "WP_Baseline",
        "tags": []
    },
    "WP_VSP_ECD": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn"],
        "reference_output": "WP_VSP_ECD.outb",
        "turbine_directory": "WP_Baseline",
        "tags": []
    },
    "WP_VSP_WTurb": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn"],
        "reference_output": "WP_VSP_WTurb.outb",
        "turbine_directory": "WP_Baseline",
        "tags": []
    },
    "WP_Stationary_Linear": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["elastodyn"],
        #  TODO!
        # "reference_output": "WP_Stationary_Linear.outb",
        "turbine_directory": "WP_Baseline",
        "tags": ["linear"]
    },

    "SWRT_YFree_VS_EDG01": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn"],
        "reference_output": "SWRT_YFree_VS_EDG01.outb",
        "turbine_directory": "5MW_Baseline",
        "tags": []
    },
    "SWRT_YFree_VS_EDC01": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn"],
        "reference_output": "SWRT_YFree_VS_EDC01.outb",
        "turbine_directory": "5MW_Baseline",
        "tags": []
    },
    "SWRT_YFree_VS_WTurb": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn14", "elastodyn", "servodyn"],
        "reference_output": "SWRT_YFree_VS_WTurb.outb",
        "turbine_directory": "5MW_Baseline",
        "tags": []
    },
    "5MW_Land_DLL_WTurb": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn"],
        "reference_output": "5MW_Land_DLL_WTurb.outb",
        "turbine_directory": "5MW_Baseline",
        "tags": []
    },
    "5MW_OC3Mnpl_DLL_WTurb_WavesIrr": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn", "hydrodyn", "subdyn"],
        "reference_output": "5MW_OC3Mnpl_DLL_WTurb_WavesIrr.outb",
        "turbine_directory": "5MW_Baseline",
        "tags": ["offshore"]
    },
    "5MW_OC3Trpd_DLL_WSt_WavesReg": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn", "hydrodyn", "subdyn"],
        "reference_output": "5MW_OC3Trpd_DLL_WSt_WavesReg.outb",
        "turbine_directory": "5MW_Baseline",
        "tags": ["offshore"]
    },
    "5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn", "hydrodyn", "subdyn"],
        "reference_output": "5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth.outb",
        "turbine_directory": "5MW_Baseline",
        "tags": ["offshore"]
    },
    "5MW_ITIBarge_DLL_WTurb_WavesIrr": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn14", "elastodyn", "servodyn", "hydrodyn", "map"],
        "reference_output": "5MW_ITIBarge_DLL_WTurb_WavesIrr.outb",
        "turbine_directory": "5MW_Baseline",
        "tags": ["offshore"]
    },
    "5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn", "hydrodyn", "map"],
        "reference_output": "5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti.outb",
        "turbine_directory": "5MW_Baseline",
        "tags": ["offshore"]
    },
    "5MW_OC3Spar_DLL_WTurb_WavesIrr": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn", "hydrodyn", "map"],
        "reference_output": "5MW_OC3Spar_DLL_WTurb_WavesIrr.outb",
        "turbine_directory": "5MW_Baseline",
        "tags": ["offshore"]
    },
    "5MW_OC4Semi_WSt_WavesWN": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "elastodyn", "servodyn", "hydrodyn", "moordyn"],
        "reference_output": "5MW_OC4Semi_WSt_WavesWN.outb",
        "turbine_directory": "5MW_Baseline",
        "tags": ["offshore"]
    },
    "5MW_Land_BD_DLL_WTurb": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "beamdyn", "servodyn"],
        "reference_output": "5MW_Land_BD_DLL_WTurb.outb",
        "turbine_directory": "5MW_Baseline",
        "tags": []
    },
    "5MW_Land_BD_Linear": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["aerodyn", "beamdyn", "servodyn"],
        "reference_output": "5MW_Land_BD_Linear.outb",
        "turbine_directory": "5MW_Baseline",
        "tags": ["linear"]
    },
    "Ideal_Beam_Fixed_Free_Linear": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["beamdyn"],
        "reference_output": "Ideal_Beam_Fixed_Free_Linear.outb",
        "turbine_directory": "5MW_Baseline",
        "tags": ["linear"]
    },
    "Ideal_Beam_Free_Free_Linear": {
        "type": "regression",
        "driver": "openfast",
        "modules": ["beamdyn"],
        "reference_output": "Ideal_Beam_Free_Free_Linear.outb",
        "turbine_directory": "5MW_Baseline",
        "tags": ["linear"]
    }
}

BEAMDYN_DRIVER = {
    "bd_5MW_dynamic": {
        "type": "regression",
        "driver": "beamdyn",
        "modules": ["beamdyn"],
        "reference_output": "bd_driver.out",
        "tags": ["dynamic"]
    },
    "bd_5MW_dynamic_gravity_Az00": {
        "type": "regression",
        "driver": "beamdyn",
        "modules": ["beamdyn"],
        "reference_output": "bd_driver.out",
        "tags": ["dynamic"]
    },
    "bd_5MW_dynamic_gravity_Az90": {
        "type": "regression",
        "driver": "beamdyn",
        "modules": ["beamdyn"],
        "reference_output": "bd_driver.out",
        "tags": ["dynamic"]
    },
    "bd_curved_beam": {
        "type": "regression",
        "driver": "beamdyn",
        "modules": ["beamdyn"],
        "reference_output": "bd_driver.out",
        "tags": ["static"]
    },
    "bd_isotropic_rollup": {
        "type": "regression",
        "driver": "beamdyn",
        "modules": ["beamdyn"],
        "reference_output": "bd_driver.out",
        "tags": ["static"]
    },
    "bd_static_cantilever_beam": {
        "type": "regression",
        "driver": "beamdyn",
        "modules": ["beamdyn"],
        "reference_output": "bd_driver.out",
        "tags": ["static"]
    },
    "bd_static_twisted_with_k1": {
        "type": "regression",
        "driver": "beamdyn",
        "modules": ["beamdyn"],
        "reference_output": "bd_driver.out",
        "tags": ["static"]
    }
}

CASE_MAP = {
    **OPENFAST,
    **BEAMDYN_DRIVER
}
