from dataclasses import dataclass

@dataclass
class WarningMessage:
    AREA_ESTIMATION_FAIL = """
    <<<< Area Estimation Fail >>>> \n
    The net areas for the x,y,z directions were not 
    estimated. This error will cause that it will 
    not be possible to generate RCGF as the 
    induced currents are weighted by the areas to the 
    total induced magnetic field estimation. \n
    However, it doesn't kill the process and one 
    can continue to work with the molecule or 
    with estimating induced magnetic field for 
    test/debugging purposes. \n
    The error is likely happening due to Polygon 
    area estimation. Check: \n 
    rcm.py --> get_net_area_weight() \n
    area.py --> area_3d() or _area()
    """


@dataclass
class LoadingMessage:
    LOADING_XYZ = """
    >> >> >> loading xyz file and converting it 
    to rdkit.Chem.mol object.
    """

    LOADING_CONN = """
    >> >> >> loading connectivity file.
    """

