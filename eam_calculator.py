from ase.calculators.eam import EAM
def set_up_eam_calculator() -> EAM:
    calculator = EAM(potential = "./Zr_2.eam.fs")
    return calculator