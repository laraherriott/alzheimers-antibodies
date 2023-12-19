#
# Class to define the model parameters - use as template to generate other parameter files
# note that these parameters are in time units of hour
#

class NoAbParameters:
    def __init__(self):
        self.k_in = 0.0792 # fixed - literature derived value
        self.k_olig_inc = 9.99929547e-01 #0.924
        self.k_olig_sep = 1.95240222e-03 #0
        self.k_plaque_inc = 9.98804313e-01 #7.0E-8 * 360 #
        self.k_plaque_sep = 5.94184817e-03 #0
        self.k_clear_Abeta_plasma = 1.35608656e-04 #0.231
        self.k_clear_Abeta_csf = 1.00692209e-03#0
        self.k_clear_Abeta_brain =  3.02817495e-03 #0.165
        self.k_clear_oligomer_plasma = 1.64901259e-03 #0.231
        self.k_clear_oligomer_csf = 1.18861133e-03 #0
        self.k_clear_oligomer_brain = 1.01256217e-02 #0.165
        self.k_monomer_brain_plasma = 0.00000000e+00 #0.0923
        self.k_monomer_plasma_brain = 9.99999801e-01 #0
        self.k_oligomer_brain_plasma = 2.14885568e-03 #0.0923
        self.k_oligomer_plasma_brain = 9.99994697e-01 #0
        self.k_monomer_brain_csf = 1.55890798e-03 #0.0982
        self.k_monomer_csf_brain = 0
        self.k_oligomer_brain_csf = 1.83697079e-04 #0.0982
        self.k_oligomer_csf_brain = 0
        self.k_monomer_csf_plasma = 9.99634308e-01 #0.065
        self.k_monomer_plasma_csf =  2.03738347e-03 #0
        self.k_oligomer_csf_plasma = 9.99211953e-01 #0.065
        self.k_oligomer_plasma_csf = 1.80688680e-03 #0

class NoAbParameters_prefit:
    def __init__(self):
        self.k_in = 0.0792 # fixed - literature derived value
        self.k_olig_inc = 0.924
        self.k_olig_sep = 0
        self.k_plaque_inc = 0.165
        self.k_plaque_sep = 0
        self.k_clear_Abeta_plasma = 0.231
        self.k_clear_Abeta_csf = 0
        self.k_clear_Abeta_brain =  0.165
        self.k_clear_oligomer_plasma = 0.231
        self.k_clear_oligomer_csf = 0
        self.k_clear_oligomer_brain = 0.165
        self.k_monomer_brain_plasma = 0.0923
        self.k_monomer_plasma_brain = 0
        self.k_oligomer_brain_plasma = 0.0923
        self.k_oligomer_plasma_brain = 0
        self.k_monomer_brain_csf = 0.0982
        self.k_monomer_csf_brain = 0
        self.k_oligomer_brain_csf = 0.0982
        self.k_oligomer_csf_brain = 0
        self.k_monomer_csf_plasma = 0.065
        self.k_monomer_plasma_csf =  0
        self.k_oligomer_csf_plasma = 0.065
        self.k_oligomer_plasma_csf = 0

class OneAbParameters: # ab-bound monomer should have slower clearance than monomer in the plasma
    def __init__(self):
        self.onPP 
        self.k_onPF 
        self.offma0
        self.offma1 
        self.offma2 
        self.k_offPF 
        self.k_clear_complex 
        self.clearance 
        self.k_mAb_csf_plasma 
        self.k_mAb_plasma_csf 
        self.k_mAb_brain_plasma 
        self.k_mAb_plasma_brain 
        self.k_mAb_brain_csf 
        self.k_mAb_csf_brain = 0
        self.k_synth_FcR 
        self.k_clear_FcR 
        self.k_ADCP 