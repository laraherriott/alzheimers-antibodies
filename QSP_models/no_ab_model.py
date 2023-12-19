from .parameters import NoAbParameters, NoAbParameters_prefit

class NoAbModel:
    def __init__(self):
        self.params = NoAbParameters_prefit()
        self.type = 'no_ab'

    def brain(self, t, y):
        brain_monomer = y[0]
        brain_oligomer = y[1]
        brain_plaque = y[2]

        plasma_monomer = y[3]
        plasma_oligomer = y[4]

        csf_monomer = y[5]
        csf_oligomer = y[6]

        d_brain_monomer = (self.params.k_in - self.params.k_olig_inc*brain_monomer +
                           self.params.k_olig_sep*brain_oligomer - self.params.k_clear_Abeta_brain*brain_monomer -
                           self.params.k_monomer_brain_plasma*brain_monomer + self.params.k_monomer_plasma_brain*plasma_monomer -
                           self.params.k_monomer_brain_csf*brain_monomer + self.params.k_monomer_csf_brain*csf_monomer)
        d_brain_oligomer = (self.params.k_olig_inc*brain_monomer - self.params.k_olig_sep*brain_oligomer -
                            self.params.k_plaque_inc*brain_oligomer + self.params.k_plaque_sep*brain_plaque -
                            self.params.k_clear_oligomer_brain*brain_oligomer -
                            self.params.k_oligomer_brain_plasma*brain_oligomer + self.params.k_oligomer_plasma_brain*plasma_oligomer -
                            self.params.k_oligomer_brain_csf*brain_oligomer + self.params.k_oligomer_csf_brain*csf_oligomer)
        d_brain_plaque = (self.params.k_plaque_inc*brain_oligomer - self.params.k_plaque_sep*brain_plaque)

        d_plasma_monomer = (self.params.k_monomer_brain_plasma*brain_monomer - self.params.k_monomer_plasma_brain*plasma_monomer +
                            self.params.k_monomer_csf_plasma*csf_monomer - self.params.k_monomer_plasma_csf*plasma_monomer -
                            self.params.k_olig_inc*plasma_monomer + self.params.k_olig_sep*plasma_oligomer - self.params.k_clear_Abeta_plasma*plasma_monomer)
        d_plasma_oligomer = (self.params.k_oligomer_brain_plasma*brain_oligomer - self.params.k_oligomer_plasma_brain*plasma_oligomer +
                             self.params.k_oligomer_csf_plasma*csf_oligomer - self.params.k_oligomer_plasma_csf*plasma_oligomer +
                             self.params.k_olig_inc*plasma_monomer - self.params.k_olig_sep*plasma_oligomer - self.params.k_clear_oligomer_plasma*plasma_oligomer)

        d_csf_monomer = (self.params.k_monomer_brain_csf*brain_monomer - self.params.k_monomer_csf_brain*csf_monomer -
                         self.params.k_monomer_csf_plasma*csf_monomer + self.params.k_monomer_plasma_csf*plasma_monomer -
                         self.params.k_olig_inc*csf_monomer + self.params.k_olig_sep*csf_oligomer - self.params.k_clear_Abeta_csf*csf_monomer)
        d_csf_oligomer = (self.params.k_oligomer_brain_csf*brain_oligomer - self.params.k_oligomer_csf_brain*csf_oligomer -
                          self.params.k_oligomer_csf_plasma*csf_oligomer + self.params.k_oligomer_plasma_csf*plasma_oligomer +
                          self.params.k_olig_inc*plasma_monomer - self.params.k_olig_sep*plasma_oligomer - self.params.k_clear_oligomer_csf*csf_oligomer)

        dYdt = [d_brain_monomer, d_brain_oligomer, d_brain_plaque, d_plasma_monomer, d_plasma_oligomer, d_csf_monomer, d_csf_oligomer]
        
        return dYdt
