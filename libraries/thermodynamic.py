import math

import libraries.constants as constants


class Thermodynamic:
    def __init__(self, mixture):
        self.mixture = mixture

    def K(self):
        return self.mixture.Cpg / self.mixture.Cvg

    def P(self):
        return (self.mixture.P * constants.ureg.Pa).to_base_units()

    def R(self):
        return (self.mixture.R_specific * (constants.ureg.J / constants.ureg.kg / constants.ureg.K)).to_base_units()

    def T(self):
        return (self.mixture.T * constants.ureg.K).to_base_units()

    def Cp(self):
        return (self.mixture.Cpg * (constants.ureg.J / constants.ureg.kg / constants.ureg.K)).to_base_units()

    def density(self):
        return (self.mixture.rhog * (constants.ureg.kg / (constants.ureg.m ** 3))).to_base_units()

    def speed_of_sound(self):
        """speed_of_sound = sqrt(K(temp) * R() * temp)
        """

        speed_of_sound = math.sqrt((self.K() * self.R() * self.T()).magnitude) * (constants.ureg.m / constants.ureg.sec)

        return speed_of_sound
