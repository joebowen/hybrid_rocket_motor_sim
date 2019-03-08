import libraries.constants as constants


class Oxidiser:
    def __init__(self, external_temp):
        self.mass = self.n2o_density(external_temp) * constants.initial_oxidiser_volume
        self.mass_flow_rate = 0 * constants.ureg.kg / constants.ureg.sec

    def mass_flow_rate_function(self):
        self.mass -= (constants.injector_mass_flow_rate * constants.num_of_injectors) * constants.time_step

        if self.mass > 0 * constants.ureg.kg:
            self.mass_flow_rate = constants.injector_mass_flow_rate * constants.num_of_injectors
        else:
            self.mass_flow_rate = 0 * constants.ureg.kg / constants.ureg.sec

        return self.mass_flow_rate
        
    def n2o_density(self, external_temp):    
        return constants.n2o_density[external_temp.to(constants.ureg.degF).magnitude]
