import libraries.constants as constants


class Oxidiser:
    def __init__(self, initial_volume, external_temp, time_step):
        self.time_step = time_step
        self.mass = self.n2o_density(external_temp) * initial_volume
        self.mass_flow_rate = 0 * constants.ureg.kg / constants.ureg.sec

    def mass_flow_rate_function(self):
        self.mass -= (constants.injector_mass_flow_rate * constants.num_of_injectors) * self.time_step

        if self.mass > 0 * constants.ureg.kg:
            self.mass_flow_rate = constants.injector_mass_flow_rate * constants.num_of_injectors
        else:
            self.mass_flow_rate = 0 * constants.ureg.kg / constants.ureg.sec

        return self.mass_flow_rate
        
    def n2o_density(self, external_temp):    
        return constants.n2o_density[external_temp.magnitude]
