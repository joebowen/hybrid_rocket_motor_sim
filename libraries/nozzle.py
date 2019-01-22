import math
import copy

import libraries.constants as constants

from libraries.thermodynamic import Thermodynamic


class Nozzle:
    def __init__(self, mixture, nozzle_dims):
        self.original_mixture = mixture

        self.combustion_thermo = Thermodynamic(copy.deepcopy(mixture))
        self.star_thermo = Thermodynamic(copy.deepcopy(mixture))
        self.exit_thermo = Thermodynamic(copy.deepcopy(mixture))

        if nozzle_dims:
            self.throat_dia = nozzle_dims['nozzle_throat_dia_avg']
            self.exit_dia = nozzle_dims['nozzle_exit_dia_avg']
            self.diffuser_len = nozzle_dims['nozzle_diffuser_len_avg']

            self.ideal = False
        else:
            self.ideal = True

    def inlet_velocity(self, total_mass_flow_rate):
        """inlet_velocity = total_mass_flow_rate / (total_inlet_density * post_chamber_area)
        """
    
        inlet_velocity = total_mass_flow_rate / (self.combustion_thermo.density() * constants.post_chamber_area)

        return inlet_velocity

    def inlet_mach(self, total_mass_flow_rate):
        """inlet_mach = inlet_velocity / speed_of_sound()
        """

        inlet_mach = self.inlet_velocity(total_mass_flow_rate) / self.combustion_thermo.speed_of_sound()

        return inlet_mach

    def throat_area(self, total_mass_flow_rate):
        """throat_area = ((inlet_area * inlet_mach) / star_mach) * sqrt(((1 + ((K() - 1) / 2) * star_mach ** 2) / (1 + ((K() - 1) / 2) * inlet_mach ** 2)) ** ((K() + 1) / (K() - 1)))
        """

        throat_area = (
            (
                (
                    constants.inlet_area * self.inlet_mach(total_mass_flow_rate)
                ) / 
                self.star_mach(total_mass_flow_rate)
            ) *
            math.sqrt(
                (
                    (1 + ((self.combustion_thermo.K() - 1) / 2) * self.star_mach(total_mass_flow_rate) ** 2) /
                    (1 + ((self.combustion_thermo.K() - 1) / 2) * self.inlet_mach(total_mass_flow_rate) ** 2)
                ) ** (
                    (self.combustion_thermo.K() + 1) / (self.combustion_thermo.K() - 1)
                )
            )
        )

        return throat_area

    def throat_diameter(self, total_mass_flow_rate):
        if not self.ideal:
            throat_diameter = self.throat_dia
        else:
            throat_diameter = constants.area_to_diameter(self.throat_area(total_mass_flow_rate))

        return throat_diameter

    def naught_temp(self, total_mass_flow_rate):
        """naught_temp =  + ((inlet_velocity ** 2) / (2 * Cp()))
        """

        naught_temp =  (
            self.combustion_thermo.T() + 
            (
                (
                    self.inlet_velocity(total_mass_flow_rate) ** 2
                ) / 
                (
                    2 * self.combustion_thermo.Cp()
                )
            )
        )

        return naught_temp

    def star_temp(self, total_mass_flow_rate):
        """star_temp = (2 * naught_temp) / (K() + 1)
        """

        star_temp = (
            (
                2 * self.naught_temp(total_mass_flow_rate)
            ) / (
                self.combustion_thermo.K() + 1
            )
        )

        return star_temp

    def star_pressure(self):
        """star_pressure = P() * ((2 / (K() + 1)) ** (K() / (K() - 1)))
        """
    
        star_pressure = (
            self.combustion_thermo.P() *
            (
                (
                    2 / (self.combustion_thermo.K() + 1)
                ) ** (
                    self.combustion_thermo.K() / (self.combustion_thermo.K() - 1)
                )
            )
        )

        return star_pressure

    def star_velocity(self, total_mass_flow_rate):
        """star_velocity = math.sqrt(K() * R() * T())
        """

        self.star_thermo.mixture.T = self.star_temp(total_mass_flow_rate).to(constants.ureg.K).magnitude
        self.star_thermo.mixture.P = self.star_pressure().to(constants.ureg.Pa).magnitude

        star_velocity = math.sqrt(
            (
                self.star_thermo.K() * self.star_thermo.R() * self.star_thermo.T()
            ).to_base_units().magnitude
        ) * (constants.ureg.m / constants.ureg.sec)

        return star_velocity

    def star_mach(self, total_mass_flow_rate):
        """star_mach = star_velocity / speed_of_sound
        """

        star_mach = self.star_velocity(total_mass_flow_rate) / self.star_thermo.speed_of_sound()

        return star_mach

    def exit_mach(self):
        """exit_mach = math...
        """

        exit_mach = (
            math.sqrt(2) * 
            math.sqrt(
                (
                    (
                        self.combustion_thermo.P() * 
                        (
                            (self.combustion_thermo.P() / constants.P_exit) ** (-1 / self.star_thermo.K())
                        )
                    ) - constants.P_exit
                ).to_base_units().magnitude
            ) / 
            math.sqrt(
                (
                    (self.star_thermo.K() * constants.P_exit) - constants.P_exit
                ).to_base_units().magnitude
            )
        )

        return exit_mach

    def exit_area(self, total_mass_flow_rate):
        """exit_area = ((throat_area * throat_mach) / exit_mach) * sqrt(((1 + ((K() - 1) / 2) * exit_mach ** 2) / (1 + ((K() - 1) / 2) * throat_mach ** 2)) ** ((K() + 1) / (K() - 1)))
        """

        exit_area = (
            (
                (self.throat_area(total_mass_flow_rate) * constants.throat_mach) / self.exit_mach()
            ) *
            math.sqrt(
                (
                    (1 + ((self.combustion_thermo.K() - 1) / 2) * self.exit_mach() ** 2) /
                    (1 + ((self.combustion_thermo.K() - 1) / 2) * constants.throat_mach ** 2)
                ) ** (
                    (self.combustion_thermo.K() + 1) / (self.combustion_thermo.K() - 1)
                )
            )
        )

        return exit_area

    def exit_diameter(self, total_mass_flow_rate):
        if not self.ideal:
            exit_diameter = self.exit_dia
        else:
            exit_diameter = constants.area_to_diameter(self.exit_area(total_mass_flow_rate))

        return exit_diameter

    def nozzle_diffuser_length(self, total_mass_flow_rate):
        """nozzle_diffuser_length = (exit_diameter - throat_diameter) / (2 * tan(nozzle_angle))
        """
        if not self.ideal:
            nozzle_diffuser_length = self.diffuser_len
        else:
            nozzle_diffuser_length = (
                (
                    self.exit_diameter(total_mass_flow_rate) - self.throat_diameter(total_mass_flow_rate)
                ) / 
                (
                    2 * math.tan(constants.nozzle_angle)
                )
            )

        return nozzle_diffuser_length

    def exit_velocity(self):
        """exit_velocity = exit_mach * exit_thermo.speed_of_sound()
        """

        self.exit_thermo.P = constants.P_exit

        exit_velocity = self.exit_mach() * self.exit_thermo.speed_of_sound()

        return exit_velocity

    def thrust(self, total_mass_flow_rate):
        """thrust = total_mass_flow_rate * (exit_velocity - inlet_velocity)
        """

        thrust = total_mass_flow_rate * (self.exit_velocity() - self.inlet_velocity(total_mass_flow_rate))

        return thrust

    
