import math

import libraries.constants as constants


class Combustion:
    def __init__(self, oxidiser, time_step):
        self.average_port_diameter = constants.initial_port_diameter
        self.oxidiser = oxidiser
        self.time_step = time_step

        self.average_total_mass_flow_rate = self.oxidiser.mass_flow_rate_function()

    def solve_for_average_total_mass_flow_rate(self, iteration_precision):
        oxi_mass_flow_rate = self.oxidiser.mass_flow_rate_function()

        self.average_total_mass_flow_rate = oxi_mass_flow_rate

        while True:
            temp_average_total_mass_flow_rate = self.average_total_mass_flow_rate

            self.average_total_mass_flow_rate = self.average_fuel_mass_flow_rate_function() + oxi_mass_flow_rate

            if abs(self.average_total_mass_flow_rate - temp_average_total_mass_flow_rate).magnitude < iteration_precision:
                break

        self.average_port_diameter += 2 * self.average_regression_rate * self.time_step / 100

        if self.average_port_diameter.to_base_units().magnitude > constants.grain_diameter.to_base_units().magnitude:
            print('Motor burn through! Fuel grain too thin!')
            raise ValueError('Motor burn through! Fuel grain too thin!')

        return self.average_total_mass_flow_rate

    def average_total_mass_flux_function(self):
        """average_total_mass_flux = average_total_mass_flow_rate / average_port_surface_area
        """

        self.average_total_mass_flux = self.average_total_mass_flow_rate / self.average_port_surface_area()

        return self.average_total_mass_flux

    def average_regression_rate_function(self):
        """average_regression_rate = a * (average_mass_flux ** n) * (port_length ** m)
        """

        average_total_mass_flux = self.average_total_mass_flux_function().to_base_units()
        average_total_mass_flux_units = average_total_mass_flux.units

        port_length = constants.port_length.to_base_units()
        port_length_units = port_length.units

        self.average_regression_rate = constants.a * (average_total_mass_flux.magnitude ** constants.n) * (port_length.magnitude ** constants.m)

        # Now fix the units
        self.average_regression_rate = self.average_regression_rate * port_length_units * average_total_mass_flux_units

        return self.average_regression_rate

    def average_fuel_mass_flow_rate_function(self):
        """average_fuel_mass_flow_rate = fuel_density * average_regression_rate * port_surface_area
        """
        self.average_fuel_mass_flow_rate = constants.fuel_density * self.average_regression_rate_function() * self.average_port_surface_area()

        return self.average_fuel_mass_flow_rate

    def average_port_surface_area(self):
        """average_port_surface_area = self.average_port_diameter * math.pi *self.port_length
        """

        average_port_surface_area = self.average_port_diameter * math.pi * constants.port_length

        return average_port_surface_area

    def oxi_fuel_ratio_function(self):
        """oxi_fuel_ratio = oxidiser.mass_flow_rate / average_fuel_mass_flow_rate
        """

        try:
            self.oxi_fuel_ratio = self.oxidiser.mass_flow_rate / self.average_fuel_mass_flow_rate
        except ZeroDivisionError:
            pass

        return self.oxi_fuel_ratio
