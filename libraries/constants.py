import math
import pint

ureg = pint.UnitRegistry()

iteration_precision = 0.01

throat_mach = 1

nozzle_angle = 15 * ureg.degree

P_exit = (1 * ureg.atmosphere).to_base_units()

initial_oxidiser_volume = 0.4139 * ureg.L

initial_port_diameter = 1.0 * ureg.inches

port_length = 13.4 * ureg.inches

inlet_dia = 1.45 * ureg.inches 

time_step = .001 * ureg.sec

# Currently based on average of RattWorks H70, I80, and I90 (N20 weight / burntime)
# injector_mass_flow_rate = 0.042 * ureg.kg / ureg.sec

# Based on hybrid motor test #1 (8.75 sec burn time)
# injector_mass_flow_rate = 0.068 * ureg.kg / ureg.sec

# Based on hybrid motor test #2 (11.77 sec burn time)
injector_mass_flow_rate = 0.0274 * ureg.kg / ureg.sec

num_of_injectors = 1
grain_diameter = 1.75 * ureg.inches

a = 0.05 * (ureg.m ** 2) / ureg.kg
n = 0.65
m = -0.2

fuel_density = 3.957 * ureg.kg / (ureg.m ** 3)

ideal_OF_ratio = 4.83

motor_codes = {
    '1/8A': 0.3125 * ureg.newton * ureg.second,
    '1/4A': 0.625 * ureg.newton * ureg.second,
    '1/2A': 1.25 * ureg.newton * ureg.second,
    'A': 2.5 * ureg.newton * ureg.second,
    'B': 5 * ureg.newton * ureg.second,
    'C': 10 * ureg.newton * ureg.second,
    'D': 20 * ureg.newton * ureg.second,
    'E': 40 * ureg.newton * ureg.second,
    'F': 80 * ureg.newton * ureg.second,
    'G': 160 * ureg.newton * ureg.second,
    'H': 320 * ureg.newton * ureg.second,
    'I': 640 * ureg.newton * ureg.second,
    'J': 1280 * ureg.newton * ureg.second,
    'K': 2560 * ureg.newton * ureg.second,
    'L': 5120 * ureg.newton * ureg.second,
    'M': 10240 * ureg.newton * ureg.second,
    'N': 20480 * ureg.newton * ureg.second,
    'O': 40960 * ureg.newton * ureg.second
}


def get_motor_code(impulse):
    for code, impulse_limit in motor_codes.items():
        if impulse < impulse_limit:
            return code
    
    return '>O'


def area_to_diameter(area):
    """diameter = sqrt(4 * area / pi)
    """

    diameter = math.sqrt((4 * area / math.pi).to_base_units().magnitude) * ureg.m

    return diameter


def diameter_to_area(diameter):
    """area = pi * diameter ** 2
    """

    area = math.pi * (diameter ** 2)

    return area


inlet_area = diameter_to_area(inlet_dia)
