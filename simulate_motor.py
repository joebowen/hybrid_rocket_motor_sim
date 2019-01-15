import pandas as pd 
import numpy as np

import matplotlib.pyplot as plt

import libraries.constants as constants

from fpdf import FPDF

from libraries.combustion import Combustion
from libraries.oxidiser import Oxidiser
from libraries.nozzle import Nozzle
from libraries.thermochemical import Thermochemical

initial_oxidiser_volume = 1 * constants.ureg.L

external_temp = constants.ureg.Quantity(85, constants.ureg.degF)

initial_port_diameter = 1.2 * constants.ureg.inches
port_length = 15 * constants.ureg.inches

time_step = .001 * constants.ureg.sec

pd.set_option('display.max_columns', 500)


def simulate():
    oxidiser = Oxidiser(initial_oxidiser_volume, external_temp, time_step)
    combustion = Combustion(initial_port_diameter, port_length, oxidiser, time_step)
    thermochemical = Thermochemical()
    ideal_nozzle = Nozzle(thermochemical.get_mixture())
    
    time = 0 * constants.ureg.sec

    count = 0

    impulse = 0 * constants.ureg.newton * constants.ureg.second

    raw_data = []

    try:
        while combustion.average_total_mass_flow_rate > 0 * constants.ureg.kg / constants.ureg.sec:
            combustion.solve_for_average_total_mass_flow_rate(iteration_precision=0.1)

            average_total_mass_flow_rate = combustion.average_total_mass_flow_rate.to_base_units()

            raw_data.append(
                {
                    'time ({0.units})'.format(time.to_base_units()): time.to_base_units().magnitude,
                    'average total mass flow rate ({0.units})'.format(average_total_mass_flow_rate.to(constants.ureg.lbs / constants.ureg.second)):round(average_total_mass_flow_rate.to(constants.ureg.lbs / constants.ureg.seconds).magnitude, 3),
                    'average port diameter ({0.units})'.format(combustion.average_port_diameter.to(constants.ureg.inches)): round(combustion.average_port_diameter.to(constants.ureg.inches).magnitude, 4),
                    'average regression rate ({0.units})'.format(combustion.average_regression_rate.to(constants.ureg.inches / constants.ureg.second)): round(combustion.average_regression_rate.to(constants.ureg.inches / constants.ureg.second).magnitude / 100, 5),
                    'average fuel mass flow rate ({0.units})'.format(combustion.average_fuel_mass_flow_rate.to(constants.ureg.lbs / constants.ureg.second)): round(combustion.average_fuel_mass_flow_rate.to(constants.ureg.lbs / constants.ureg.second).magnitude, 7),
                    'oxidiser mass flow rate ({0.units})'.format(oxidiser.mass_flow_rate.to(constants.ureg.lbs / constants.ureg.second)): round(oxidiser.mass_flow_rate.to(constants.ureg.lbs / constants.ureg.second).magnitude, 4),
                    'average total mass flux ({0.units})'.format(combustion.average_total_mass_flux.to(constants.ureg.lbs / ((constants.ureg.inches **2) * constants.ureg.second))): round(combustion.average_total_mass_flux.to(constants.ureg.lbs / ((constants.ureg.inches **2) * constants.ureg.second)).magnitude, 4),
                    'oxi fuel ratio': round(combustion.oxi_fuel_ratio_function().to_base_units().magnitude, 4),
                    'inlet velocity ({0.units})'.format(ideal_nozzle.inlet_velocity(average_total_mass_flow_rate).to(constants.ureg.mph)): round(ideal_nozzle.inlet_velocity(average_total_mass_flow_rate).to(constants.ureg.mph).magnitude, 4),
                    'inlet mach': round(ideal_nozzle.inlet_mach(average_total_mass_flow_rate).magnitude, 4),
                    'ideal nozzle throat area ({0.units})'.format(ideal_nozzle.throat_area(average_total_mass_flow_rate).to(constants.ureg.inches ** 2)): round(ideal_nozzle.throat_area(average_total_mass_flow_rate).to(constants.ureg.inches ** 2).magnitude, 4),
                    'ideal nozzle throat diameter ({0.units})'.format(ideal_nozzle.throat_diameter(average_total_mass_flow_rate).to(constants.ureg.inches)): round(ideal_nozzle.throat_diameter(average_total_mass_flow_rate).to(constants.ureg.inches).magnitude, 4),
                    'ideal nozzle naught temp ({0.units})'.format(ideal_nozzle.naught_temp(average_total_mass_flow_rate).to(constants.ureg.degF)): round(ideal_nozzle.naught_temp(average_total_mass_flow_rate).to(constants.ureg.degF).magnitude, 4),
                    'ideal nozzle star temp ({0.units})'.format(ideal_nozzle.star_temp(average_total_mass_flow_rate).to(constants.ureg.degF)): round(ideal_nozzle.star_temp(average_total_mass_flow_rate).to(constants.ureg.degF).magnitude, 4),
                    'ideal nozzle star pressure ({0.units})'.format(ideal_nozzle.star_pressure().to(constants.ureg.psi)): round(ideal_nozzle.star_pressure().to(constants.ureg.psi).magnitude, 4),
                    'ideal nozzle star velocity ({0.units})'.format(ideal_nozzle.star_velocity(average_total_mass_flow_rate).to(constants.ureg.mph)): round(ideal_nozzle.star_velocity(average_total_mass_flow_rate).to(constants.ureg.mph).magnitude, 4),
                    'ideal nozzle star mach': round(ideal_nozzle.star_mach(average_total_mass_flow_rate).to_base_units().magnitude, 4),
                    'ideal nozzle exit mach': round(ideal_nozzle.exit_mach(), 4),
                    'ideal nozzle exit area ({0.units})'.format(ideal_nozzle.exit_area(average_total_mass_flow_rate).to(constants.ureg.inches ** 2)): round(ideal_nozzle.exit_area(average_total_mass_flow_rate).to(constants.ureg.inches ** 2).magnitude, 4),
                    'ideal nozzle exit diameter ({0.units})'.format(ideal_nozzle.exit_diameter(average_total_mass_flow_rate).to(constants.ureg.inches)): round(ideal_nozzle.exit_diameter(average_total_mass_flow_rate).to(constants.ureg.inches).magnitude, 4),
                    'ideal nozzle nozzle diffuser length ({0.units})'.format(ideal_nozzle.nozzle_diffuser_length(average_total_mass_flow_rate).to(constants.ureg.inches)): round(ideal_nozzle.nozzle_diffuser_length(average_total_mass_flow_rate).to(constants.ureg.inches).magnitude, 4),
                    'ideal nozzle exit velocity ({0.units})'.format(ideal_nozzle.exit_velocity().to(constants.ureg.mph)): round(ideal_nozzle.exit_velocity().to(constants.ureg.mph).magnitude, 4),
                    'ideal nozzle thrust ({0.units})'.format(ideal_nozzle.thrust(average_total_mass_flow_rate).to(constants.ureg.newton)): round(ideal_nozzle.thrust(average_total_mass_flow_rate).to(constants.ureg.newton).magnitude, 4)
                }
            )

            impulse += ideal_nozzle.thrust(average_total_mass_flow_rate).to(constants.ureg.newton) * time_step

            time += time_step

            count += 1
    except Exception as e:
        print(e)

        raise e
    finally:
        data = pd.DataFrame(raw_data)
        
        data.set_index(
            'time (second)',
            drop=True,
            inplace=True
        )

        data.to_csv('results/motor_data.csv')

        pdf = FPDF()

        pdf.add_page('P')

        pdf.set_font('Arial', '', 14)

        pdf.write(10, 'Simulation Inputs:\n')

        pdf.write(5, 'a: {0:.4f}\n'.format(constants.a.to((constants.ureg.inches ** 2)/constants.ureg.lbs)))
        pdf.write(5, 'n: {}\n'.format(constants.n))
        pdf.write(5, 'm: {}\n'.format(constants.m))
        
        pdf.write(5, 'Initial Oxidiser Volume: {}\n'.format(initial_oxidiser_volume.to(constants.ureg.liter)))
        pdf.write(5, 'External Temp: {}\n'.format(external_temp.to(constants.ureg.degF)))

        pdf.write(5, 'Grain Diameter: {0:.4f}\n'.format(constants.grain_diameter.to(constants.ureg.inches)))
        pdf.write(5, 'Initial Port Diameter: {}\n'.format(initial_port_diameter.to(constants.ureg.inches)))
        pdf.write(5, 'Port Length: {}\n'.format(port_length.to(constants.ureg.inches)))
        pdf.write(5, 'Fuel Density: {0:.6f}\n'.format(constants.fuel_density.to(constants.ureg.lbs / (constants.ureg.inches ** 3))))
        
        pdf.write(5, 'Injector Mass Flow Rate: {0:.4f}\n'.format(constants.injector_mass_flow_rate.to(constants.ureg.lbs / constants.ureg.second)))
        pdf.write(5, 'Number of Injectors: {}\n'.format(constants.num_of_injectors))
        
        pdf.write(5, 'Ideal O/F Ratio: {}\n'.format(constants.ideal_OF_ratio))

        pdf.write(5, 'Time Step: {}\n'.format(time_step.to_base_units()))

        pdf.write(10, 'Simulation Results:\n')

        pdf.write(5, 'Total Burn Time: {}\n'.format(round(time, 3)))
        pdf.write(5, 'Impulse: {}\n'.format(round(impulse, 2)))

        average_trust = data['ideal nozzle thrust (newton)'].mean()

        pdf.write(5, 'Average Thrust: {} newton\n'.format(round(average_trust, 2)))

        motor_code = constants.get_motor_code(impulse)

        pdf.write(5, 'Motor Code: {}\n'.format(motor_code))

        pdf.write(5, 'Motor: {}{}\n'.format(motor_code, int(round(average_trust))))

        for column in data.columns.values:
            fig = data.plot(
                y=column,
                title='{} vs time (second)'.format(column),
                grid=True,
                use_index=True
            ).get_figure()

            fig.savefig('results/' + column.split(' ')[0])

            plt.close(fig=fig)

            pdf.add_page('L')
            pdf.image('results/' + column.split(' ')[0] + '.png')

        pdf.output('results/results.pdf', 'F')
      

if __name__ == "__main__":
    simulate()
