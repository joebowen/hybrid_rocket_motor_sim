from thermo.mixture import Mixture


class Thermochemical:
    def __init__(self):
        self.mixture = Mixture(['n2', 'h2o', 'co2'], zs=[.599, .224, .178], T=3000, P=3447000)

    def get_mixture(self):
        return self.mixture