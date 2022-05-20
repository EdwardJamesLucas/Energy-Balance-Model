def myfunc(x):
    print(x + 5)

class Media:
    def __init__(self, name):
        if name == "water":
            self.density = float(1000)
            self.conduction_coeff = float(0.6)
            self.heat_capacity = float(4200)
        else:
            print("Only water supported")

class Storage:
    def __init__(self, height, footprint, no_of_layers, time_step):
        self.height = height
        self.footprint = footprint
        self.no_of_layers = no_of_layers
        self.time_step = time_step

    def my_storage_func(self):
        print("Hi Eddie")

class Dog:
    def __init__(self, name):
        self.name = name
        
    def bark(self):
        print("bark bark")