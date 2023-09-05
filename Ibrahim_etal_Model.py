#This script is used for the model in the manuscript 'From Dome to Duplex: Convergent Gravitational Collapse Explains Intracratonic Dome and Nappe Tectonics in  Central Australia"
#Authors: Ibrahim, Y.; Rey, P.F.; Whitney, D.L.; Teyssier, C.; Cenki, B.; Roger, F.
#This script was written by Youseph Ibrahim and Patrice F. Rey. 

#We use Underworld version 2.13

import underworld as uw
from underworld import UWGeodynamics as GEO
import underworld.function as fn
import numpy as np
import math


# Scaling
u = GEO.UnitRegistry
resolution = (768,384) # 500m resolution
half_rate = 0.5 * u.millimeter / u.year
model_length = 384e3 * u.meter
surfaceTemp = 293.15 * u.degK
baseModelTemp = 1673.15 * u.degK
bodyforce = 3150 * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2

KL = model_length
Kt = KL / half_rate
KM = bodyforce * KL**2 * Kt**2
KT = (baseModelTemp - surfaceTemp)

GEO.scaling_coefficients["[length]"] = KL
GEO.scaling_coefficients["[time]"] = Kt
GEO.scaling_coefficients["[mass]"]= KM
GEO.scaling_coefficients["[temperature]"] = KT

# rcParams


GEO.rcParams["CFL"] = 0.2
GEO.rcParams["shear.heating"] = True
GEO.rcParams["popcontrol.split.threshold"] = 0.95
GEO.rcParams["popcontrol.max.splits"] = 100
GEO.rcParams["swarm.particles.per.cell.2D"] = 60
GEO.rcParams["popcontrol.particles.per.cell.2D"] = 60
GEO.rcParams["advection.diffusion.method"] = "SLCN"

# Model

Model = GEO.Model(elementRes=resolution, 
                  minCoord=(0. * u.kilometer, -168. * u.kilometer), 
                  maxCoord=(384. * u.kilometer, 24. * u.kilometer), 
                  gravity=(0.0, -9.81 * u.meter / u.second**2))

Model.outputDir="CA_200"

Model.diffusivity = 9e-7 * u.metre**2 / u.second 
Model.capacity    = 1000. * u.joule / (u.kelvin * u.kilogram)


# Radiogenic Heat Production

IrindinaRHP = 2 * 1.67e-6 * u.watt / u.meter**3
UpperBradyGneissRHP = 2 * 2.34e-6 * u.watt / u.meter**3
LowerBradyGneissRHP = 2 * 2.51e-6 * u.watt / u.meter**3
StavonsGneissRHP = 2 * 1.67e-6 * u.watt / u.meter**3 
SedimentRHP = 2 * 0.88e-6 * u.watt / u.meter**3 

# Materials 

air = Model.add_material(name="Air", shape=GEO.shapes.Layer(top=Model.top, bottom=0. * u.kilometer))
air.density = 1. * u.kilogram / u.metre**3
air.diffusivity = 1e-5 * u.metre**2 / u.second
air.capacity = 100. * u.joule / (u.kelvin * u.kilogram)
air.compressibility = 0.  # Compressibility should be zero when using Lecode isostasy

# Sediments

Sediment1 = Model.add_material(name="Sediment1", shape=GEO.shapes.Layer(top=0. * u.kilometer, bottom=-2. * u.kilometer))
Sediment1.radiogenicHeatProd =SedimentRHP
Sediment1.density  = GEO.LinearDensity(reference_density  = 2400. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)
Sediment2 = Model.add_material(name="Sediment2", shape=GEO.shapes.Layer(top=-2. * u.kilometer, bottom=-4. * u.kilometer))
Sediment2.radiogenicHeatProd =SedimentRHP
Sediment2.density  = GEO.LinearDensity(reference_density  = 2420. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)
Sediment3 = Model.add_material(name="Sediment3", shape=GEO.shapes.Layer(top=-4. * u.kilometer, bottom=-6. * u.kilometer))
Sediment3.radiogenicHeatProd =SedimentRHP
Sediment3.density  = GEO.LinearDensity(reference_density  = 2440. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)
Sediment4 = Model.add_material(name="Sediment4", shape=GEO.shapes.Layer(top=-6. * u.kilometer, bottom=-8. * u.kilometer))
Sediment4.radiogenicHeatProd =SedimentRHP
Sediment4.density  = GEO.LinearDensity(reference_density  = 2460. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

Sediment5 = Model.add_material(name="Sediment5", shape=GEO.shapes.Layer(top=-8. * u.kilometer, bottom=-10. * u.kilometer))
Sediment5.radiogenicHeatProd =SedimentRHP
Sediment5.density  = GEO.LinearDensity(reference_density  = 2480. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

Sediment6 = Model.add_material(name="BradyGneissSediment6", shape=GEO.shapes.Layer(top=-10. * u.kilometer, bottom=-12. * u.kilometer))
Sediment6.radiogenicHeatProd = UpperBradyGneissRHP
Sediment6.density  = GEO.LinearDensity(reference_density  = 2500. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

Sediment7 = Model.add_material(name="BradyGneissSediment7", shape=GEO.shapes.Layer(top=-12. * u.kilometer, bottom=-14. * u.kilometer))
Sediment7.radiogenicHeatProd = LowerBradyGneissRHP
Sediment7.density  = GEO.LinearDensity(reference_density  = 2520. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

Sediment8 = Model.add_material(name="IrindinaGneiss8", shape=GEO.shapes.Layer(top=-14. * u.kilometer, bottom=-16. * u.kilometer))
Sediment8.radiogenicHeatProd = IrindinaRHP
Sediment8.density  = GEO.LinearDensity(reference_density  = 2540. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

Sediment9 = Model.add_material(name="IrindinaGneiss9", shape=GEO.shapes.Layer(top=-16. * u.kilometer, bottom=-18. * u.kilometer))
Sediment9.radiogenicHeatProd = IrindinaRHP
Sediment9.density  = GEO.LinearDensity(reference_density  = 2560. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

Sediment10 = Model.add_material(name="SIrindinaGneiss10", shape=GEO.shapes.Layer(top=-18. * u.kilometer, bottom=-20. * u.kilometer))
Sediment10.radiogenicHeatProd = IrindinaRHP
Sediment10.density  = GEO.LinearDensity(reference_density  = 2580. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

Sediment11 = Model.add_material(name="IrindinaGneiss11", shape=GEO.shapes.Layer(top=-20. * u.kilometer, bottom=-22. * u.kilometer))
Sediment11.radiogenicHeatProd = IrindinaRHP
Sediment11.density  = GEO.LinearDensity(reference_density  = 2600. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

Sediment12 = Model.add_material(name="IrindinaGneiss12", shape=GEO.shapes.Layer(top=-22. * u.kilometer, bottom=-24. * u.kilometer))
Sediment12.radiogenicHeatProd = IrindinaRHP
Sediment12.density  = GEO.LinearDensity(reference_density  = 2620. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)
Sediment13 = Model.add_material(name="IrindinaGneiss13", shape=GEO.shapes.Layer(top=-24. * u.kilometer, bottom=-26. * u.kilometer))
Sediment13.radiogenicHeatProd = IrindinaRHP
Sediment13.density  = GEO.LinearDensity(reference_density  = 2640. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)
Sediment14 = Model.add_material(name="Stanovos_Gneiss14", shape=GEO.shapes.Layer(top=-26. * u.kilometer, bottom=-28. * u.kilometer))
Sediment14.radiogenicHeatProd = StavonsGneissRHP
Sediment14.density  = GEO.LinearDensity(reference_density  = 2660. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)
Sediment15 = Model.add_material(name="Stanovos_Gneiss15", shape=GEO.shapes.Layer(top=-28. * u.kilometer, bottom=-30. * u.kilometer))
Sediment15.radiogenicHeatProd = StavonsGneissRHP
Sediment15.density  = GEO.LinearDensity(reference_density  = 2680. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)
Sediment16 = Model.add_material(name="Stanovos_Gneiss16", shape=GEO.shapes.Layer(top=-30. * u.kilometer, bottom=-32. * u.kilometer))
Sediment16.radiogenicHeatProd = StavonsGneissRHP
Sediment16.density  = GEO.LinearDensity(reference_density  = 2700. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)
Sediment17 = Model.add_material(name="Stanovos_Gneiss17", shape=GEO.shapes.Layer(top=-32. * u.kilometer, bottom=-34. * u.kilometer))
Sediment17.radiogenicHeatProd = StavonsGneissRHP
Sediment17.density  = GEO.LinearDensity(reference_density  = 2720. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

# Basement Shape, overprints the sediments defined above

BasementShape = [(0. * u.kilometer, -6. * u.kilometer), (140 * u.kilometer, -6 * u.kilometer),
                (161. * u.kilometer, -12. * u.kilometer), (167 * u.kilometer, -12 * u.kilometer),
                (173. * u.kilometer, -18. * u.kilometer), (176 * u.kilometer, -18 * u.kilometer),
                (197 * u.kilometer, -34 * u.kilometer), (204 * u.kilometer, -34 * u.kilometer),
                (217 * u.kilometer, -2 * u.kilometer), (384 * u.kilometer, -2 * u.kilometer),
                (384 * u.kilometer, -40. * u.kilometer),(0 * u.kilometer, -40. * u.kilometer)]

Basement = Model.add_material(name="Basement", shape=GEO.shapes.Polygon(BasementShape))
Basement.radiogenicHeatProd = 1.5 * 0.88e-6 * u.watt / u.meter**3
Basement.density  = GEO.LinearDensity(reference_density  = 2820. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)


# Mantle
mantle = Model.add_material(name="Upper Mantle", shape=GEO.shapes.Layer(top=-40. * u.kilometer, bottom=-120. * u.kilometer))
mantle.density = GEO.LinearDensity(reference_density=3370. * u.kilogram / u.metre**3,
                                        thermalExpansivity= 2.8e-5 * u.kelvin**-1)

asthenosphere = Model.add_material(name="Asthenosphere", shape=GEO.shapes.Layer(top=-120. * u.kilometer, bottom=-168. * u.kilometer))
asthenosphere.density = GEO.LinearDensity(reference_density=3370. * u.kilogram / u.metre**3, thermalExpansivity= 2.8e-5 * u.kelvin**-1)

#Faults

fault = GEO.shapes.Polygon([(140 * u.kilometer, -6 * u.kilometer), (161. * u.kilometer, -12. * u.kilometer), (167 * u.kilometer, -12 * u.kilometer),
                (173. * u.kilometer, -18. * u.kilometer), (176 * u.kilometer, -18 * u.kilometer),
                (197 * u.kilometer, -34 * u.kilometer), (204 * u.kilometer, -34 * u.kilometer),
                (204 * u.kilometer, -32 * u.kilometer), (199 * u.kilometer, -32 * u.kilometer),
                (178. * u.kilometer, -16. * u.kilometer),(176 * u.kilometer, -16. * u.kilometer),
                (169. * u.kilometer, -10. * u.kilometer),(163 * u.kilometer, -10. * u.kilometer),(141. * u.kilometer, -5. * u.kilometer)])


# # Rheology

rh = GEO.ViscousCreepRegistry()
Model.minViscosity = 1e18 * u.pascal * u.second
Model.maxViscosity = 1e23 * u.pascal * u.second

air.viscosity = 1e18 * u.pascal * u.second

#Sediments

# Viscosity for sediments
sediment_stronger_viscosity = 0.1 * rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
sediment_weaker_viscosity = 0.01 * rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
decollement_viscosity = 1e19 * u.pascal * u.second

#Plasticity for sediments
sediment_plasticity = GEO.DruckerPrager(name="Sediment_Plasticity",
                                                cohesion=5. * u.megapascal,
                                                cohesionAfterSoftening=1. * u.megapascal,
                                                frictionCoefficient=0.54,
                                                frictionAfterSoftening=0.011,
                                                epsilon1=0.0, epsilon2=0.15)

sediment_layers = [Sediment5, Sediment6, Sediment7, Sediment9, Sediment10, Sediment11, Sediment13, Sediment14, Sediment15, Sediment17]
weaker_sediment_layers =[Sediment1, Sediment2, Sediment3]
decollement_sediments = [Sediment4,Sediment8,Sediment12,Sediment16]
deep_sediment_layers = [Sediment13, Sediment14, Sediment15, Sediment16, Sediment17]

for layer in weaker_sediment_layers:
    layer.viscosity = sediment_weaker_viscosity

for layer in sediment_layers:
    layer.viscosity = sediment_stronger_viscosity

for layer in decollement_sediments:
    layer.viscosity = decollement_viscosity

for layer_set in [weaker_sediment_layers, sediment_layers, decollement_sediments]:
     for layer in layer_set: 
        layer.plasticity = sediment_plasticity

# Basement

Basement.viscosity = 1 * rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
Basement.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=15. * u.megapascal,
                                                cohesionAfterSoftening=1.5 * u.megapascal,
                                                frictionCoefficient=0.54,
                                                frictionAfterSoftening=0.011,
                                                epsilon1=0.0, epsilon2=0.15)


# Mantle

mantle.viscosity = rh.Dry_Olivine_Dislocation_Karato_and_Wu_1993

asthenosphere.viscosity = 1e20 * u.pascal * u.second

# Faults

fault.viscosity = 1e19 * u.pascal * u.second
fault.plasticity = GEO.DruckerPrager(name="Continental Crust",
                                                cohesion=5. * u.megapascal,
                                                cohesionAfterSoftening=0.5 * u.megapascal,
                                                frictionCoefficient=0.54,
                                                frictionAfterSoftening=0.011,
                                                epsilon1=0.0, epsilon2=0.1)


# Melt Properties

solidii = GEO.SolidusRegistry()
my_crust_solidus = GEO.Solidus(A1=923 * u.kelvin, A2=-1.2e-07 * u.kelvin / u.pascal, A3=1.2e-16 * u.kelvin / u.pascal**2, A4=0.0 * u.kelvin / u.pascal**3)
mid_crust_solidus = GEO.Solidus(A1=1263 * u.kelvin, A2=-1.2e-07 * u.kelvin / u.pascal, A3=1.2e-16 * u.kelvin / u.pascal**2, A4=0.0 * u.kelvin / u.pascal**3)
mantle_solidus = solidii.Mantle_Solidus

liquidii = GEO.LiquidusRegistry()
my_crust_liquidus = GEO.Liquidus(A1=1423 * u.kelvin, A2=-1.2e-07 * u.kelvin / u.pascal, A3=1.6e-16 * u.kelvin / u.pascal**2, A4=0.0 * u.kelvin / u.pascal**3)
mid_crust_liquidus = GEO.Liquidus(A1=1763 * u.kelvin, A2=-1.2e-07 * u.kelvin / u.pascal, A3=1.6e-16 * u.kelvin / u.pascal**2, A4=0.0 * u.kelvin / u.pascal**3)
mantle_liquidus = liquidii.Mantle_Liquidus

for layer in deep_sediment_layers:
    layer.add_melt_modifier(my_crust_solidus, my_crust_liquidus,
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.2,
                         meltExpansion=0.13,
                         viscosityChangeX1 = 0.05,
                         viscosityChangeX2 = 0.20,
                         viscosityChange = 1e-3
                         )


Basement.add_melt_modifier(my_crust_solidus, my_crust_liquidus,
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.2,
                         meltExpansion=0.13,
                         viscosityChangeX1 = 0.05,
                         viscosityChangeX2 = 0.20,
                         viscosityChange = 1e-3
)


mantle.add_melt_modifier(mantle_solidus, mantle_liquidus,
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.03,
                         meltExpansion=0.13,
                         viscosityChangeX1 = 0.001,
                         viscosityChangeX2 = 0.03,
                         viscosityChange = 1e-2
                        )

mantle.phase_changes = GEO.PhaseChange((Model.y < GEO.nd(-120. * u.kilometer)), asthenosphere.index)

# Temperature Setup. We create a non-sustained thermal anomaly that is left to decay through time

# For the thermal anomaly
thermal_anomaly_shape = GEO.shapes.Disk((162 * u.kilometer, -140 * u.kilometer), (75. * u.kilometer))
Model.set_temperatureBCs(top=293.15 * u.degK, bottom = 1273.15 * u.degK,
                         nodeSets = [(air.shape, 293.15 * u.degK), (thermal_anomaly_shape,  1633.15 * u.degK), (asthenosphere.shape,  1633.15 * u.kelvin)])

Model.init_model()

#Resets temperature conditions without the thermal anomaly and without a defined asthenosphere temperature
Model.set_temperatureBCs(top=293.15 * u.degK, bottom = 1673.15 * u.degK, nodeSets = [(air.shape, 293.15 * u.degK)])

#Velocity Boundary Conditions
velocity = 0.5 * u.mm / u.year

Model.set_velocityBCs(left=[-velocity, None],
                      right=[velocity, None],
                      top=[None, None],
                      bottom = GEO.LecodeIsostasy(reference_mat=asthenosphere, average=False),
                      order_wall_conditions=['left','right','top','bottom'])


#Tracers for surface and Moho, and Finite Strain Ellipsoids
npoints_surface = 1536
coords_surface = np.ndarray((npoints_surface, 2))
coords_surface[:, 0] = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), npoints_surface)
coords_surface[:, 1] = GEO.nd(0. * u.kilometre)
surface_tracers = Model.add_passive_tracers(name="Surface", vertices=coords_surface)

npoints_moho = 1536
coords_moho = np.ndarray((npoints_moho, 2))
coords_moho[:, 0] = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), npoints_moho)
coords_moho[:, 1] = GEO.nd(0. * u.kilometre) - GEO.nd(Basement.bottom)
moho_tracers = Model.add_passive_tracers(name="Moho", vertices=coords_moho)

coords_FSE_Crust = GEO.circles_grid(radius = 1.0 * u.kilometer,
                           minCoord=[Model.minCoord[0], mantle.top],
                           maxCoord=[Model.maxCoord[0], 0.*u.kilometer])

FSE_Crust = Model.add_passive_tracers(name="FSE_Crust", vertices=coords_FSE_Crust)

#Erosion

Model.surfaceProcesses = GEO.surfaceProcesses.ErosionThreshold(air=[air],threshold = 2.0 * u.kilometer, surfaceTracers=surface_tracers)

#Solver Options
solver = Model.solver

solver.set_inner_method("superludist")
solver.options.A11.use_previous_guess = True
solver.options.scr.use_previous_guess = True
solver.set_penalty(1.0e3) #1e6
solver.options.scr.ksp_type="cg"






Model.run_for(20000000 * u.years, checkpoint_interval=50000 * u.year, restartStep=-1, restartDir="CA_200")


