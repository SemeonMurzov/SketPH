import numpy as np


# ideal gasdynamics
nodetype          = "nodetype"
material          = "material"
coords            = "coords"
velocity          = "velocity"
density           = "density"
size              = "size"
mass              = "mass"
energy            = "energy"
pressure          = "pressure"
force             = "force"
soundSpeed        = "soundSpeed"
energyRate        = "energyRate"
volumeRate        = "volumeRate"
neibsCell         = "neibsCell"
densityRate       = "angularVelocity"
KinEnergyPair     = "KinEnergyPair"

# for elastoplaticity
deviatorStress    = "deviatorStress"
deviatorSRate     = "deviatorStrainRate"
soundSpeedL       = "soundSpeedL"
soundSpeedT       = "soundSpeedT"
angularRate       = "angularVelocity"
shearModulus      = "shearModulus"
elpli             = "elpli"

ndim=2
fields_dtype = {
	coords:         np.dtype((str(ndim)+'f8')),
	velocity:       np.dtype((str(ndim)+'f8')),
	force:          np.dtype((str(ndim)+'f8')),
	density:        np.float64,
	energy:         np.float64,
	pressure:       np.float64,
	deviatorStress: np.dtype((str(ndim*ndim)+'f8')),
	deviatorSRate:  np.dtype((str(ndim*ndim)+'f8')),
	size:           np.float64,
	nodetype:       "20S",
	material:       np.int64,
	mass:           np.float64,
	soundSpeed:     np.float64,
	soundSpeedL:    np.float64,
	soundSpeedT:    np.float64,
	energyRate:     np.float64,
	volumeRate:     np.float64,
	angularRate:    np.float64,
	KinEnergyPair:  np.float64,
	shearModulus:	np.float64,
	neibsCell:      np.int64,
	elpli:          bool,
}

fields_list = [
	coords,
	velocity,
	force,
	density,
	energy,
	pressure,
	deviatorStress,
	deviatorSRate,
	size,
	nodetype,
	material,
	mass,
	soundSpeed,
	soundSpeedL,
	soundSpeedT,
	energyRate,
	volumeRate,
	angularRate,
	KinEnergyPair,
	shearModulus,
	elpli
]
# sub list of that list which are
fields_list_in_out = [
	coords,
	velocity,
	force,
	density,
	energy,
	pressure,
	size,
	deviatorStress,
	soundSpeed,
	soundSpeedL,
	soundSpeedT
]