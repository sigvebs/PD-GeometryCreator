Gemoetry generator for Peridynamics
==============
For use with the Peridyn library for Peridynamics simulations.

External libraries:
- CImg
- Boost filesystem
- libconfig++
- armadillo

Example usage:
"./meshGenerator path_to_configuration_file"

Example configuration file
```
nParticles = 5000
threshold = 1000
multiplicationFactor = 50
savePath = "/save/path"
imgPath = "path/to/image.tif"
periodic_x = true
periodic_y = true
```
