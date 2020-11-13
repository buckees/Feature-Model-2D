"""Read in species and their fluxes."""

import yaml

yaml_flux = open("flux.yaml", 'r')
temp = yaml.load(yaml_flux, Loader=yaml.FullLoader)
print(temp)

species = temp['Species']
flux = temp['Flux']

print(species)
print(flux)

from Species import sp_full_list

sp_run_list = {k:sp_full_list[k] for k in species}
print(sp_run_list)
    

