"""Read in species and their fluxes."""

import yaml

yaml_flux = open("flux.yaml", 'r')
# temp = yaml.load(yaml_flux, Loader=yaml.FullLoader)
temp = yaml.safe_load(yaml_flux)

sp_name = temp['Species']
sp_flux = temp['Flux']
sp_flux.update((k, float(v)) for k, v in sp_flux.items())

mat_name = temp['Materials']

from Species import sp_full_list

sp_run_list = {k:sp_full_list[k] for k in sp_name}
# print(species)
# print(flux)
# print(sp_run_list)


sp_weight = [sp_flux[k] for k in sp_name]

if __name__ == '__main__':
    from random import choices
    sp_temp = choices(sp_name, weights=sp_weight, k=1)
    print(*sp_temp)
    print(sp_name)
    

