import os 
import os.path as ospath
from lib.general_helpers.generate_fastq import multiplex

fake_multiplex_path = ospath.join("data","fake_multiplexed")

folder_paths = []
expedition = "expedition_jardin_botanique"
species_to_multiplex= "Photinia_davidiana"
for root,dirs,files in os.walk(ospath.join("data",expedition)):
    for dir in dirs:
        if species_to_multiplex in dir:
            folder_paths.append(ospath.join(root,dir))
print(folder_paths)

dst= ospath.join(fake_multiplex_path, species_to_multiplex+"_mutiplexed.fastq")
multiplex(folder_paths,dst)