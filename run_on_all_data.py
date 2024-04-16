#%%
import os
import os.path as ospath
import pandas as pd
from standard_pipeline import pipeline_single
import numpy as np
from lib.identification.identification_pipelines.identification_processing import get_best_species_from_xml
#%%
expedition_folder = "expedition_jardin_botanique"

expedition_folder_path = ospath.join("data", expedition_folder)
init_db = pd.read_csv(ospath.join(expedition_folder_path, "sample_info.csv"))
df = init_db[["Barcode", "Species", "gene"]]
#df.insert(3,"test", np.zeros(27))
print(df)
#%%
for barcode in df["Barcode"]:
    if barcode != 70 and barcode != 72 and barcode != 71 and barcode != 75 and barcode != 92 and barcode != 77 and barcode != 83:
        input_expedition_path= ospath.join("data", expedition_folder)

        for root, dirs, files in os.walk(input_expedition_path):
            if root.endswith(str(f"barcode{barcode}")):
                input_folder_path = root

        for _,_, files in os.walk(input_folder_path):
            for file in files:
                if file.endswith(".fastq"):
                    input_fastq_path = ospath.join(input_folder_path, file)
                    input_fastq_filename=file
                
            base_name = ospath.splitext(input_fastq_filename)[0]
        
        pipeline_single(f"barcode{barcode}", input_fastq_path= input_fastq_path, db=None, windows=False, expedition_name="expedition_jardin_botanique")

#%%
df.insert(3, "matK", [None]*27)
df.insert(4, "rbcL", [None]*27)
df.insert(5, "psbA-trnH", [None]*27)
df.insert(6, "ITS", [None]*27)
print(df)
# %%
for barcode in df["Barcode"]:
    result_path = ospath.join("output",expedition_folder,f"barcode{barcode}")
    if ospath.exists(result_path):
        identification_results_path = ospath.join(result_path,"identification")
        for gene in ["matK", "psbA-trnH", "rbcL", "ITS"]:
            result_tuple = get_best_species_from_xml(ospath.join(identification_results_path, f"{gene}.txt"))
            df.loc[df["Barcode"]==barcode, gene] = f"species: {result_tuple[0]}, with score {result_tuple[1][0]} and evalue {result_tuple[1][1]})"
# %%
df.to_excel("results_pipeline_jardin_botanique.xlsx")
# %%
