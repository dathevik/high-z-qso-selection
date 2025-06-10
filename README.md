# Selection of high redshift quasars (z > 4.5)

This repository is for the 4MOST/CHANGES project, focusing on high-redshift quasar selection using SED fitting and color-based cuts. The workflow is modular, with each major step in a dedicated folder.

---

## Algorithmic Workflow

```
Input Data (~420 million)
        |
        v
+-------------------------------+
| Removing Initial Contaminants | (removing_initial_contaminants/)
+-------------------------------+
        |
        v
+-------------------+
|   SED Fitting     |   (sed_script/)
|   (~13 million)   |
+-------------------+
        |
        v
+-------------------+
|    Color Cut      |   (color_cut/)
+-------------------+
        |
        v
+-------------------+
| Statistical Cuts  |
|   (~10 million)   | 
+-------------------+
        |
        v
+-------------------+
|  Prioritization   |
|   (~930,000)      |
+-------------------+
        |
        v
+-------------------+
|  Initial Output   |  
|    (~19,000)      |  
+-------------------+
        |
        v
+---------------------------------------------+
| Crossmatch with DECaLS DR10: g_decal Det.?  |
+---------------------------------------------+
        |                          |
      Yes                          No
        |                          |
        v                          v
+-------------------+    +-------------------+
|  SNR Cut          |    |  SNR Cut          |
|  (~17,000)        |    |  (~17,000)        |
+-------------------+    +-------------------+
        |                          |
        v                          v
+-------------------+    +-------------------+
|   SED Fitting     |    |   SED Fitting     |
|   (~3,200)        |    |   (~3,200)        |
+-------------------+    +-------------------+
        |                          |
        v                          v
+-------------------+    +-------------------+
|   Color Cuts      |    |   Color Cuts      |
+-------------------+    +-------------------+
        |                          |
        v                          v
+-------------------+    +-------------------+
| Statistical Cut   |    | Statistical Cut   |
|   (~2,600)        |    |   (~2,600)        |
+-------------------+    +-------------------+
        |                          |
        v                          v
+-------------------------------+ +-------------------------------+
| Sample with DECaLS DR10       | | Sample with DELVE DR2         |
| photometry (1840 sources)     | | photometry (1314 sources)     |
+-------------------------------+ +-------------------------------+
        |                          |
        +-----------+--------------+
                    |
                    v
              +-------------+
              |   Merging   |
              +-------------+
                    |
                    v
         +------------------------+
         | Final Catalog: 3154    | 
         | sources                |      
         +------------------------+
```

---

## Folder Structure & Purposes

- **removing_initial_contaminants/**  
  Contains scripts and submodules for the first step: removing initial contaminants from the input data. This is the starting point of the workflow. It containts test results on each of the criteria in a separate folder and sed_cut.py which is a script of the initial cuts.

- **sed_script/**  
  Contains all SED fitting scripts and their dependencies. This is where the SED fitting step is performed after contaminants are removed and the magnitudes are converted to fluxes (mJy).  
  - `input_test/` and `output_test/` subfolders are used for input and output data, respectively, to keep the logic and file paths consistent.

  **Key scripts:**

  - `sed_calculation_command.py`:  
    This script performs SED fitting on input catalogs.  
    **How to run:**  
    ```
    python sed_calculation_command.py -i <input_file>
    ```
    where `<input_file>` is a file (e.g., `F23_fluxes_test.dat`) located in the `input_test/` subfolder.  
    **Test command:**  
    - Example: python sed_calculation_command.py -i F23_fluxes_test.dat (High-z qso samples from Fan+2023)
    
    **Outputs:**
    - The main output is a results file saved in the `output_test/` subfolder. The output file contains SED fitting results for each object in the input catalog, including best-fit parameters, chi-squared values, and statistical metrics for each template.


  - `sed_plot.py`:  
    This script generates plots from the SED fitting results and template files.  
    **How to run:**  
    ```
    python sed_plot.py
    ```
    The script will prompt for the object name, input catalog name, and list type (1 for QSO, 2 for BD, 0 for test).   

      **Test run:**  
   - Input the name of object J000009.99-041626.09 
   - Input the name of catalog F23_fluxes_test.dat
   - Input the type of list (1 for QSO, 2 for BD, 0 for test): 1

    **Outputs:**
    - The script generates plots visualizing the SED fitting results and template comparisons. These plots are saved as image files (e.g., PNG) in the `output_test/` subfolder.

- **color_cut/**  
  Contains scripts for applying color-based cuts to the data after SED fitting. This step further refines the candidate selection.

- **paper/**  
  Contains scripts, plots, and results used in publications and documentation of the project (Mkrtchyan et al. 2025 in prep.) Provided only with access.

- **query/**  
  Contains scripts for querying external catalogs (e.g., Gaia, DELVE, DECaLS) to crossmatch and supplement the main dataset.

- **output_catalogs/**  
  Stores the final and intermediate output catalogs generated at various steps of the workflow. Provided only with access.

- **tests/**  
  Contains test scripts and data for validating the pipeline and its components.

- **qso_bd_samples/**  
  Contains reference samples of known quasars and brown dwarfs for validation and comparison.Provided only with access.

---

## Notes

- The workflow is modular; you can run each step independently.
- The schematic above reflects the logical flow and data reduction at each stage.
- For more details on each script, see the comments and docstrings within the code.

