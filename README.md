# Reproducibility Package for Section 4.2
This repository provides the Python code and instance generation procedure used to reproduce the computational results reported in Section 4.2 of our study, including Tables 5 and 6. The repository also includes the generated Excel files that document the instance data and results used in the reported experiment.
## Scope
This package reproduces the experiment in Section 4.2. 
The four models solved in this experiment are:
- `LT-TMM`: the proposed integrated model, referred to as `Our model` in the paper;
- `TMM`: the model without lateral transshipment;
- `LT-ERT-T`: the model without multi-modal transportation matching;
- `NSU-ERT-T`: the model without lateral transshipment and multi-modal transportation matching.
  
All four models use the same disaster relief network, deterministic parameters, samples, ambiguity set parameters, and solution procedure. They differ only in whether the corresponding strategic component is enabled.
## Software and hardware
The reported results were obtained using:
- Python 3.8.8 (Anaconda base environment);
- Gurobi Optimizer 10.0.1;
- gurobipy 10.0.1;
- NumPy 1.19.5;
- pandas 1.2.4;
- openpyxl 3.0.7;
- AMD Ryzen 7 5800U processor (1.90 GHz);
- 16 GB RAM.
  
A valid Gurobi license is required to run the optimization models.
The script explicitly sets `OutputFlag=0` and `TimeLimit=300`. Other Gurobi algorithmic parameters retain their default values. Consequently, small differences in runtime, alternative optimal decisions, or displayed decimal values may occur on different hardware or solver configurations. With the stated NumPy seed, Python environment, package versions, and code version, the reported results should be reproducible within normal numerical tolerances.
## Problem dimensions
The computational instance contains:
- 3 relief centers;
- 10 affected areas;
- 3 time periods;
- 100 samples for each uncertain parameter.
  
The relief center to affected area distances and relief center to relief center distances are specified directly in the script. These values correspond to the distance data reported in Tables 3 and 4 of the paper.
## Parameter generation
The NumPy random seed is set to `42` before generating the instance data.
The simulated parameters are generated according to the following distributions and values:
- basic resource demand: $U[200,1000]$;
- emergency resource demand: $U[100,500]$;
- shortage cost: $U[10,50]$;
- road condition parameter: $U[1,10]$;
- truck transportation time: $U[D_{ij},5D_{ij}]$;
- risk level: $\alpha_{jt}=0.05$;
- truck capacity: `750` pieces;
- helicopter capacity: `500` pieces;
- available trucks in each period: `15`;
- available helicopters in each period: `6`;
- delay cost coefficient: `0.2` Yuan/(piece·minute);
- basic resource transshipment cost coefficient: `0.008` Yuan/(piece·km);
- emergency resource transshipment cost coefficient: `0.02` Yuan/(piece·km).

The transportation cost coefficients are calculated as:

$$
c_{ij}=0.02D_{ij}, \qquad c'_{ij}=0.1D_{ij}+1.5.
$$

The initial inventories and replenishment quantities are fixed as follows:
```text
Initial basic resource inventory:     [5000, 2000, 2500]
Initial emergency resource inventory: [2500, 1500, 1000]
Basic resource replenishment b:
[[0, 2000, 4000],
 [50, 300, 2000],
 [500, 1500, 1500]]
Emergency resource replenishment b':
[[0, 3000, 3000],
 [20, 20, 1000],
 [800, 1500, 1500]]
```
## Ambiguity set construction
For each uncertain parameter, the script calculates the following quantities from the 100 samples:
- sample mean;
- sample mean absolute deviation (MAD) around the sample mean;
- sample minimum as the lower support bound;
- sample maximum as the upper support bound.
  
These quantities define the marginal-information ambiguity set used in the distributionally robust optimization model.
The data processing function is retained exactly as used to generate the reported results. In Section 4.2, its perturbation standard deviation is set to zero. Therefore, the processed samples are numerically equal to the initially generated uniform samples. The function is retained because it preserves the random number consumption sequence used when generating the reported samples and results.
## Running the experiment
Install the required Python packages:
```bash
pip install numpy==1.19.5 pandas==1.2.4 openpyxl==3.0.7 gurobipy==10.0.1
```
Ensure that Gurobi 10.0.1 is installed and that a valid Gurobi license is available. Then run:
```bash
python experiment_4_2_en.py
```
The script generates the instance data, constructs the ambiguity set parameters, solves the four models, calculates the cost components and synergy benefits, and writes the results to Excel files in the current working directory.
## Output files
The script generates the following files:
- `experiment_1_results.xlsx`: objective values and cost components for the four models, together with the generated deterministic and DRO parameters;
- `road condition_1.xlsx`, `road condition_2.xlsx`, and `road condition_3.xlsx`: road condition samples;
- `transportation time_1.xlsx`, `transportation time_2.xlsx`, and `transportation time_3.xlsx`: transportation time samples;
- `emergency demand.xlsx`: emergency demand samples.
## Reported results
The script reports the following cost components for each model:
- total cost;
- direct transportation cost, including truck and helicopter transportation;
- lateral transshipment cost;
- basic resource shortage cost;
- emergency resource delay cost.
  
It also calculates the following quantities for the proposed model relative to each benchmark:
- lateral transshipment investment;
- effectiveness gain from shortage and delay cost reductions;
- efficiency gain from direct transportation cost reduction;
- net spatial capacity synergy benefit.
  
The resulting values correspond to Tables 5 and 6 of the paper, subject to normal rounding differences.
## Reproducibility notes
The script uses the same generated samples for all four models. The random seed is fixed for NumPy data generation, while Gurobi uses its default algorithmic parameters apart from the output and time limit settings described above.
The generated Excel files are archival records of the instance data and reported results; the current script does not read them as model inputs and instead generates the samples during execution.
## Citation
If you use this code or data, please cite the corresponding paper.
