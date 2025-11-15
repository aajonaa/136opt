# ExplorExploi - Independent Experiment Folder

This folder contains a self-contained version of the exploration-exploitation balance experiment framework.

## Main Script
- `exp_CCMWOA1.m` - Main experimental script for running algorithm comparisons with CEC2017 benchmark functions

## Algorithm Implementations
- `LRRIME/LRRIME.m` - Learning + Reflection RIME algorithm
- `RMRIME/RIME.m` - Standard RIME algorithm

## Utility Functions
### subFunction/
- `Get_Functions.m` - Benchmark function definitions (F1-F146)
- `initialization.m` - Population initialization
- `MySampling.m` - Data sampling for convergence curves

### statics/
- `Orderhao.m` - Statistical ranking of algorithms
- `pValueToExcelhao.m` - P-value calculations for statistical comparison
- `FridTest3.m` - Friedman test (per function)
- `FridTest4.m` - Friedman test (overall)

### Root Level
- `Divergence_calculation.m` - Diversity metric calculation

## Benchmark Functions
- `cec05_func.m` - CEC2005 benchmark functions
- `cec13_func.mexw64` - CEC2013 benchmark functions (compiled)
- `cec14_func.mexw64` - CEC2014 benchmark functions (compiled)
- `cec17_func.mexw64` - CEC2017 benchmark functions (compiled)
- `cec19_func.mexw64` - CEC2019 benchmark functions (compiled)

## Data Folders
- `input_data_cec2017/` - Input data files for CEC2017 benchmarks (328 .txt files)
- `cec05/` - Data files for CEC2005 benchmarks (.mat files)
- `exp_result/` - Output directory for experimental results

## How to Run
1. Open MATLAB
2. Navigate to the ExplorExploi folder
3. Run the main script:
   ```matlab
   exp_CCMWOA1
   ```

## Current Configuration
- Algorithms compared: LRRIME vs RIME
- Benchmark suite: CEC2017 (F107, F109-F136)
- Dimensions: 30
- Population size: 30
- Function evaluations: 300,000 (dim × 5000 × 2)
- Independent runs: 30

## Output Files
Results are saved to `exp_result/[AlgorithmName]-[timestamp]/`:
- Excel file with convergence data, statistics, and rankings
- Convergence curve figures (.fig and .tif)
- Balance analysis figures (.fig and .tif)
- Diversity analysis figures (.fig and .tif)

## Notes
- This folder is self-contained and can be moved/copied independently
- All dependencies are included
- No additional MATLAB toolboxes required beyond standard functions
- The script uses parallel processing (parfor) for faster execution

