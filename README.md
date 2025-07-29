# PsychTimeCost

## About The Project
**PsychTimeCost repository contains MATLAB code for the following tasks**
<br />
&emsp; (1) Fit a non-parametric Drift-Diffusion Model (nDDM) to behavioral data
<br />
&emsp; (2) Find the optimal bound shape by Dynamic Programming with the fitted parameters from (1)

## To clone from GitHub
Open the Terminal and change the directory to where you want to clone the repo.
```sh
cd YOUR_FOLDER
```
Then execute this command to clone the repo.
```sh
git clone https://github.com/shinkira/PsychTimeCostDemo.git
```
## Prerequisites

Install `MATLAB Coder` `Parallel Computing Toolbox` `Statistics and Machine Learning Toolbox` in MATLAB.

MEX function should be used to run Dynamic Programming faster. Compile the MEX function with this command.
```sh
compile_d2b_dp
```

## Usage

Run demo code to fit nDDM to simulated data and compute the optimal bounds.
```sh
demo
```
Run on the actual data with the following MATLAB script.
```sh
run_pipeline
```
Please request the data from Shin Kira (skira@fsu.edu). The data should be saved in `SubjectData` folder under the main repository (PsychTimeCost).
<br />
Hint! &ensp; Set `debug = 1` to run everything faster (with lesss accuracy) and grasp the overall process.

## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

