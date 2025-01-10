# Predictor for ps concentrations and critical temperatures (WIP)

## How to Install

- Install conda if you haven't:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
conda init
source ~/.bashrc
```

- `git pull` or `git clone` this repo
- Create virtual environment and install packages with conda:

```
conda create -n harryml python=3.10 -y
conda activate harryml
pip install -r requirements.txt
```

## How to Run

Right now, the main ML pipeline is implemented in the `ml_pipeline.ipynb` notebook.
This notebook was built from refactoring the code inside `predictor_densities_and_ct.ipynb`,
`predictor_density_ray.ipynb` and `predictor_grid_search.ipynb`.
It takes about 40 minutes to run the whole notebook.

- Copy the data to `data/FINISHED_SIMULATIONS`
- `conda activate harryml`
- `jupyter notebook`
- Open and run `ml_pipeline.ipynb`

## MLflow Experiments

ML training runs are stored as MLflow experiments. Experiments are stored
and can be viewed locally - i.e. they are not shared between computers.

To check the experiments registered so far after running `ml_pipeline.ipynb`:

Run `mlflow server --backend-store-uri ./experiments/`, and open the url shown
(like `http://127.0.0.1:5000`).
