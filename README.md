# ConnectomeDB
This repository contains the scripts for creating the database website for [ConnectomeDB](https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/). Follow the instructions below to get started with the project.

## Prerequisites
Before setting up the project, ensure you have the following installed:

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) 
- Git
- [jupyter lab](https://jupyter.org/install)

## Clone the Repository

To get started, clone the repository to your local machine:

via HTTPS:

```bash
git clone https://github.com/bioinfo-YCU/ConnectomeDB.git
```

or via ssh:

```bash
git clone git@github.com:bioinfo-YCU/ConnectomeDB.git
```
then change directory to "ConnectomeDB":

```bash
cd ConnectomeDB
```

## Setting Up the Conda Environment
The project requires specific dependencies, which can be installed using the requirements.txt file. Follow the steps below:

Create a new conda environment with the necessary dependencies:

```bash
conda env create --file environment.yml
```
Activate the conda environment:

```bash
conda activate quarto_llm_env
```
This will ensure that the correct versions of all required libraries are installed.

## Data Directory
If you are a collaborator working on this project, please request the "data" directory from Sakura to access the necessary credentials/data. Once you have access to the data, place it within the project directory.

```bash
project-directory/
  ├── data/                # Place your requested data here
  ├── environment.yml
  ├── [other project files]
```

## Running the Project
You can render and preview the website to see changes you've made.

```bash
quarto render
quarto preview
```

### Note: to use as an env in Jupyter, make sure you add it as a new kernel

```bash
conda activate quarto_llm_env
python -m ipykernel install --user --name=quarto_llm_env --display-name "quarto_llm_env)"
```


