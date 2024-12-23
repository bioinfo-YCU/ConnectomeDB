# ConnectomeDB
This repository contains the scripts for creating the database website for ConnectomeDB. Follow the instructions below to get started with the project.

## Prerequisites
Before setting up the project, ensure you have the following installed:

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/individual)
- Git

## Clone the Repository

To get started, clone the repository to your local machine:

```bash
git clone https://github.com/bioinfo-YCU/ConnectomeDB.git
cd ConnectomeDB
```
Here’s a template for a README.md that provides instructions for pulling the repository, setting up a conda environment, and requesting a "data" directory from Sakura (if collaborators need it):

# Project Name

This repository contains [brief description of the project]. Follow the instructions below to get started with the project.

## Prerequisites

Before setting up the project, ensure you have the following installed:

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/individual)
- Git

## Clone the Repository

To get started, clone the repository to your local machine:

```bash
git clone https://github.com/yourusername/repository-name.git
cd repository-name
Setting Up the Conda Environment
The project requires specific dependencies, which can be installed using the requirements.txt file. Follow the steps below:
```

## Setting Up the Conda Environment
The project requires specific dependencies, which can be installed using the requirements.txt file. Follow the steps below:

Create a new conda environment with the necessary dependencies:

```bash
conda create --name quarto_env --file requirements.txt
```
Activate the conda environment:

```bash
conda activate quarto_env
```
This will ensure that the correct versions of all required libraries are installed.

## Data Directory
If you are a collaborator working on this project, please request the "data" directory from Sakura to access the necessary credentials/data. Once you have access to the data, place it within the project directory.

```bash
project-directory/
  ├── data/                # Place your requested data here
  ├── requirements.txt
  ├── [other project files]
```

## Running the Project
You can render and preview the website to see changes you've made.

```bash
quarto render
quarto preview
```