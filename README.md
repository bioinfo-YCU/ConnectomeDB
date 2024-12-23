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

## Contributing
We welcome contributions to the project! If you would like to contribute, please follow these steps when making changes:
Assuming you have cloned the repository and have completed the set-up as instructed above,

```bash
git checkout -b your-branch-name
```
Make your changes (e.g., fix a bug, add a feature).

Commit your changes:

```bash
git add .
git commit -m "Description of your changes"
```

Push your changes to your branch:
```bash
git push origin your-branch-name
```

Create a pull request (PR) from your forked repository to the main repository. You can do this by navigating to the "Pull Requests" tab in the original repository and clicking "New Pull Request."

Be sure to describe your changes clearly in the PR description, so the maintainers can understand the purpose and scope of your changes.
