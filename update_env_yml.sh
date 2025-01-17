#!/bin/bash
conda env export --name quarto_llm_env --no-builds | grep -v "^prefix:" > environment.yml
