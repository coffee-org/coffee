[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

Latest Version (off main branch): [![latesttag](https://img.shields.io/github/tag/coffee-org/coffee.svg)](https://github.com/coffee-org/coffee/tree/main)

| Branch    | Build   | Docker Deployment    | Travis-CI    | Activity   |
|-------------|-------------|-------------|-------------|-------------|
**main**|[![CMake badge](https://github.com/coffee-org/coffee/actions/workflows/cmake.yml/badge.svg?branch=main)](https://github.com/coffee-org/coffee/actions/workflows/cmake.yml)|[![CMake badge](https://github.com/coffee-org/coffee/actions/workflows/docker-image.yml/badge.svg?branch=main)](https://github.com/coffee-org/coffee/actions/workflows/docker-image.yml)|[![Build Status](https://www.travis-ci.com/coffee-org/coffee.svg?branch=main)](https://www.travis-ci.com/coffee-org/coffee)|![lastcommit](https://img.shields.io/github/last-commit/coffee-org/coffee/main.svg)|
**dev**|[![CMake badge](https://github.com/coffee-org/coffee/actions/workflows/cmake.yml/badge.svg?branch=dev)](https://github.com/coffee-org/coffee/actions/workflows/cmake.yml)|[![CMake badge](https://github.com/coffee-org/coffee/actions/workflows/docker-image.yml/badge.svg?branch=dev)](https://github.com/coffee-org/coffee/actions/workflows/docker-image.yml)|[![Build Status dev](https://www.travis-ci.com/coffee-org/coffee.svg?branch=dev)](https://www.travis-ci.com/coffee-org/coffee)|![lastcommit](https://img.shields.io/github/last-commit/coffee-org/coffee/dev.svg)|


Code metrics (dev branch) :
[![CodeScene Code Health](https://codescene.io/projects/14780/status-badges/code-health)](https://codescene.io/projects/14781)
[![CodeScene System Mastery](https://codescene.io/projects/14780/status-badges/system-mastery)](https://codescene.io/projects/14781)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/6eefa0e1c1254889b1e2f6fda55930ca)](https://www.codacy.com/gh/coffee-org/coffee/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=coffee-org/coffee&amp;utm_campaign=Badge_Grade)
[![CodeFactor](https://www.codefactor.io/repository/github/coffee-org/coffee/badge)](https://www.codefactor.io/repository/github/coffee-org/coffee)



# coffee
Coronagraph Optimization For Fast Exoplanets Exploration


## Installing coffee

Install pre-requisite packages as needed. For example, under Ubuntu :

	sudo apt install git make dpkg-dev libc6-dev cmake pkg-config python3-dev libcfitsio-dev pybind11-dev python3-pybind11 libgsl-dev libfftw3-dev libncurses-dev libbison-dev libfl-dev libreadline-dev pkg-config gcc-10 g++-10 


Coffee is a milk plugin. To install coffee, follow the milk installation steps, adding extra coffee plugins.
To compile and install executable in local directory :

    git clone https://github.com/milk-org/milk.git
    cd milk
    ./fetch_coffee_dev.sh
    ./compile.sh $PWD/local

## Run coffee

Executable is located in subdirectory of path specified as argument to compile.sh script:

    ./local/milk-1.01.02/bin/coffee

## Tutorial

Milk tutorial:

    ./local/milk-1.01.02/bin/milk-tutorial

Coffee tutorial to be done
