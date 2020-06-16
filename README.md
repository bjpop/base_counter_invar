[![travis](https://travis-ci.org/bjpop/base_counter.svg?branch=master)](https://travis-ci.org/bjpop/base_counter)

# Overview 

This program counts DNA bases in targeted sequencing data.

In the examples below, `$` indicates the command line prompt.

# Licence

This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/bjpop/base_counter/master/LICENSE).

# Installing

You can install base_counter directly from the source code or build and run it from within Docker container.

## Installing directly from source code

Clone this repository: 
```
$ git clone https://github.com/bjpop/base_counter
```

Move into the repository directory:
```
$ cd base_counter
```

Python 3 is required for this software.

Base_counter can be installed using `pip` in a variety of ways (`$` indicates the command line prompt):

1. Inside a virtual environment:
```
$ python3 -m venv base_counter_dev
$ source base_counter_dev/bin/activate
$ pip install -U /path/to/base_counter
```
2. Into the global package database for all users:
```
$ pip install -U /path/to/base_counter
```
3. Into the user package database (for the current user only):
```
$ pip install -U --user /path/to/base_counter
```


