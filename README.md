# Table of contents

# Short introductory guide
`radexi-cmd` is an extension of RADEX computer program, which is used for atomic and molecular lines strength calculation. You can read about the original program on it's [official page](https://personal.sron.nl/~vdtak/radex/index.shtml) and you may even try to use it's [online variant](http://var.sron.nl/radex/radex.php). I also provided some information about it below, but fell free to correct me if im wrong.

Remember that when you use either official RADEX or radexi-cmd -- always refer to this [article](https://ui.adsabs.harvard.edu/abs/2007A%26A...468..627V/abstract).

## Download and installation

#### Dependencies
In order to download and install `radexi-cmd` you'll need
- gcc
- make
- gsl
- git (optional)
- curl (optional)
 
Install them with your packet manager.

#### Proceed with installation
It's a good practice to use your `$(HOME)/build` directory for downloading and building packages, so let's create one:

```bash
$ mkdir $HOME/build
$ cd $HOME/build
```

Now let's clone the repository with `git`, but you should choose the right branch for it. I try to keep `main` branch stable, but it's recommended to choose one of the `release/*` branches (easier to cope with bugs if you'll report them to me):

``` bash
$ git clone https://github.com/doedre/radexi-cmd.git -b release/*
$ cd radexi-cmd
``` 

If you don't have git you may download any of the release archives from this link. Just download it and unpack in any directory you want. 
Finally we can build this package either with the specified `build.sh` script or manually. If you want to use script, enter the program's directory and run it

``` bash
$ cd radexi-cmd
$ ./build.sh
``` 
#### Result of the installation
New folder called `radexi` should appear in your `$(HOME)` directory. Here the program stores generated databases for molecules and settings. The executable called `radexi` is also stored here.

In order to use the program freely you should create your local `bin` directory (or use `/usr/local/bin` or anything if you know how it works) and store the executable here

```bash
$ mkdir $HOME/bin
$ cp $HOME/radexi/radexi $HOME/bin/
$ echo "export PATH=$PATH:$HOME/bin" >> $HOME/.bashrc
```

Now you'll be able to execute `radexi` from terminal. Check it with

```bash
$ radexi --version
```

If you see no errors then all was installed correctly.

## Usage
Here I will describe the most common usage examples. Text enclosed in < > should be changed with your parameters. If you want to know more features follow the full usage guide (~ minutes).

##### Running the dialogue 

```bash
$ radexi
```

##### Adding molecule to the database
In order to use any molecule from LAMDA database you can add it to `radexi`'s local database

```bash
$ radexi --add-molecule <name> <path to the database file>
```

So if you want to add methanol molecule and call it *my_favourite_one*

```bash
$ radexi --add-molecule my_favourite_one <path to the file>/ch3oh.dat
```

##### Using input files
The `-r` (or `--result`) flag specifies the location of the output file.

```bash
$ radexi -r <path>/result.txt
```

If no path given, results will be stored in the current directory in `radexi_output.txt`. The same file is created if you've only specified the folder. 
---
# Full guide
Will appear 
## Notes on the installation
## Usage
### Starting guide
## When you should use RADEX

# For interested in this program
### How to report a bug
### How to ask for adding/suggesting a new feature
### How to contribute
## My thanks
## Plans for future releases

# References
