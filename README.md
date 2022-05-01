# radexi-cmd

I've made this program in order to improve usage experience of RADEX, which is more of a script than a real program.
With this you won't need to recompile RADEX every time you want to experiment with cloud geometries or restart the
calculations if you've entered something wrong. Some of the working features:
- Better dialogue system with history (use arrows) and autocompletion (use TAB);
- Choose escape probability calculation method as one of the paramaters;
- Add atoms and molecules with custom names;
- Calculate intensities for several molecules in one run;
- See if line overlap happen.

In future plans:
- Make stable library with some kind of API;
- Being able to use H2 in place of pH2 and oH2 only if the user specifies it (in RADEX it is done automatically and I don't allow it);
- Improving dialogues
- Statistical analysis for your hypothesis on line intensities
- ...

# Short introductory guide

`radexi-cmd` is an extension of RADEX computer program, which is used for atomic and molecular lines strength
calculation. You can read about the original program on it's [official page](https://personal.sron.nl/~vdtak/radex/index.shtml)
and you may even try to use it's [online variant](http://var.sron.nl/radex/radex.php).

Remember that when you use either official RADEX or radexi-cmd -- always refer to this [article](https://ui.adsabs.harvard.edu/abs/2007A%26A...468..627V/abstract).

## Download and installation

#### Dependencies
In order to download and install `radexi-cmd` you'll need
- `gcc`
- `make`
- `gsl`
- `git` (optional)

Install them with your packet manager.

#### Proceed with installation
It's a good practice to use your `$(HOME)/build` directory for downloading and building packages, so let's create one:

```bash
$ mkdir $HOME/build
$ cd $HOME/build
```

Now let's clone the repository with `git`, but you should choose the right branch for it. I try to keep `main` branch
stable, but it's recommended to choose one of the `release/*` branches (easier to cope with bugs if you'll report them to me):

``` bash
$ git clone https://github.com/doedre/radexi-cmd.git -b release/v0.2-alpha
$ cd radexi-cmd
```

If you don't have git you may download any of the release archives from this link. Just download it and unpack in any directory you want.
Finally we can build this package with `make`

``` bash
$ make
```

This will compile a program and place a binary into local `bin` directory. If you want to add global executable (be able
to start it by simple `radexi` command) instead enter

```bash
$ sudo make install
```

It will put executable into `/usr/local/bin` directory.

## Usage
Here I will describe the most common usage examples. Text enclosed in < > should be changed with your parameters.

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

##### Using output files
The `-r` (or `--result`) flag specifies the location of the output file.

```bash
$ radexi -r <path>/result.txt
```

If no path given, results will be stored in the current directory in `radexi_output.txt`. The same file is created if you've only specified the folder.

---
# Full guide
Will appear
## Plans for future releases
v0.2-alpha:
- Ability to read input files
- Download database files directly from LAMDA
- Improve dialogue system
- Add logs and correct output on internal errors
- Verbose console output
# References
