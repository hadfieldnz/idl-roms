# The IDL-ROMS Library

Mark Hadfield

## Synopsis

IDL-ROMS is a library of IDL code written largely by me at NIWA for processing input and output from
the [ROMS ocean model](https://www.myroms.org/) and its siblings. IDL-ROMS is published under the
[MIT Open Source License](http://www.opensource.org/licenses/mit-license.php).
It is hosted on
GitHub in project [hadfieldnz/idl-roms](https://github.com/hadfieldnz/idl-roms).

I have a number of other ROMS-related routines in my code at NIWA, and am disentangling them and moving them
to IDL-ROMS as time permits.

## Installation

### Installation Method 1: IDL Package Manager

If you have IDL 8.7.1 (due out in September 2018) or later, you can install IDL-ROMS with the IDL Package Manager, eg:

```
IDL> ipm, /INSTALL, 'https://github.com/hadfieldnz/idl-roms'
```

This will install a package named IDL-ROMS in the !PACKAGE_PATH directory, typically ${HOME}/.idl/idl/packages.
The relevant subdirectories will also be added to the [!PATH](https://www.harrisgeospatial.com/docs/Managing_IDL_Paths.html).

### Installation Method 2: Cloning the source

If you don't have the IDL Package Manager the recommended method for installing IDL-ROMS is to clone the repository, eg:

```
$ cd ${HOME}/IDL
$ git clone https://github.com/hadfieldnz/idl-roms.git
```

You will then need to add the three code subdirectories (roms, san and examples) to your [!PATH](https://www.harrisgeospatial.com/docs/Managing_IDL_Paths.html).
The simplest way to do this is to add an entry like the following to the IDL path preferences dialogue:

```
'+/home/hadfield/IDL/idl-roms'
```

### Installation Method 3: Downloading a Zip archive

GitHub also allows you to download a snapshot of the code [as a ZIP archive](https://github.com/hadfieldnz/idl-roms/archive/master.zip).
You then need to extract the code into a suitable directory and modify the !PATH as for Method 2.

## Dependencies

IDL-ROMS requires the IDL-Motley library, also on GitHub in project [hadfieldnz/idl-motley](https://github.com/hadfieldnz/idl-motley).
If you are installing with the IDL Package Manager, then this dependency will be handled automatically.

________________________________________
Mark Hadfield 2018-08-06

