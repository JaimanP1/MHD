# Workflow

Some general tips and tricks to running and visualizing a simulation

## MHD

For binary data

### I/O

Always ensure that the `constants.f90` file in both `MHD/` and `VAPOR` is correctly set, otherwise will overwrite data.

`mv` the `NAMELIST` file to the output directory. Ensure that its timesteps are correct prior to doing so.

### Grid Size

Ensure that `num_procs` and dimension sizes are compatible.

### Job Submission

Set up a `.out` and `.err` directory; the messages are quite helpful

## VAPOR

### I/O

Same as previous step.

### Debug

Use `main.f90` and `subroutines.f90` for outputting debug files

## Visualizing

### GPU partition

API requires GPU to plot images, must `salloc` to a GPU node

### BOV

Run the `file_creation.sh` script to generate the BOV format

### Environment

Activate conda using `/project/cstr/Jaiman/sp25/env2/`

### API

Run `attempt10.py` to interact with the API for images

### Output

Use a linux desktop with `image magick` to view the images.

Use the ffmpeg script to stitch images into a video

`scp` the video to local, using `hostname -I` for the IP address 

# Current

The project I am currently working on

## Continuity equation

* Added a diffusion term to the continuity equation and re-ran the simulation of the confined eruption.
* Using density cross section to determine if continuity is satisfied
	* Currently, density is just 0 everywhere
	* Likely a problem with VAPOR API and not data itself, given that all other parameters seem fine
