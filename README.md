# PSYCHO - Pure Data Library

### Version 0.0.1 alpha-0 (Unreleased)

--------------------------------------------------------------------------

This project is still in a very experimental/initial phase.

This is a library of Pd externals that includes some objects that deals with psychoacoustic measures and conversions. This was part of my masters and phd research, when I first started coding and now I'm finally revisiting it and getting it more presentable and out there (this is still quite far). 

The work I'm bringing back to life here had been presented as a "Psychoacoustic/Dissonance Toolbox for Pure Data", but now I'm just calling it the "PSYCHO" library!

--------------------------------------------------------------------------

### Acknowldegdements

This work started in my master's research in 2005. I need to thank people that got me started with coding in Pure Data at that time at NICS/Unicamp: Jônatas Manzolli, Fábio Furlanete. Later, in my PhD, André Pires helped me write a roughness external and we wrote a paper together. Marcelo Queiroz, my PhD supervisor from the computer department a USP was also very helpful. During my research internship at CIRMMT/McGill in 2010, I took coding lessons with Mathieu Bouchard who helped me code many of these.

--------------------------------------------------------------------------

####Building PSYCHO for Pd Vanilla:

PSYCHO relies on the build system called "pd-lib-builder" by Katja Vetter (check the project in: <https://github.com/pure-data/pd-lib-builder>). PdLibBuilder tries to find the Pd source directory at several common locations, but when this fails, you have to specify the path yourself using the pdincludepath variable. Example:

<pre>make pdincludepath=~/pd-0.49-0/src/  (for Windows/MinGW add 'pdbinpath=~/pd-0.49-0/bin/)</pre>

* Installing with pdlibbuilder

go to the pd-psycho folder and use "objectsdir" to set a relative path for your build, something like:

<pre>make install objectsdir=../psycho-build</pre>

Then move it to your preferred install folder for Pd and add it to the path.

Cross compiling is also possible with something like this

<pre>make CC=arm-linux-gnueabihf-gcc target.arch=arm7l install objectsdir=../</pre>

--------------------------------------------------------------------------

##Objects:

- [bark2hz]
- [barks~]
- [centroid~]
- [curve.gen]
- [db2phon]
- [dec2frac]
- [flunson]
- [harmonicity]
- [hz2bark]
- [hz2mel]
- [indigestibility]
- [iso226]
- [mel2hz]
- [phon2db]
- [phon2sone]
- [roughcurve]
- [roughness]
- [sharpness~]
- [sone2phon]
