# PSYCHO - Pure Data Library

### Version 1.0 alpha-1

--------------------------------------------------------------------------

About: 

This project is being released under the GNU GENERAL PUBLIC LICENSE, see LICENSE.

This is a library of Pd externals that includes some objects that deals with psychoacoustic measures and conversions. This was part of my masters and phd research. I'm now finally revisiting it and getting it more presentable and out there. 

The work I'm bringing back to life here had been presented as a "Psychoacoustic/Dissonance Toolbox for Pure Data", but now I'm just calling it the "PSYCHO" library! The work has been deeply revised and rewritten from scratch since last I touched it in 2012 (when I finished my PhD). I learned a lot in how to program externals since then, and Pd also changed a lot. Most of my patches were based on Pd-Extended, which is now long gone and abandoned. Hence, please cope with me and hang in there. 

There are just a few abstractions in the library. I have mostly compiled objects, but even so the code has also been  implemented as Pd Vanilla patches and presented as such in the help file!

This work is still in a very experimental/initial alpha phase. I'm still putting other objetcs and tools up until the final release. Sta tuned for tools that implements techniques such as the spectral mapping (change a sound's spectrum to match a given scale) and more.  

--------------------------------------------------------------------------

### Acknowldegdements

This work started in my master's research in 2005. I need to thank people that got me started with coding in Pure Data at that time at NICS/Unicamp: Jônatas Manzolli, Fábio Furlanete. Later, in my PhD, Marcelo Queiroz, my PhD supervisor from the computer department at USP was also very helpful. André Pires helped me write a roughness external and we wrote a paper together in 2009.  During my research internship at CIRMMT/McGill in 2010, I took coding lessons with Mathieu Bouchard who helped me get started with most of these. 

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

##Objects (25 in total):

- [al]
- [bark2hz]
- [barks~]
- [centroid~]
- [commonality]
- [curve.gen]
- [db2phon]
- [dec2frac]
- [distance]
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
- [salience]
- [sharpness~]
- [sone2phon]
- [tonalness]
- [yl]

--------------------------------------------------------------------------

### References for Main objetcs and Models:

- The main references I have for this work are:

-Porres, A. T. (2012). Modelos psicoacústicos de dissonância para eletrônica ao vivo. Tese de Doutorado, Escola de Comunicações e Artes, Universidade de São Paulo, São Paulo. 

-Porres, A.T. “A Dissonance Model for Live Electronics” Proceedings of the 4th International Conference of Students of Systematic Musicology (SysMus11), Cologne, Germany, 2011 

-Porres, A.T. “Dissonance Model Toolbox in Pure Data” Proceedings of the 4th Pure Data International Convention (PdCon11), Berlin, Germany, 2011  

- Roughness:

Most of my resaearch work was in the revision of roughness models. Roughness curves can be used to derive musical scales based on a spectrum and to measure dissonance (I'm yet to put up objects that can perform spectral changes to match a musical scale).

The main reference for my Roughness model is the work by Clarence Barlow.

- Indigestibility / Harmonicity:

These are two concepts developed by Clarence Barlow in: Barlow, Clarence & Lohner, Henning. (1988). Two Essays on Theory. Computer Music Journal. 11. 44. 10.2307/3680177. 

- Tonalness / Tone Salience  /Multiplicity / Pitch Commonality / Pitch Distance: 

Reference: Parncutt, Richard & Strasburger, Hans. (1994). Applying psychoacoustics in composition: Harmonic progressions of non-harmonic sonorities. Perspectives of New Music. 32. 88-129. 10.2307/833600. 

These are from a Pitch Model by Parncutt & Strasburger (1994), based on the theoretical work of Ernst Terhardt in Pitch & Consonance, specially Pure/Complex tone Sensation and Sonorousness (a measure of Pitch clarity). The main outputs of the model are: Tonalness (same as sonorousness), Pitch Commonality and Pitch Distance. The Pitch Commonality is a measure of the degree in which two sonorities evoke common pitches. The Pitch Distance is a similar concept that considers the probability of noticing a pitch from one sonority in another, but it's more pertinent in the melodic context. The model considers an input spectrum and is not generalized in octaves.

The input spectrum are lists of frequency values and linear amplitudes. The frequency list needs to first be converted to Pitch Category (which is basically a logarythmic MIDI scale). The amplitude list needs to first be converted to Auditory Levels Y(L) - defined as dB over the partial’s Level Threshold L(TH). The next stage converts to Audible Level AL(P), which is the Audible level after a masking stage. This is the used to calculate Pure Tone Audibilities Ap(P) also defined as the Spectral Pitch Weigth.

With this new input we can recognice a harmonic pattern. The Complex Tone Audibility Ac(P) (or Virtual Pitch Weight) is a degree to which a Harmonic Series or part of it is present in the spectrum. For that, a template of 10 harmonics is used, each with a different weigth. When there’s a match, the Complex Tone Audibility value is registered according to the Pure Tone Audibilities. If a Pure and Complex Tones Audibilities have the same Pitch Category, the greater Audibility value is registered.

The Sonorounsess/Tonalness is a psychoacoustic measure of Pitch clarity that can be combined with Roughness to compute Sensory Consonance. The Pure Tonalness is a normalized qudratic sum of Pure Tone Audibilities Ap(P). And the Complex Tonalness is derived from the maximum Complex Tone Audibility Ac(P).

Pitch Multiplicity is the number of tones consciously perceived in a sonority. It is calculated from Pitch Audibility A(P), which is given by the maximum value of Pure and Complex Tone Audibilities for each Pitch Category. The Tone Salience is the probability of consciously perceiving a pitch and depends on Pitch Audibility and Pitch Multiplicity. 

Successive Pitch Relationships are given by Pitch Commonality and Pitch Distance, both dependant on the Tone Salience output. The Pitch Commonality is a Pearson correlation coefficient between the Tone Salience's pitch categories. The relationship increases according to common Tone Saliences and is equal to 1 in the case of equal spectra and -1 in the case of complementary sonorities. 

The Pitch Distance considers the probability of noticing a pitch from one sonority in another, but takes into account all possible intervals between perceived pitches in two sonorities.  The distance is zero in the case of identical sonorities, and exceeds zero otherwise. In the case of Pure Tones, the Pitch Distance is equal to the interval between them in semitones.
