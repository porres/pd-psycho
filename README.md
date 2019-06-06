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

About: Tonalness / Tone Salience  /Multiplicity / Pitch Commonality / Pitch Distance: 

These are from a Pitch Model by Parncutt & Strasburger (1994), based on the theoretical work of Ernst Terhardt in Pitch & Consonance, specially Pure/Complex tone Sensation and Sonorousness (a measure of Pitch clarity). The main outputs of the model are: Tonalness (same as sonorousness), Pitch Commonality and Pitch Distance. The Pitch Commonality is a measure of the degree in which two sonorities evoke common pitches. The Pitch Distance is a similar concept, but more pertinent in the melodic context. As the model considers an input spectrum so is not purely theoretical. It is also not generalized in octaves.

The input spectrum are lists of frequency values and linear amplitudes.  The frequency list needs to be converted to Pitch Category (which is basically a logarythmic MIDI scale). The amplitude lists needs to fist be converted to Auditory Levels Y(L) - defined as dB over the partial’s Level Threshold L(TH). The next stage converts to Audible Level AL(P), which is the Audible level after a masking stage. This is the used to calculate Pure Tone Audibilities Ap(P) also defined as the Spectral Pitch Weigth.

With this new input we can recognice a harmonic pattern. The Complex Tone Audibility Ac(P) (or Virtual Pitch Weight) is a degree to which a Harmonic Series or part of it is present in the spectrum. For that, a template of 10 harmonics is used, each with a different weigth. When there’s a match, the Complex Tone Audibility value is registered according to the Pure Tone Audibilities. If a Pure and Complex Tones Audibilities have the same Pitch Category, the greater Audibility value is registered.

The Sonorounsess/Tonalness is a psychoacoustic measure of Pitch clarity that can be combined with Roughness to compute Sensory Consonance. The Pure Tonalness is a normalized qudratic sim of Pure Tone Audibilities Ap(P). And the Complex Tonalness is derived from the maximum Complex Tone Audibility Ac(P).

Pitch Multiplicity is the number of tones consciously perceived in a sonority.  It is calculated from Pitch Audibility A(P), which is given by the maximum value of Pure and Complex Tone Audibilities for each Pitch Category. The Tone Salience is the probability of consciously perceiving a pitch and depends on Pitch Audibility and Pitch Multiplicity. 

Successive Pitch Relationships are given by Pitch Commonality and Pitch Distance, both depend on the Tone Salience output. The Pitch Commonality is a Pearson correlation coefficient between Tone Salience profiles. The relationship increases according to common Tone Saliences. It is equal to 1 in the case of equal spectra and -1 in the case of complementary sonorities. 

The Pitch Distance considers the probability of noticing a pitch from one sonority in another, but takes into account all possible intervals between pitches perceived in two sonorities. In the case of identical sonorities, the distance is zero, otherwise it always exceeds zero. In the case of Pure Tones, the Pitch Distance is equal to the interval between them in semitones.

Reference: 
- Parncutt, Richard & Strasburger, Hans. (1994). Applying psychoacoustics in composition: Harmonic progressions of non-harmonic sonorities. Perspectives of New Music. 32. 88-129. 10.2307/833600. 
