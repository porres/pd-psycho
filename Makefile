############
# Makefile #
############

# library name

lib.name = psycho

# Sources:

al.class.sources := Classes/al.c
bark2hz.class.sources := Classes/bark2hz.c
commonality.class.sources := Classes/commonality.c
db2phon.class.sources := Classes/db2phon.c
distance.class.sources := Classes/distance.c
flunson.class.sources := Classes/flunson.c
harmonicity.class.sources := Classes/harmonicity.c
hz2bark.class.sources := Classes/hz2bark.c
hz2mel.class.sources := Classes/hz2mel.c
indigestibility.class.sources := Classes/indigestibility.c
iso226.class.sources := Classes/iso226.c
# iso226b.class.sources := Classes/iso226b.c
mel2hz.class.sources := Classes/mel2hz.c
phon2db.class.sources := Classes/phon2db.c
phon2sone.class.sources := Classes/phon2sone.c
roughness.class.sources := Classes/roughness.c
salience.class.sources := Classes/salience.c
sone2phon.class.sources := Classes/sone2phon.c
tonalness.class.sources := Classes/tonalness.c
yl.class.sources := Classes/yl.c

# extra files

datafiles = \
$(wildcard Help-files/*.pd) \
psycho-meta.pd \
README.md

# include Makefile.pdlibbuilder from submodule directory 'pd-lib-builder'
PDLIBBUILDER_DIR=pd-lib-builder/
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
