############
# Makefile #
############

# library name

lib.name = psycho

# Sources:

db2phon.class.sources := Classes/db2phon.c
flunson.class.sources := Classes/flunson.c
harmonicity.class.sources := Classes/harmonicity.c
hz2bark.class.sources := Classes/hz2bark.c
hz2mel.class.sources := Classes/hz2mel.c
indigestibility.class.sources := Classes/indigestibility.c
iso226.class.sources := Classes/iso226.c
mel2hz.class.sources := Classes/mel2hz.c
# iso226b.class.sources := Classes/iso226b.c
phon2db.class.sources := Classes/phon2db.c
phon2sone.class.sources := Classes/phon2sone.c
sone2phon.class.sources := Classes/sone2phon.c
roughcurve.class.sources := Classes/roughcurve.c
roughness.class.sources := Classes/roughness.c

# extra files

datafiles = \
$(wildcard Help-files/*.pd) \
psycho-meta.pd \
README.md

# include Makefile.pdlibbuilder from submodule directory 'pd-lib-builder'
PDLIBBUILDER_DIR=pd-lib-builder/
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
