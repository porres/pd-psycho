############
# Makefile #
############

# library name

lib.name = psycho

# Sources:

psycho.class.sources := Classes/psycho.c
hz2bark.class.sources := Classes/hz2bark.c
phon2sone.class.sources := Classes/phon2sone.c
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
