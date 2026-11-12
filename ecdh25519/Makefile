# ECDH25519 Project Makefile
# Modify the source lists below to add your own C and assembly files
#
# e.g. make run-qemu PLATFORM=qemu WINDOW=4 FE_IMPL=r8

# field arithmetic args
FE_IMPL ?= r25p5
FE_DEFINE_r8 = -DFE25519_R8
FE_DEFINE_r25p5 = -DFE25519_R25P5

# Add your C source files here:
PROJECT_C_SOURCES = test.c smult.c helper.c \
					fe25519_common.c fe25519_$(FE_IMPL).c group.c
# Add your assembly source files here:
PROJECT_ASM_SOURCES = 

# Convert sources to object file paths
PROJECT_C_OBJS = $(addprefix obj/,$(PROJECT_C_SOURCES:.c=.c.o))
PROJECT_ASM_OBJS = $(addprefix obj/,$(PROJECT_ASM_SOURCES:.S=.S.o))
PROJECT_OBJS = $(PROJECT_C_OBJS) $(PROJECT_ASM_OBJS)

# Include the common build infrastructure
include ../common/common.mk

# choose the field arithmetic impl by the argument
CFLAGS	+= $(FE_DEFINE_$(FE_IMPL)) -MMD -MP
ASFLAGS	+= -MMD -MP

-include $(PROJECT_OBJS:.o=.d)

# fixed window table generation args
PYTHON 	   := python3
WINDOW 	   ?= 4
LIMB_RADIX ?= 25.5
BASE_TABLE ?= table.h

TABLE_ARGS = --window $(WINDOW) --limb-radix $(LIMB_RADIX)

$(BASE_TABLE): gen_table.py
	$(PYTHON) $< $(TABLE_ARGS) --out $@

$(PROJECT_OBJS): | $(BASE_TABLE)

# project level clean up
.PHONY: clean-gen
clean: clean-gen

clean-gen:
	rm -f $(BASE_TABLE)
