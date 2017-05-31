# Created by and for Qt Creator. This file was created for editing the project sources only.
# You may attempt to use it for building too, by modifying this file here.

#TARGET = CCD_benchmark

HEADERS = \
   $$PWD/abstractSPbasis.hpp \
   $$PWD/ccdClasses.hpp \
   $$PWD/electronGasInteraction.hpp \
   $$PWD/electronGasSPBasis.hpp \
   $$PWD/infMatterSPBasis.hpp \
   $$PWD/minnesotaInfMatter.hpp \
   $$PWD/pairingInteraction.hpp \
   $$PWD/pairingSPBasis.hpp \
   $$PWD/qDots2DSPBasis.hpp \
   $$PWD/SymBlock.hpp

SOURCES = \
   $$PWD/ccdClasses.cpp \
   $$PWD/electronGasInteraction.cpp \
   $$PWD/electronGasSPBasis.cpp \
   $$PWD/infMatterSPBasis.cpp \
   $$PWD/main.cpp \
   $$PWD/minnesotaInfMatter.cpp \
   $$PWD/pairingInteraction.cpp \
   $$PWD/pairingSPBasis.cpp \
   $$PWD/qDots2DSPBasis.cpp \
   $$PWD/SymBlock.cpp

LIBS += -L/usr/local/Cellar/OpenBLAS -lOpenBLAS

INCLUDEPATH += /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers


