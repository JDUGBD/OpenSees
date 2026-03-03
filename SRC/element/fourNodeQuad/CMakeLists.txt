#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Element
    PRIVATE
    ConstantPressureVolumeQuad.cpp
    EnhancedQuad.cpp
    FourNodeQuad.cpp
    FourNodeQuad3d.cpp
    FourNodeQuadWithSensitivity.cpp
    EightNodeQuad.cpp
    NineNodeMixedQuad.cpp
    NineNodeQuad.cpp
    SixNodeTri.cpp
    quadEmbedded.cpp
    PUBLIC
    ConstantPressureVolumeQuad.h
    EnhancedQuad.h
    FourNodeQuad.h
    FourNodeQuad3d.h
    FourNodeQuadWithSensitivity.h
    EightNodeQuad.h
    NineNodeMixedQuad.h
    NineNodeQuad.h
    SixNodeTri.h
    quadEmbedded.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
