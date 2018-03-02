# meshReader
Mesh reader library for OpenCMISS-iron

## File formats

This library can convert the following mesh formats to OpenCMISS-iron node/element numbering schemes:

1. CHeart X/T/B file format
..* linear/quadratic triangles
..* linear/quadratic quadrilaterals
..* linear/quadratic tetrahedra
..* linear/quadratic hexahedra
2. Cubit file format
..* TBD

Further mesh formats, interpolation orders and mesh types will follow, e.g.,

1. exnode/exelem file format
2. UGRID file format
3. VTK file format

## Installation

Simply run the following commands:

1. cmake .
2. make
3. make install

