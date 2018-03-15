MODULE meshReader
  USE ISO_C_BINDING
  IMPLICIT NONE
  INTEGER, PARAMETER :: INTG  = C_INT
  INTEGER, PARAMETER :: DP    = C_DOUBLE
CONTAINS
  !=============================================================================
  !
  ! This routine reads input mesh files in
  ! - CHeart X/T/B or X/T/S format
  ! - OpenCMISS NODE/ELEM/NSET format
  ! and stores them in arrays
  ! Note: this routine makes use of the Fortran 2008 feature to obtain a file unit that is free
  SUBROUTINE ReadMesh(Filename, Nodes, Elements, BoundaryPatches, Method, Err)
    ! IN / OUT variables
    CHARACTER(LEN=*),                   INTENT(IN)  :: Filename             !< The file name to import the mesh data from
    REAL(DP),           ALLOCATABLE,    INTENT(OUT) :: Nodes(:,:)           !< The coordinates of the mesh nodes
    INTEGER(INTG),      ALLOCATABLE,    INTENT(OUT) :: Elements(:,:)        !< The node IDs for each element
    INTEGER(INTG),      ALLOCATABLE,    INTENT(OUT) :: BoundaryPatches(:)   !< The boundary patch labels for all boundary nodes
    CHARACTER(LEN=*),                   INTENT(IN)  :: Method               !< The method to use, e.g, CHeart, Cubit, etc.
    INTEGER(INTG),                      INTENT(OUT) :: Err                  !< The error code.
    ! Local variables
    CHARACTER(LEN=256)                              :: VFileName,VMethod
    CHARACTER(LEN=256)                              :: VElementType
    CHARACTER(LEN=256)                              :: ElementType
    INTEGER(INTG)                                   :: FilenameLength,MethodLength,ElementTypeLength
    INTEGER(INTG)                                   :: NodeFileUnit,ElementFileUnit,BoundaryFileUnit
    INTEGER(INTG)                                   :: NumberOfNodes,NumberOfDimensions
    INTEGER(INTG)                                   :: NumberOfElements,NumberOfNodesPerElement
    INTEGER(INTG)                                   :: NumberOfBoundaryPatches,NumberOfBoundaryPatchComponents
    INTEGER(INTG)                                   :: IntValue,CurrentIdx,ComponentIdx,PatchIdx,CurrentPatchID,Offset,Idx
    INTEGER(INTG)                                   :: CharacterIdx,FirstIdx,PreviousIdx
    INTEGER(INTG), ALLOCATABLE                      :: Permutation(:),IntValuesT(:),IntValuesB(:),BoundaryPatchesTemp(:,:)
    INTEGER(INTG), ALLOCATABLE                      :: LengthOfNodesets(:)
    INTEGER(INTG)                                   :: NumberOfNodesets
    INTEGER(INTG)                                   :: NumberOfPatchIDs
    INTEGER(INTG)                                   :: PatchIDs(25),NumberOfNodesPerPatchID(25),CurrentFirstPatchIdx(25)
    INTEGER(INTG)                                   :: NumberOfNodesT,NumberOfNodesB
    INTEGER(INTG)                                   :: FirstLineOfElementFile(27)
    INTEGER(INTG)                                   :: i, j
    LOGICAL                                         :: FileExists
    CHARACTER(LEN=256)                              :: Line

    !===========================================================================
    ! Check, if intent(out) variables are already allocated
    IF(ALLOCATED(Nodes)) THEN
      WRITE(*,*) "Nodes already allocated."
      STOP
    END IF
    IF(ALLOCATED(Elements)) THEN
      WRITE(*,*) "Elements already allocated."
      STOP
    END IF
    IF(ALLOCATED(BoundaryPatches)) THEN
      WRITE(*,*) "BoundaryPatches already allocated."
      STOP
    END IF
    !===========================================================================
    ! Initialize variables
    PatchIDs                = -1_INTG ! default, not present
    NumberOfNodesPerPatchID =  0_INTG
    NumberOfPatchIDs        =  0_INTG
    NumberOfDimensions      = 0_INTG
    NumberOfNodes           = 0_INTG
    NumberOfElements        = 0_INTG
    NumberOfBoundaryPatches = -1_INTG
    NumberOfNodesPerElement = 0_INTG
    NumberOfNodesT          = 0_INTG
    NumberOfNodesets        = 0_INTG
    !===========================================================================
    ! Get file name and method name
    FilenameLength  = LEN_TRIM(Filename)
    VFilename       = Filename(1:FilenameLength)
    MethodLength    = LEN_TRIM(Method)
    VMethod         = Method(1:MethodLength)
    !===========================================================================
    ! Reading the X/T/B or X/T/S files from CHeart file format
    IF(VMethod=="CHeart") THEN
      !=========================================================================
      ! Get available file units (Fortran 2008 feature) and read first line
      OPEN(NEWUNIT=NodeFileUnit,     FILE=TRIM(VFilename)//".X", ACTION="read")
      OPEN(NEWUNIT=ElementFileUnit,  FILE=TRIM(VFilename)//".T", ACTION="read")
      READ(NodeFileUnit, *)     NumberOfNodes,NumberOfDimensions
      READ(ElementFileUnit, *)  NumberOfElements,NumberOfNodesT
      INQUIRE(FILE=TRIM(VFilename)//".B",EXIST=FileExists)
      IF(FileExists) THEN
        OPEN(NEWUNIT=BoundaryFileUnit, FILE=TRIM(VFilename)//".B", ACTION="read")
        READ(BoundaryFileUnit, *) NumberOfBoundaryPatches
      ELSE
        INQUIRE(FILE=TRIM(VFilename)//".S",EXIST=FileExists)
        IF(.NOT.FileExists) THEN
          WRITE(*,*) "File does not exist: "
          STOP
        END IF
        OPEN(NEWUNIT=BoundaryFileUnit, FILE=TRIM(VFilename)//".S", ACTION="read")
        READ(BoundaryFileUnit, *)  ElementType
        ! Get element name
        ElementTypeLength  = LEN_TRIM(ElementType)
        VElementType       = ElementType(1:ElementTypeLength)
      END IF
      !=========================================================================
      ! Do some sanity checks before reading mesh data
      IF(.NOT.(NumberOfNodes==NumberOfNodesT)) THEN
        WRITE(*,*) "X and T files have different number of nodes."
        STOP
      END IF
      IF(NumberOfNodes<0_INTG) THEN
        WRITE(*,*) "Invalid number of nodes."
        STOP
      END IF
      IF(NumberOfDimensions<0_INTG) THEN
        WRITE(*,*) "Invalid number of dimensions."
        STOP
      END IF
      IF(NumberOfElements<0_INTG) THEN
        WRITE(*,*) "Invalid number of elements."
        STOP
      END IF
      IF(NumberOfBoundaryPatches<0_INTG) THEN
        IF(.NOT.(VElementType=="TRI3"    .OR. &
          &      VElementType=="TRI6"    .OR. &
          &      VElementType=="QUAD4"   .OR. &
          &      VElementType=="QUAD9"   .OR. &
          &      VElementType=="TETRA4"  .OR. &
          &      VElementType=="TETRA10" .OR. &
          &      VElementType=="HEX8"    .OR. &
          &      VElementType=="HEX27")) THEN
          WRITE(*,*) "Invalid number of boundary patches or element type."
          STOP
        END IF
      END IF
      ! Read X/T/B-files
      IF(NumberOfBoundaryPatches>0) THEN
          !=========================================================================
          ! To get the number of nodes per element, we have to do a little trick:
          ! We read the first non-header line of the element file and read each integer value
          READ(ElementFileUnit,'(A)') Line
          PreviousIdx   = 1
          FirstIdx      = 1
          DO CharacterIdx=1,LEN(Line)
            CurrentIdx = INDEX('0123456789', Line(CharacterIdx:CharacterIdx))
            IF((CurrentIdx==0_INTG).AND.(PreviousIdx>0_INTG)) THEN
              READ(Line(FirstIdx:CharacterIdx-1), *) IntValue
              NumberOfNodesPerElement                           = NumberOfNodesPerElement + 1_INTG
              FirstLineOfElementFile(NumberOfNodesPerElement)   = IntValue
            ELSE IF((CurrentIdx>0_INTG).AND.(PreviousIdx==0_INTG)) THEN
              FirstIdx = CharacterIdx
            END IF
            PreviousIdx = CurrentIdx
          END DO
          ALLOCATE(IntValuesT(NumberOfNodesPerElement),STAT=Err)
          IF(Err/=0) THEN
            WRITE(*,*) "Could not allocate memory."
            STOP
          END IF
          ! Close files and re-open them to start reading at beginning
          CLOSE(NodeFileUnit)
          CLOSE(ElementFileUnit)
          CLOSE(BoundaryFileUnit)
          OPEN(NEWUNIT=NodeFileUnit,     FILE=TRIM(VFilename)//".X", ACTION="read")
          OPEN(NEWUNIT=ElementFileUnit,  FILE=TRIM(VFilename)//".T", ACTION="read")
          OPEN(NEWUNIT=BoundaryFileUnit, FILE=TRIM(VFilename)//".B", ACTION="read")
          !=========================================================================
          ! Figure out:
          ! - how many components the boundary file has
          ! - the mapping between different node order conventions CHeart -> OpenCMISS-iron
          ALLOCATE(Permutation(NumberOfNodesPerElement),STAT=Err)
          IF(Err/=0) THEN
            WRITE(*,*) "Could not allocate memory."
            STOP
          END IF
          SELECT CASE(NumberOfDimensions)
          CASE(2)
            SELECT CASE(NumberOfNodesPerElement)
            CASE(3)
              ! Linear triangle
              Permutation       = 0
              DO CurrentIdx=1,NumberOfNodesPerElement
                Permutation(CurrentIdx) = CurrentIdx
              END DO
              NumberOfBoundaryPatchComponents = 4_INTG
            CASE(4)
              ! Linear quadrilateral
              Permutation       = 0
              DO CurrentIdx=1,NumberOfNodesPerElement
                Permutation(CurrentIdx) = CurrentIdx
              END DO
              NumberOfBoundaryPatchComponents = 4_INTG
            CASE(6)
              ! Quadratic triangle
              Permutation       = 0
              Permutation(1)    = 1
              Permutation(2)    = 2
              Permutation(3)    = 3
              Permutation(4)    = 4
              Permutation(5)    = 6
              Permutation(6)    = 5
              NumberOfBoundaryPatchComponents = 5_INTG
            CASE(9)
              ! Quadratic quadrilateral
              Permutation       = 0
              Permutation(1)    = 1
              Permutation(2)    = 5
              Permutation(3)    = 2
              Permutation(4)    = 6
              Permutation(5)    = 7
              Permutation(6)    = 8
              Permutation(7)    = 3
              Permutation(8)    = 9
              Permutation(9)    = 4
              NumberOfBoundaryPatchComponents = 5_INTG
            CASE DEFAULT
              WRITE(*,*) "Unknown 2D mesh type for method: "//TRIM(Method)
              STOP
            END SELECT
          CASE(3)
            SELECT CASE(NumberOfNodesPerElement)
            CASE(4)
              ! Linear tetrahedron
              Permutation       = 0
              DO CurrentIdx=1,NumberOfNodesPerElement
                Permutation(CurrentIdx) = CurrentIdx
              END DO
              NumberOfBoundaryPatchComponents = 5_INTG
            CASE(8)
              ! Linear hexahedron
              Permutation       = 0
              DO CurrentIdx=1,NumberOfNodesPerElement
                Permutation(CurrentIdx) = CurrentIdx
              END DO
              NumberOfBoundaryPatchComponents = 6_INTG
            CASE(10)
              ! Quadratic tetrahedron
              Permutation       =  0
              Permutation(1)    =  1
              Permutation(2)    =  2
              Permutation(3)    =  3
              Permutation(4)    =  4
              Permutation(5)    =  5
              Permutation(6)    =  6
              Permutation(7)    =  8
              Permutation(8)    =  7
              Permutation(9)    = 10
              Permutation(10)   =  9
              NumberOfBoundaryPatchComponents = 8_INTG
            CASE(27)
              ! Quadratic hexahedron
              Permutation       =  0
              Permutation(1)    =  1
              Permutation(2)    =  9
              Permutation(3)    =  2
              Permutation(4)    = 10
              Permutation(5)    = 11
              Permutation(6)    = 12
              Permutation(7)    =  3
              Permutation(8)    = 13
              Permutation(9)    =  4
              Permutation(10)   = 14
              Permutation(11)   = 15
              Permutation(12)   = 16
              Permutation(13)   = 17
              Permutation(14)   = 18
              Permutation(15)   = 19
              Permutation(16)   = 20
              Permutation(17)   = 21
              Permutation(18)   = 22
              Permutation(19)   =  5
              Permutation(20)   = 23
              Permutation(21)   =  6
              Permutation(22)   = 24
              Permutation(23)   = 25
              Permutation(24)   = 26
              Permutation(25)   =  7
              Permutation(26)   = 27
              Permutation(27)   =  8
              NumberOfBoundaryPatchComponents = 11_INTG
            CASE DEFAULT
              WRITE(*,*) "Unknown 3D mesh type for method: "//TRIM(Method)
              STOP
            END SELECT
          CASE DEFAULT
            WRITE(*,*) "1D mesh import not supported for method: "//TRIM(Method)
            STOP
        END SELECT
        !=========================================================================
        ! Now, allocate node, elements, boundary patch variables
        ALLOCATE(Nodes(NumberOfNodes,NumberOfDimensions),STAT=Err)
        IF(Err/=0) THEN
          WRITE(*,*) "Can not allocate memory."
          STOP
        END IF
        ALLOCATE(Elements(NumberOfElements,NumberOfNodesPerElement),STAT=Err)
        IF(Err/=0) THEN
          WRITE(*,*) "Can not allocate memory."
          STOP
        END IF
        ALLOCATE(IntValuesB(NumberOfBoundaryPatchComponents),STAT=Err)
        IF(Err/=0) THEN
          WRITE(*,*) "Could not allocate memory."
          STOP
        END IF
        ALLOCATE(BoundaryPatchesTemp(NumberOfBoundaryPatches,NumberOfBoundaryPatchComponents),STAT=Err)
        IF(Err/=0) THEN
          WRITE(*,*) "Could not allocate memory."
          STOP
        END IF
        !=========================================================================
        ! Skip header line and read all other lines in the node file
        READ(NodeFileUnit,*) IntValue
        DO CurrentIdx=1,NumberOfNodes
          READ(NodeFileUnit,*) Nodes(CurrentIdx,:)
        END DO
        ! Skip header line and read all other lines in the element file
        READ(ElementFileUnit,*) IntValue
        DO CurrentIdx=1,NumberOfElements
          READ(ElementFileUnit,*) IntValuesT(:)
          DO Idx=1,NumberOfNodesPerElement
            Elements(CurrentIdx,Idx)  = IntValuesT(Permutation(Idx))
          END DO
        END DO
        !=========================================================================
        ! Skip header line and read all other lines in the boundary file
        READ(BoundaryFileUnit,*) IntValue
        ! First of all, let's read all the boundary patches in a temporary variable and count the number of unique boundary patch IDs
        DO CurrentIdx=1,NumberOfBoundaryPatches
          READ(BoundaryFileUnit,*) BoundaryPatchesTemp(CurrentIdx,:)
          DO PatchIdx=1,SIZE(PatchIDs,1)
            ! Check whether it is an existing patch ID or if we found a new one
            IF(PatchIDs(PatchIdx)==BoundaryPatchesTemp(CurrentIdx,NumberOfBoundaryPatchComponents)) THEN
              NumberOfNodesPerPatchID(PatchIdx)   = NumberOfNodesPerPatchID(PatchIdx) + NumberOfBoundaryPatchComponents - 2_INTG
              EXIT
            ELSE IF(PatchIDs(PatchIdx)==-1_INTG) THEN
              NumberOfPatchIDs                    = NumberOfPatchIDs + 1_INTG
              PatchIDs(PatchIdx)                  = BoundaryPatchesTemp(CurrentIdx,NumberOfBoundaryPatchComponents)
              NumberOfNodesPerPatchID(PatchIdx)   = NumberOfNodesPerPatchID(PatchIdx) + NumberOfBoundaryPatchComponents - 2_INTG
              EXIT
            ELSE
              ! Do nothing
            END IF
          END DO
        END DO
        ! Now transform boundary patches to a common format
        ALLOCATE(BoundaryPatches(1+2*NumberOfPatchIDs+NumberOfBoundaryPatches*(NumberOfBoundaryPatchComponents-2_INTG)),STAT=Err)
        IF(Err/=0) THEN
          WRITE(*,*) "Could not allocate memory."
          STOP
        END IF
        ! Set the number of patch IDs
        BoundaryPatches(1)    = NumberOfPatchIDs
        CurrentFirstPatchIdx  = 1_INTG+2*NumberOfPatchIDs+1_INTG
        DO CurrentIdx=1,NumberOfPatchIDs
          ! Get number of nodes for each patch ID
          BoundaryPatches(1+CurrentIdx)                   = NumberOfNodesPerPatchID(CurrentIdx)
          ! Get the current patch ID
          BoundaryPatches(1+NumberOfPatchIDs+CurrentIdx)  = PatchIDs(CurrentIdx)
          ! Define starting indices for boundary patches
          IF(CurrentIdx>1) CurrentFirstPatchIdx(CurrentIdx) = CurrentFirstPatchIdx(CurrentIdx-1_INTG) + BoundaryPatches(CurrentIdx)
        END DO
        DO CurrentIdx=1,NumberOfBoundaryPatches
          ! Get the current patch ID to match
          CurrentPatchID=BoundaryPatchesTemp(CurrentIdx,NumberOfBoundaryPatchComponents)
          DO PatchIdx=1,NumberOfPatchIDs
            IF(CurrentPatchID==PatchIDs(PatchIdx)) THEN
              Offset=CurrentFirstPatchIdx(PatchIdx)
              EXIT
            END IF
          END DO
          DO ComponentIdx=2,NumberOfBoundaryPatchComponents-1
            BoundaryPatches(Offset+ComponentIdx-2)    = BoundaryPatchesTemp(CurrentIdx,ComponentIdx)
          END DO
         CurrentFirstPatchIdx(PatchIdx)=CurrentFirstPatchIdx(PatchIdx)+NumberOfBoundaryPatchComponents-2_INTG
        END DO
      ! Read X/T/S-files
      ELSE
        ! Specifing number of nodes per element, interpolation type (linear/quadratic) and the permutation vector based on element type
        IF (VElementType=="TRI3") THEN
          NumberOfNodesPerElement = 3
          ALLOCATE(Permutation(NumberOfNodesPerElement),STAT=Err)
          IF(Err/=0) THEN
            WRITE(*,*) "Could not allocate memory."
            STOP
          END IF
          Permutation       = 0
          DO CurrentIdx=1,NumberOfNodesPerElement
            Permutation(CurrentIdx) = CurrentIdx
          END DO
        ELSE IF (VElementType=="TRI6") THEN
          NumberOfNodesPerElement = 6
          ALLOCATE(Permutation(NumberOfNodesPerElement),STAT=Err)
          IF(Err/=0) THEN
            WRITE(*,*) "Could not allocate memory."
            STOP
          END IF
          Permutation       = 0
          Permutation(1)    = 1
          Permutation(2)    = 2
          Permutation(3)    = 3
          Permutation(4)    = 4
          Permutation(5)    = 6
          Permutation(6)    = 5
        ELSE IF (VElementType=="QUAD4") THEN
          NumberOfNodesPerElement = 4
          ALLOCATE(Permutation(NumberOfNodesPerElement),STAT=Err)
          IF(Err/=0) THEN
            WRITE(*,*) "Could not allocate memory."
            STOP
          END IF
          Permutation       = 0
          DO CurrentIdx=1,NumberOfNodesPerElement
            Permutation(CurrentIdx) = CurrentIdx
          END DO
        ELSE IF (VElementType=="QUAD9") THEN
          NumberOfNodesPerElement = 9
          ALLOCATE(Permutation(NumberOfNodesPerElement),STAT=Err)
          IF(Err/=0) THEN
            WRITE(*,*) "Could not allocate memory."
            STOP
          END IF
          Permutation       = 0
          Permutation(1)    = 1
          Permutation(2)    = 5
          Permutation(3)    = 2
          Permutation(4)    = 6
          Permutation(5)    = 7
          Permutation(6)    = 8
          Permutation(7)    = 3
          Permutation(8)    = 9
          Permutation(9)    = 4
        ELSE IF (VElementType=="TETRA4") THEN
          NumberOfNodesPerElement = 4
          ALLOCATE(Permutation(NumberOfNodesPerElement),STAT=Err)
          IF(Err/=0) THEN
            WRITE(*,*) "Could not allocate memory."
            STOP
          END IF
          Permutation       = 0
          DO CurrentIdx=1,NumberOfNodesPerElement
            Permutation(CurrentIdx) = CurrentIdx
          END DO
        ELSE IF (VElementType=="TETRA10") THEN
          NumberOfNodesPerElement = 10
          ALLOCATE(Permutation(NumberOfNodesPerElement),STAT=Err)
          IF(Err/=0) THEN
            WRITE(*,*) "Could not allocate memory."
            STOP
          END IF
          Permutation       =  0
          Permutation(1)    =  1
          Permutation(2)    =  2
          Permutation(3)    =  3
          Permutation(4)    =  4
          Permutation(5)    =  5
          Permutation(6)    =  6
          Permutation(7)    =  8
          Permutation(8)    =  7
          Permutation(9)    = 10
          Permutation(10)   =  9
        ELSE IF (VElementType=="HEX8") THEN
          NumberOfNodesPerElement = 8
          ALLOCATE(Permutation(NumberOfNodesPerElement),STAT=Err)
          IF(Err/=0) THEN
            WRITE(*,*) "Could not allocate memory."
            STOP
          END IF
          Permutation       = 0
          DO CurrentIdx=1,NumberOfNodesPerElement
            Permutation(CurrentIdx) = CurrentIdx
          END DO
        ELSE IF (VElementType=="HEX27") THEN
          NumberOfNodesPerElement = 27
          ALLOCATE(Permutation(NumberOfNodesPerElement),STAT=Err)
          IF(Err/=0) THEN
            WRITE(*,*) "Could not allocate memory."
            STOP
          END IF
          Permutation       =  0
          Permutation(1)    =  1
          Permutation(2)    =  9
          Permutation(3)    =  2
          Permutation(4)    = 10
          Permutation(5)    = 11
          Permutation(6)    = 12
          Permutation(7)    =  3
          Permutation(8)    = 13
          Permutation(9)    =  4
          Permutation(10)   = 14
          Permutation(11)   = 15
          Permutation(12)   = 16
          Permutation(13)   = 17
          Permutation(14)   = 18
          Permutation(15)   = 19
          Permutation(16)   = 20
          Permutation(17)   = 21
          Permutation(18)   = 22
          Permutation(19)   =  5
          Permutation(20)   = 23
          Permutation(21)   =  6
          Permutation(22)   = 24
          Permutation(23)   = 25
          Permutation(24)   = 26
          Permutation(25)   =  7
          Permutation(26)   = 27
          Permutation(27)   =  8
        END IF
        ALLOCATE(IntValuesT(NumberOfNodesPerElement),STAT=Err)
        IF(Err/=0) THEN
          WRITE(*,*) "Could not allocate memory."
          STOP
        END IF
        !=========================================================================
        ! Now, allocate node, elements variables
        ALLOCATE(Nodes(NumberOfNodes,NumberOfDimensions),STAT=Err)
        IF(Err/=0) THEN
          WRITE(*,*) "Can not allocate memory."
          STOP
        END IF
        ALLOCATE(Elements(NumberOfElements,NumberOfNodesPerElement),STAT=Err)
        IF(Err/=0) THEN
          WRITE(*,*) "Can not allocate memory."
          STOP
        END IF
        !=========================================================================
        ! Read nodes, elements, node sets
        DO i = 1, NumberOfNodes
          READ (NodeFileUnit,*) Nodes(i,:)
        END DO
        DO i = 1, NumberOfElements
          READ (ElementFileUnit,*)  IntValuesT(:)
          DO j=1,NumberOfNodesPerElement
            Elements(i,j)  = IntValuesT(Permutation(j))
          END DO
        END DO
        ! Reading the nodesets
        READ(BoundaryFileUnit, *)  NumberOfNodesets
        IF (NumberOfNodesets > 0) THEN
          ALLOCATE(LengthOfNodesets(NumberOfNodesets))
          READ (BoundaryFileUnit,*) LengthOfNodesets(:)
          ALLOCATE(BoundaryPatches(1+2*NumberOfNodesets+SUM(LengthOfNodesets)))

          BoundaryPatches(1)                    = NumberOfNodesets
          BoundaryPatches(2:1+NumberOfNodesets) = LengthOfNodesets(:)
          DO i = 1, NumberOfNodesets
            BoundaryPatches(1+NumberOfNodesets+i) = i
          END DO
          READ (BoundaryFileUnit,*)  BoundaryPatches(2+2*NumberOfNodesets:)
        END IF
      END IF
      !=========================================================================
      ! Close files
      CLOSE(NodeFileUnit)
      CLOSE(ElementFileUnit)
      CLOSE(BoundaryFileUnit)
    !===========================================================================
    ! Reading the NODE/ELEM/NSET files from OpenCMISS file format
    ELSE IF(VMethod=="OpenCMISS") THEN
      ! Get available file units (Fortran 2008 feature)
      OPEN(NEWUNIT=NodeFileUnit,     FILE=TRIM(VFilename)//".NODE", ACTION="read")
      OPEN(NEWUNIT=ElementFileUnit,  FILE=TRIM(VFilename)//".ELEM", ACTION="read")
      OPEN(NEWUNIT=BoundaryFileUnit, FILE=TRIM(VFilename)//".NSET", ACTION="read")

      ! Reading the headers with mesh data
      READ(NodeFileUnit, *)     NumberOfNodes, NumberOfDimensions
      READ(ElementFileUnit, *)  NumberOfElements, NumberOfNodesT
      READ(BoundaryFileUnit, *) ElementType

      ! Get element name
      ElementTypeLength  = LEN_TRIM(ElementType)
      VElementType       = ElementType(1:ElementTypeLength)

      ! Do some sanity checks before reading mesh data
      IF(.NOT.(NumberOfNodes==NumberOfNodesT)) THEN
        WRITE(*,*) "X and T files have different number of nodes."
        STOP
      END IF
      IF(NumberOfNodes<0_INTG) THEN
        WRITE(*,*) "Invalid number of nodes."
        STOP
      END IF
      IF(NumberOfDimensions<0_INTG) THEN
        WRITE(*,*) "Invalid number of dimensions."
        STOP
      END IF
      IF(NumberOfElements<0_INTG) THEN
        WRITE(*,*) "Invalid number of elements."
        STOP
      END IF
      IF(.NOT.(VElementType=="TRI3"    .OR. &
               VElementType=="TRI6"    .OR. &
               VElementType=="QUAD4"   .OR. &
               VElementType=="QUAD9"   .OR. &
               VElementType=="TETRA4"  .OR. &
               VElementType=="TETRA10" .OR. &
               VElementType=="HEX8"    .OR. &
               VElementType=="HEX27")) THEN
        WRITE(*,*) "Invalid element type."
        STOP
      END IF
      ! Reading the nodes
      ALLOCATE(Nodes(NumberOfNodes,NumberOfDimensions),STAT=Err)
      DO i = 1, NumberOfNodes
        READ (NodeFileUnit,*) Nodes(i,:)
      END DO

      ! Specifing number of nodes per element and interpolation type (linear/quadratic) based on element type
      IF (VElementType=="TRI3") THEN
        NumberOfNodesPerElement = 3
      ELSE IF (VElementType=="TRI6") THEN
        NumberOfNodesPerElement = 6
      ELSE IF (VElementType=="QUAD4") THEN
        NumberOfNodesPerElement = 4
      ELSE IF (VElementType=="QUAD9") THEN
        NumberOfNodesPerElement = 9
      ELSE IF (VElementType=="TETRA4") THEN
        NumberOfNodesPerElement = 4
      ELSE IF (VElementType=="TETRA10") THEN
        NumberOfNodesPerElement = 10
      ELSE IF (VElementType=="HEX8") THEN
        NumberOfNodesPerElement = 8
      ELSE IF (VElementType=="HEX27") THEN
        NumberOfNodesPerElement = 27
      END IF

      ALLOCATE(Elements(NumberOfElements,NumberOfNodesPerElement))
      DO i = 1, NumberOfElements
        READ (ElementFileUnit,*)  Elements(i,:)
      END DO

      ! Reading the nodesets
      READ(BoundaryFileUnit, *)  NumberOfNodesets
      IF (NumberOfNodesets > 0) THEN
        ALLOCATE(LengthOfNodesets(NumberOfNodesets))
        READ (BoundaryFileUnit,*) LengthOfNodesets(:)
        ALLOCATE(BoundaryPatches(1+2*NumberOfNodesets+SUM(LengthOfNodesets)))

        BoundaryPatches(1)                    = NumberOfNodesets
        BoundaryPatches(2:1+NumberOfNodesets) = LengthOfNodesets(:)
        DO i = 1, NumberOfNodesets
          BoundaryPatches(1+NumberOfNodesets+i) = i
        END DO
        READ (BoundaryFileUnit,*)  BoundaryPatches(2+2*NumberOfNodesets:)
      END IF

      ! Close files
      CLOSE(NodeFileUnit)
      CLOSE(ElementFileUnit)
      CLOSE(BoundaryFileUnit)
    !===========================================================================
    ! Unknown format, mesh reader not implemented
    ELSE
      WRITE(*,*) "Invalid mesh import type. Valid types are: CHeart, OpenCMISS."
      STOP
    END IF

    !===========================================================================
    ! Deallocate temporary variables
    IF(ALLOCATED(Permutation))          DEALLOCATE(Permutation)
    IF(ALLOCATED(IntValuesT))           DEALLOCATE(IntValuesT)
    IF(ALLOCATED(IntValuesB))           DEALLOCATE(IntValuesB)
    IF(ALLOCATED(BoundaryPatchesTemp))  DEALLOCATE(BoundaryPatchesTemp)
    IF(ALLOCATED(LengthOfNodesets))     DEALLOCATE(LengthOfNodesets)

    RETURN

  END SUBROUTINE ReadMesh

  !
  !=============================================================================
  !
  SUBROUTINE ImportedMesh_SurfaceWeightsGet(nodalWeights,numNodesPerElement,Err)
    ! argument variables
    REAL(DP),       ALLOCATABLE,  INTENT(INOUT) :: nodalWeights(:)      !< on return, the consistent nodal weights
    INTEGER(INTG),                INTENT(IN)    :: numNodesPerElement   !< the number of nodes of surface element
    INTEGER(INTG),                INTENT(OUT)   :: Err                  !< the error code.

    SELECT CASE(numNodesPerElement)
      ! tri-quadratic hexahedron
      CASE(27)
        IF(ALLOCATED(nodalWeights)) DEALLOCATE(nodalWeights)
        ! note: surface element has 9 nodes
        ALLOCATE(nodalWeights(9))
        nodalWeights(1) = 1.0_DP / 36.0_DP
        nodalWeights(5) = 1.0_DP /  9.0_DP
        nodalWeights(2) = 1.0_DP / 36.0_DP
        nodalWeights(6) = 1.0_DP /  9.0_DP
        nodalWeights(7) = 4.0_DP /  9.0_DP
        nodalWeights(8) = 1.0_DP /  9.0_DP
        nodalWeights(3) = 1.0_DP / 36.0_DP
        nodalWeights(9) = 1.0_DP /  9.0_DP
        nodalWeights(4) = 1.0_DP / 36.0_DP
      CASE DEFAULT
        WRITE(*,*) "Not implemented."
    END SELECT

  END SUBROUTINE ImportedMesh_SurfaceWeightsGet

  !
  !=============================================================================
  !
  SUBROUTINE GeneratedMesh_SurfaceWeightsGet(nodalWeights,faceNumber,numberGlobalXelements, &
    & numberGlobalYelements,numberGlobalZelements,interpolationType,Err)
    ! argument variables
    REAL(DP),       ALLOCATABLE,  INTENT(INOUT) :: nodalWeights(:)        !< on return, the consistent nodal weights
    INTEGER(INTG),      INTENT(IN)    :: faceNumber             !< the face number
    INTEGER(INTG),      INTENT(IN)    :: numberGlobalXelements  !< number of elements in x-direction
    INTEGER(INTG),      INTENT(IN)    :: numberGlobalYelements  !< number of elements in y-direction
    INTEGER(INTG),      INTENT(IN)    :: numberGlobalZelements  !< number of elements in z-direction
    INTEGER(INTG),      INTENT(IN)    :: interpolationType      !< interpolation type
    INTEGER(INTG),      INTENT(OUT)   :: Err                    !< the error code.
    ! local variables
    REAL(DP)      :: weight
    INTEGER(INTG) :: node_idx,node_idx_2
    INTEGER(INTG) :: numNodesX,numNodesY,numNodesZ,numNodesXY,numNodesXZ,numNodesYZ
    ! the surface numbers correspond to OpenCMISS-iron face numbers!
    INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_LEFT_SURFACE   =  8
    INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_RIGHT_SURFACE  =  9
    INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_TOP_SURFACE    = 10
    INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_BOTTOM_SURFACE = 11
    INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_FRONT_SURFACE  = 12
    INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_BACK_SURFACE   = 13

    ! deallocate nodal weights if already allocated
    IF(ALLOCATED(nodalWeights)) THEN
      DEALLOCATE(nodalWeights)
    END IF

    ! get the number of nodes in each coordinate direction
    IF(interpolationType == 1) THEN
      numNodesX     = 1 * numberGlobalXElements + 1
      numNodesY     = 1 * numberGlobalYElements + 1
      numNodesZ     = 1 * numberGlobalZElements + 1
    ELSE IF(interpolationType == 2) THEN
      numNodesX     = 2 * numberGlobalXElements + 1
      numNodesY     = 2 * numberGlobalYElements + 1
      numNodesZ     = 2 * numberGlobalZElements + 1
    ELSE
      WRITE(*,*) "Interpolation type ",interpolationType," not implemented."
    END IF
    numNodesXY    = numNodesX * numNodesY
    numNodesXZ    = numNodesX * numNodesZ
    numNodesYZ    = numNodesY * numNodesZ
    ! compute the nodal weights depending on which surface we are on
    !
    !         template
    !    __ __ __ __ __ __
    !   |     |     |     |
    !
    !   |     |     |     |
    !    __ __ __ __ __ __
    !   |     |     |     |
    !
    !   |     |     |     |
    !    __ __ __ __ __ __
    !   |     |     |     |
    !
    !   |     |     |     |
    !    __ __ __ __ __ __
    !
    SELECT CASE(faceNumber)
      CASE(GENERATED_MESH_REGULAR_LEFT_SURFACE,GENERATED_MESH_REGULAR_RIGHT_SURFACE)
        ! allocate and initialise nodal weights
        ALLOCATE(nodalWeights(numNodesYZ))
        nodalWeights  = 0.0_DP
        IF(interpolationType == 1) THEN
          IF(numberGlobalZElements<=0) THEN
            weight  = 1.0_DP / 2.0_DP
            nodalWeights(1)           = weight
            nodalWeights(numNodesYZ)  = weight
            weight  = 2.0_DP / 2.0_DP
            DO node_idx=2,numNodesYZ-1
              nodalWeights(node_idx)  = weight
            END DO
          ELSE
            !
            !   X__ __ __ __ __ __X
            !   |     |     |     |
            !
            !   |     |     |     |
            !    __ __ __ __ __ __
            !   |     |     |     |
            !
            !   |     |     |     |
            !    __ __ __ __ __ __
            !   |     |     |     |

            !   |     |     |     |
            !   X__ __ __ __ __ __X
            !
            weight                                = 1.0_DP / 4.0_DP
            nodalWeights(1)                       = weight
            nodalWeights(numNodesY)               = weight
            nodalWeights(numNodesYZ-numNodesY+1)  = weight
            nodalWeights(numNodesYZ)              = weight
            !
            !    __ __X__ __X__ __
            !   |     |     |     |
            !
            !   |     |     |     |
            !   X__ __ __ __ __ __X
            !   |     |     |     |
            !
            !   |     |     |     |
            !   X__ __ __ __ __ __X
            !   |     |     |     |
            !
            !   |     |     |     |
            !    __ __X__ __X__ __
            !
            weight                                = 2.0_DP / 4.0_DP
            ! bottom
            DO node_idx=2,numNodesY-1,1
              nodalWeights(node_idx)  = weight
            END DO
            ! top
            DO node_idx=numNodesYZ-numNodesY+2,numNodesYZ-1,1
              nodalWeights(node_idx)  = weight
            END DO
            ! left/right
            DO node_idx=1,1*numberGlobalZelements-1,1
             nodalWeights(node_idx*numNodesY+1)   = weight
             nodalWeights((node_idx+1)*numNodesY) = weight
            END DO
            !    __ __ __ __ __ __
            !   |     |     |     |
            !
            !   |     |     |     |
            !    __ __X__ __X__ __
            !   |     |     |     |
            !
            !   |     |     |     |
            !    __ __X__ __X__ __
            !   |     |     |     |
            !
            !   |     |     |     |
            !    __ __ __ __ __ __
            !
            weight                                = 4.0_DP / 4.0_DP
            DO node_idx=1,1*numberGlobalZelements-1,1
              DO node_idx_2=node_idx*numNodesY+2,(node_idx+1)*numNodesY-1,1
                nodalWeights(node_idx_2)  = weight
              END DO
            END DO
            ! we are done.. do a sanity check
            IF(ABS(SUM(nodalWeights)-numberGlobalYelements*numberGlobalZelements)>1.0e-11) THEN
              WRITE(*,*) ""
              WRITE(*,*) "Invalid nodal weights for computing consistent nodal forces."
              WRITE(*,*) " SIZE(nodalWeights),numNodesY,numNodesZ,SUM(nodalWeights) = ", &
                & SIZE(nodalWeights),numNodesY,numNodesZ,SUM(nodalWeights)
              WRITE(*,*) ""
              DO node_idx=1,numNodesZ
                WRITE(*,*) " ",nodalWeights(1+(node_idx-1)*numNodesY:node_idx*numNodesY:1)
              END DO
              WRITE(*,*) ""
              STOP
            END IF
          END IF
        ELSE IF(interpolationType == 2) THEN
          IF(numberGlobalZElements<=0) THEN
            weight  = 1.0_DP / 6.0_DP
            nodalWeights(1)           = weight
            nodalWeights(numNodesYZ)  = weight
            weight  = 2.0_DP / 3.0_DP
            DO node_idx=2,numNodesYZ-1,2
              nodalWeights(node_idx)  = weight
            END DO
            weight  = 1.0_DP / 3.0_DP
            DO node_idx=3,numNodesYZ-2,2
              nodalWeights(node_idx)  = weight
            END DO
          ELSE
            !
            !   X__ __ __ __ __ __X
            !   |     |     |     |
            !
            !   |     |     |     |
            !    __ __ __ __ __ __
            !   |     |     |     |
            !
            !   |     |     |     |
            !    __ __ __ __ __ __
            !   |     |     |     |

            !   |     |     |     |
            !   X__ __ __ __ __ __X
            !
            weight                                = 1.0_DP / 36.0_DP
            nodalWeights(1)                       = weight
            nodalWeights(numNodesY)               = weight
            nodalWeights(numNodesYZ-numNodesY+1)  = weight
            nodalWeights(numNodesYZ)              = weight
            !
            !    __ __X__ __X__ __
            !   |     |     |     |
            !
            !   |     |     |     |
            !   X__ __ __ __ __ __X
            !   |     |     |     |
            !
            !   |     |     |     |
            !   X__ __ __ __ __ __X
            !   |     |     |     |
            !
            !   |     |     |     |
            !    __ __X__ __X__ __
            !
            weight                                = 2.0_DP / 36.0_DP
            ! bottom
            DO node_idx=3,numNodesY-2,2
              nodalWeights(node_idx)  = weight
            END DO
            ! top
            DO node_idx=numNodesYZ-numNodesY+3,numNodesYZ-2,2
              nodalWeights(node_idx)  = weight
            END DO
            ! left/right
            DO node_idx=2,2*numberGlobalZelements-2,2
             nodalWeights(node_idx*numNodesY+1)   = weight
             nodalWeights((node_idx+1)*numNodesY) = weight
            END DO
            !    __ __ __ __ __ __
            !   |     |     |     |
            !
            !   |     |     |     |
            !    __ __X__ __X__ __
            !   |     |     |     |
            !
            !   |     |     |     |
            !    __ __X__ __X__ __
            !   |     |     |     |
            !
            !   |     |     |     |
            !    __ __ __ __ __ __
            !
            weight                                = 4.0_DP / 36.0_DP
            DO node_idx=2,2*numberGlobalZelements-2,2
              DO node_idx_2=node_idx*numNodesY+3,(node_idx+1)*numNodesY-2,2
                nodalWeights(node_idx_2)  = weight
              END DO
            END DO
            !    __ __ __ __ __ __
            !   |     |     |     |
            !      X     X     X
            !   |     |     |     |
            !    __ __ __ __ __ __
            !   |     |     |     |
            !      X     X     X
            !   |     |     |     |
            !    __ __ __ __ __ __
            !   |     |     |     |
            !      X     X     X
            !   |     |     |     |
            !    __ __ __ __ __ __
            !
            weight                                = 4.0_DP / 9.0_DP
            DO node_idx=1,2*numberGlobalZelements-1,2
              DO node_idx_2=node_idx*numNodesY+2,(node_idx+1)*numNodesY-1,2
                nodalWeights(node_idx_2)  = weight
              END DO
            END DO
            !
            !    __X__ __X__ __X__
            !   |     |     |     |
            !   X                 X
            !   |     |     |     |
            !    __ __ __ __ __ __
            !   |     |     |     |
            !   X                 X
            !   |     |     |     |
            !    __ __ __ __ __ __
            !   |     |     |     |
            !   X                 X
            !   |     |     |     |
            !    __X__ __X__ __X__
            !
            weight                                = 1.0_DP / 9.0_DP
            ! top/bottom
            DO node_idx=2,numNodesY-1,2
              nodalWeights(node_idx)                      = weight
              nodalWeights(numNodesYZ-numNodesY+node_idx) = weight
            END DO
            ! left/right
            DO node_idx=1,2*numberGlobalZelements-1,2
              nodalWeights(node_idx*numNodesY+1)    = weight
              nodalWeights((node_idx+1)*numNodesY)  = weight
            END DO
            !    __ __ __ __ __ __
            !   |     |     |     |
            !         X     X
            !   |     |     |     |
            !    __X__ __X__ __X__
            !   |     |     |     |
            !         X     X
            !   |     |     |     |
            !    __X__ __X__ __X__
            !   |     |     |     |
            !         X     X
            !   |     |     |     |
            !    __X__ __X__ __X__
            !
            weight                                = 2.0_DP / 9.0_DP
            ! vertical edges
            DO node_idx=1,2*numberGlobalZelements-1,2
              DO node_idx_2=node_idx*numNodesY+3,(node_idx+1)*numNodesY-2,2
                nodalWeights(node_idx_2)  = weight
              END DO
            END DO
            ! horizontal edges
            DO node_idx=2,2*numberGlobalZelements-2,2
              DO node_idx_2=node_idx*numNodesY+2,(node_idx+1)*numNodesY-1,2
                nodalWeights(node_idx_2)  = weight
              END DO
            END DO
            ! we are done.. do a sanity check
            IF(ABS(SUM(nodalWeights)-numberGlobalYelements*numberGlobalZelements)>1.0e-11) THEN
              WRITE(*,*) ""
              WRITE(*,*) "Invalid nodal weights for computing consistent nodal forces."
              WRITE(*,*) " SIZE(nodalWeights),numNodesY,numNodesZ,SUM(nodalWeights) = ", &
                & SIZE(nodalWeights),numNodesY,numNodesZ,SUM(nodalWeights)
              WRITE(*,*) ""
              DO node_idx=1,numNodesZ
                WRITE(*,*) " ",nodalWeights(1+(node_idx-1)*numNodesY:node_idx*numNodesY:1)
              END DO
              WRITE(*,*) ""
              STOP
            END IF
          END IF
        END IF
        nodalWeights=nodalWeights/SUM(nodalWeights)
      CASE(GENERATED_MESH_REGULAR_TOP_SURFACE,GENERATED_MESH_REGULAR_BOTTOM_SURFACE)
        ! make sure that there are some nodes on this surface
        IF((numberGlobalXElements<=0).OR.(numberGlobalYElements<=0)) THEN
          WRITE(*,*) ">>>GeneratedMesh_SurfaceWeightsGet: There are no nodes on this surface. Face number: ",faceNumber
          STOP
        END IF
        ! allocate and initialise nodal weights
        ALLOCATE(nodalWeights(numNodesXY))
        nodalWeights  = 0.0_DP
        !
        !   X__ __ __ __ __ __X
        !   |     |     |     |
        !
        !   |     |     |     |
        !    __ __ __ __ __ __
        !   |     |     |     |
        !
        !   |     |     |     |
        !    __ __ __ __ __ __
        !   |     |     |     |

        !   |     |     |     |
        !   X__ __ __ __ __ __X
        !
        weight                                = 1.0_DP / 36.0_DP
        nodalWeights(1)                       = weight
        nodalWeights(numNodesX)               = weight
        nodalWeights(numNodesXY-numNodesX+1)  = weight
        nodalWeights(numNodesXY)              = weight
        !
        !    __ __X__ __X__ __
        !   |     |     |     |
        !
        !   |     |     |     |
        !   X__ __ __ __ __ __X
        !   |     |     |     |
        !
        !   |     |     |     |
        !   X__ __ __ __ __ __X
        !   |     |     |     |
        !
        !   |     |     |     |
        !    __ __X__ __X__ __
        !
        weight                                = 2.0_DP / 36.0_DP
        ! bottom
        DO node_idx=3,numNodesX-2,2
          nodalWeights(node_idx)  = weight
        END DO
        ! top
        DO node_idx=numNodesXY-numNodesX+3,numNodesXY-2,2
          nodalWeights(node_idx)  = weight
        END DO
        ! left/right
        DO node_idx=2,2*numberGlobalYelements-2,2
         nodalWeights(node_idx*numNodesX+1)   = weight
         nodalWeights((node_idx+1)*numNodesX) = weight
        END DO
        !    __ __ __ __ __ __
        !   |     |     |     |
        !
        !   |     |     |     |
        !    __ __X__ __X__ __
        !   |     |     |     |
        !
        !   |     |     |     |
        !    __ __X__ __X__ __
        !   |     |     |     |
        !
        !   |     |     |     |
        !    __ __ __ __ __ __
        !
        weight                                = 4.0_DP / 36.0_DP
        DO node_idx=2,2*numberGlobalYelements-2,2
          DO node_idx_2=node_idx*numNodesX+3,(node_idx+1)*numNodesX-2,2
            nodalWeights(node_idx_2)  = weight
          END DO
        END DO
        !    __ __ __ __ __ __
        !   |     |     |     |
        !      X     X     X
        !   |     |     |     |
        !    __ __ __ __ __ __
        !   |     |     |     |
        !      X     X     X
        !   |     |     |     |
        !    __ __ __ __ __ __
        !   |     |     |     |
        !      X     X     X
        !   |     |     |     |
        !    __ __ __ __ __ __
        !
        weight                                = 4.0_DP / 9.0_DP
        DO node_idx=1,2*numberGlobalYelements-1,2
          DO node_idx_2=node_idx*numNodesX+2,(node_idx+1)*numNodesX-1,2
            nodalWeights(node_idx_2)  = weight
          END DO
        END DO
        !
        !    __X__ __X__ __X__
        !   |     |     |     |
        !   X                 X
        !   |     |     |     |
        !    __ __ __ __ __ __
        !   |     |     |     |
        !   X                 X
        !   |     |     |     |
        !    __ __ __ __ __ __
        !   |     |     |     |
        !   X                 X
        !   |     |     |     |
        !    __X__ __X__ __X__
        !
        weight                                = 1.0_DP / 9.0_DP
        ! top/bottom
        DO node_idx=2,numNodesX-1,2
          nodalWeights(node_idx)                      = weight
          nodalWeights(numNodesXY-numNodesX+node_idx) = weight
        END DO
        ! left/right
        DO node_idx=1,2*numberGlobalYelements-1,2
          nodalWeights(node_idx*numNodesX+1)    = weight
          nodalWeights((node_idx+1)*numNodesX)  = weight
        END DO
        !    __ __ __ __ __ __
        !   |     |     |     |
        !         X     X
        !   |     |     |     |
        !    __X__ __X__ __X__
        !   |     |     |     |
        !         X     X
        !   |     |     |     |
        !    __X__ __X__ __X__
        !   |     |     |     |
        !         X     X
        !   |     |     |     |
        !    __X__ __X__ __X__
        !
        weight                                = 2.0_DP / 9.0_DP
        ! vertical edges
        DO node_idx=1,2*numberGlobalYelements-1,2
          DO node_idx_2=node_idx*numNodesX+3,(node_idx+1)*numNodesX-2,2
            nodalWeights(node_idx_2)  = weight
          END DO
        END DO
        ! horizontal edges
        DO node_idx=2,2*numberGlobalYelements-2,2
          DO node_idx_2=node_idx*numNodesX+2,(node_idx+1)*numNodesX-1,2
            nodalWeights(node_idx_2)  = weight
          END DO
        END DO
        ! we are done.. do a sanity check
        IF(ABS(SUM(nodalWeights)-numberGlobalXelements*numberGlobalYelements)>1.0e-11) THEN
          WRITE(*,*) ""
          WRITE(*,*) "Invalid nodal weights for computing consistent nodal forces."
          WRITE(*,*) " SIZE(nodalWeights),numNodesX,numNodesY,SUM(nodalWeights) = ", &
            & SIZE(nodalWeights),numNodesX,numNodesY,SUM(nodalWeights)
          WRITE(*,*) ""
          DO node_idx=1,numNodesY
            WRITE(*,*) " ",nodalWeights(1+(node_idx-1)*numNodesX:node_idx*numNodesX:1)
          END DO
          WRITE(*,*) ""
          STOP
        END IF
        nodalWeights=nodalWeights/SUM(nodalWeights)
      CASE(GENERATED_MESH_REGULAR_FRONT_SURFACE,GENERATED_MESH_REGULAR_BACK_SURFACE)
        ! make sure that there are some nodes on this surface
        IF((numberGlobalXElements<=0).OR.(numberGlobalZElements<=0)) THEN
          WRITE(*,*) ">>>GeneratedMesh_SurfaceWeightsGet: There are no nodes on this surface. Face number: ",faceNumber
          STOP
        END IF
        ! allocate and initialise nodal weights
        ALLOCATE(nodalWeights(numNodesXZ))
        nodalWeights  = 0.0_DP
        !
        !   X__ __ __ __ __ __X
        !   |     |     |     |
        !
        !   |     |     |     |
        !    __ __ __ __ __ __
        !   |     |     |     |
        !
        !   |     |     |     |
        !    __ __ __ __ __ __
        !   |     |     |     |

        !   |     |     |     |
        !   X__ __ __ __ __ __X
        !
        weight                                = 1.0_DP / 36.0_DP
        nodalWeights(1)                       = weight
        nodalWeights(numNodesX)               = weight
        nodalWeights(numNodesXZ-numNodesX+1)  = weight
        nodalWeights(numNodesXZ)              = weight
        !
        !    __ __X__ __X__ __
        !   |     |     |     |
        !
        !   |     |     |     |
        !   X__ __ __ __ __ __X
        !   |     |     |     |
        !
        !   |     |     |     |
        !   X__ __ __ __ __ __X
        !   |     |     |     |
        !
        !   |     |     |     |
        !    __ __X__ __X__ __
        !
        weight                                = 2.0_DP / 36.0_DP
        ! bottom
        DO node_idx=3,numNodesX-2,2
          nodalWeights(node_idx)  = weight
        END DO
        ! top
        DO node_idx=numNodesXZ-numNodesX+3,numNodesXZ-2,2
          nodalWeights(node_idx)  = weight
        END DO
        ! left/right
        DO node_idx=2,2*numberGlobalZelements-2,2
         nodalWeights(node_idx*numNodesX+1)   = weight
         nodalWeights((node_idx+1)*numNodesX) = weight
        END DO
        !    __ __ __ __ __ __
        !   |     |     |     |
        !
        !   |     |     |     |
        !    __ __X__ __X__ __
        !   |     |     |     |
        !
        !   |     |     |     |
        !    __ __X__ __X__ __
        !   |     |     |     |
        !
        !   |     |     |     |
        !    __ __ __ __ __ __
        !
        weight                                = 4.0_DP / 36.0_DP
        DO node_idx=2,2*numberGlobalZelements-2,2
          DO node_idx_2=node_idx*numNodesX+3,(node_idx+1)*numNodesX-2,2
            nodalWeights(node_idx_2)  = weight
          END DO
        END DO
        !    __ __ __ __ __ __
        !   |     |     |     |
        !      X     X     X
        !   |     |     |     |
        !    __ __ __ __ __ __
        !   |     |     |     |
        !      X     X     X
        !   |     |     |     |
        !    __ __ __ __ __ __
        !   |     |     |     |
        !      X     X     X
        !   |     |     |     |
        !    __ __ __ __ __ __
        !
        weight                                = 4.0_DP / 9.0_DP
        DO node_idx=1,2*numberGlobalZelements-1,2
          DO node_idx_2=node_idx*numNodesX+2,(node_idx+1)*numNodesX-1,2
            nodalWeights(node_idx_2)  = weight
          END DO
        END DO
        !
        !    __X__ __X__ __X__
        !   |     |     |     |
        !   X                 X
        !   |     |     |     |
        !    __ __ __ __ __ __
        !   |     |     |     |
        !   X                 X
        !   |     |     |     |
        !    __ __ __ __ __ __
        !   |     |     |     |
        !   X                 X
        !   |     |     |     |
        !    __X__ __X__ __X__
        !
        weight                                = 1.0_DP / 9.0_DP
        ! top/bottom
        DO node_idx=2,numNodesX-1,2
          nodalWeights(node_idx)                      = weight
          nodalWeights(numNodesXZ-numNodesX+node_idx) = weight
        END DO
        ! left/right
        DO node_idx=1,2*numberGlobalZelements-1,2
          nodalWeights(node_idx*numNodesX+1)    = weight
          nodalWeights((node_idx+1)*numNodesX)  = weight
        END DO
        !    __ __ __ __ __ __
        !   |     |     |     |
        !         X     X
        !   |     |     |     |
        !    __X__ __X__ __X__
        !   |     |     |     |
        !         X     X
        !   |     |     |     |
        !    __X__ __X__ __X__
        !   |     |     |     |
        !         X     X
        !   |     |     |     |
        !    __X__ __X__ __X__
        !
        weight                                = 2.0_DP / 9.0_DP
        ! vertical edges
        DO node_idx=1,2*numberGlobalZelements-1,2
          DO node_idx_2=node_idx*numNodesX+3,(node_idx+1)*numNodesX-2,2
            nodalWeights(node_idx_2)  = weight
          END DO
        END DO
        ! horizontal edges
        DO node_idx=2,2*numberGlobalZelements-2,2
          DO node_idx_2=node_idx*numNodesX+2,(node_idx+1)*numNodesX-1,2
            nodalWeights(node_idx_2)  = weight
          END DO
        END DO
        ! we are done.. do a sanity check
        IF(ABS(SUM(nodalWeights)-numberGlobalXelements*numberGlobalZelements)>1.0e-11) THEN
          WRITE(*,*) ""
          WRITE(*,*) "Invalid nodal weights for computing consistent nodal forces."
          WRITE(*,*) " SIZE(nodalWeights),numNodesX,numNodesY,SUM(nodalWeights) = ", &
            & SIZE(nodalWeights),numNodesX,numNodesZ,SUM(nodalWeights)
          WRITE(*,*) ""
          DO node_idx=1,numNodesZ
            WRITE(*,*) " ",nodalWeights(1+(node_idx-1)*numNodesX:node_idx*numNodesX:1)
          END DO
          WRITE(*,*) ""
          STOP
        END IF
        nodalWeights=nodalWeights/SUM(nodalWeights)
      CASE DEFAULT
        WRITE(*,*) "Not implemented: Face number ",faceNumber
        STOP
    END SELECT

  END SUBROUTINE GeneratedMesh_SurfaceWeightsGet

  !
  !=============================================================================
  !

  !>Returns starting and stopping index of nodes belonging to a surface of given patch ID
  SUBROUTINE ImportedMesh_SurfaceGet(BoundaryPatches,PatchID,StartIdx,StopIdx,Err)

    !Argument variables
    INTEGER(INTG),      INTENT(IN)  :: BoundaryPatches(:)   !< The boundary patch labels for all boundary nodes
    INTEGER(INTG),      INTENT(IN)  :: PatchID              !< The desired boundary patch label
    INTEGER(INTG),      INTENT(OUT) :: StartIdx             !< On return, first index for corresponding PatchID
    INTEGER(INTG),      INTENT(OUT) :: StopIdx              !< On return, last index for corresponding PatchID
    INTEGER(INTG),      INTENT(OUT) :: Err                  !< The error code.
    !Local variables
    LOGICAL                         :: BoundaryFound
    INTEGER(INTG)                   :: NodeIdx,NumberOfPatchIDs


    BoundaryFound=.FALSE.
    NumberOfPatchIDs=BoundaryPatches(1)
    ! Check for minimum length of BoundaryPatches variables
    IF(SIZE(BoundaryPatches)<1+3*NumberOfPatchIDs) THEN
      WRITE(*,*) "BoundaryPatches variable has incompatible size."
      STOP
    END IF
    ! Set index to first non-header index
    StartIdx=2+NumberOfPatchIDs*2
    StopIdx=1+NumberOfPatchIDs*2
    ! Set StartIdx and StopIdx for given PatchID
    DO NodeIdx=2,NumberOfPatchIDs+1
      StopIdx=StopIdx+BoundaryPatches(NodeIdx)
      IF(BoundaryPatches(NodeIdx+NumberOfPatchIDs)==PatchID) THEN
        BoundaryFound=.TRUE.
        EXIT
      ELSE
        StartIdx=StartIdx+BoundaryPatches(NodeIdx)
      END IF
    END DO
    IF(.NOT.BoundaryFound) THEN
      WRITE(*,*) "Could not find boundary patch ID."
      STOP
    END IF
    RETURN
  END SUBROUTINE ImportedMesh_SurfaceGet

END MODULE meshReader
