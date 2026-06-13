#include "Geo_Element.h"
//#include "Constants.h"
#include "Glb_Dimension.h"
using namespace std;

namespace PHSPACE
{

BasicElement::BasicElement()
{
    
}

BasicElement::~BasicElement()
{
}

int BasicElement::GetElementNodeNumber(int elem_type)
{
    int node_number = 0;
    if (elem_type == NODE)
    {
        node_number = 1;
    }
    else if (elem_type == BAR_2)
    {
        node_number = 2;
    }
    else if (elem_type == BAR_3)
    {
        node_number = 3;
    }
    else if (elem_type == TRI_3)
    {
        node_number = 3;
    }
    else if (elem_type == TRI_6)
    {
        node_number = 6;
    }
    else if (elem_type == QUAD_4)
    {
        node_number = 4;
    }
    else if (elem_type == QUAD_6)
    {
        node_number = 6;
    }
    else if (elem_type == QUAD_8)
    {
        node_number = 8;
    }
    else if (elem_type == QUAD_9)
    {
        node_number = 9;
    }
    else if (elem_type == TETRA_4)
    {
        node_number = 4;
    }
    else if (elem_type == TETRA_10)
    {
        node_number = 10;
    }
    else if (elem_type == PYRA_5)
    {
        node_number = 5;
    }
    else if (elem_type == PYRA_14)
    {
        node_number = 14;
    }
    else if (elem_type == PENTA_6)
    {
        node_number = 6;
    }
    else if (elem_type == PENTA_12)
    {
        node_number = 12;
    }
    else if (elem_type == PENTA_15)
    {
        node_number = 15;
    }
    else if (elem_type == PENTA_18)
    {
        node_number = 18;
    }
    else if (elem_type == HEXA_8)
    {
        node_number = 8;
    }
    else if (elem_type == HEXA_20)
    {
        node_number = 20;
    }
    else if (elem_type == HEXA_27)
    {
        node_number = 27;
    }

    return node_number;
}

void BasicElement::Init(int elem_type)
{
    this->elem_type = elem_type;

    int node_max = GetElementNodeNumber(elem_type);
    mp_struct.resize(node_max);

    if (elem_type == NODE)
    {
        //Face 1
        // 1

        face_list.resize(1);
        face_list[0].push_back(1-1);
    }
    else if (elem_type == BAR_2)
    {
        //Face 1
        // 1 2

        face_list.resize(1);
        face_list[0].push_back(1-1);
        face_list[0].push_back(2-1);
    }
    else if (elem_type == BAR_3)
    {
        //Face 2
        // 1 3
        // 3 2

        face_element_type.resize(2);
        face_list.resize(2);
        face_element_type[0] = BAR_2;
        face_list[0].push_back(1-1);
        face_list[0].push_back(3-1);

        face_element_type[1] = BAR_2;
        face_list[1].push_back(3-1);
        face_list[1].push_back(2-1);


        //! BAR_3
        //! 1.........3..........2
        //! For the midpoint 3, the adjacent points are 1and 2,change them to subscript which starts from zero is (0,1).
        int_array &mp = mp_struct[2];
        mp.push_back(0);
        mp.push_back(1);

        //! subcell
        //! BAR_2 1 3
        //! BAR_2 3 2

        child_element_type.resize(2);
        child_element_index.resize(2);

        child_element_type[0] = BAR_2;
        child_element_index[0].push_back(1-1);
        child_element_index[0].push_back(3-1);

        child_element_type[1] = BAR_2;
        child_element_index[1].push_back(3-1);
        child_element_index[1].push_back(2-1);
    }
    else if (elem_type == TRI_3)
    {
        //! TRI_3
        //! 1.........2.........3

        //! Face 3
        //! 1 2
        //! 2 3
        //! 3 1
        face_element_type.resize(3);
        face_list.resize(3);

        face_element_type[0] = BAR_2;
        face_list[0].push_back(1-1);
        face_list[0].push_back(2-1);

        face_element_type[1] = BAR_2;
        face_list[1].push_back(2-1);
        face_list[1].push_back(3-1);

        face_element_type[2] = BAR_2;
        face_list[2].push_back(3-1);
        face_list[2].push_back(1-1);
    }
    else if (elem_type == TRI_6)
    {
        //! TRI_6
        //! 1....4....2 ...5....3....6....1

        //! Composite Face 3
        //! BAR_3 1 2 4
        //! BAR_3 2 3 5
        //! BAR_3 3 1 6
        composite_face_element_type.resize(3);
        composite_face_list.resize(3);

        composite_face_element_type[0] = BAR_3;
        composite_face_list[0].push_back(1-1);
        composite_face_list[0].push_back(2-1);
        composite_face_list[0].push_back(4-1);

        composite_face_element_type[1] = BAR_3;
        composite_face_list[1].push_back(2-1);
        composite_face_list[1].push_back(3-1);
        composite_face_list[1].push_back(5-1);

        composite_face_element_type[2] = BAR_3;
        composite_face_list[2].push_back(3-1);
        composite_face_list[2].push_back(1-1);
        composite_face_list[2].push_back(6-1);

        face_element_type.resize(6);
        face_list.resize(6);
        face_element_type[0] = BAR_2;
        face_list[0].push_back(1-1);
        face_list[0].push_back(4-1);

        face_element_type[1] = BAR_2;
        face_list[1].push_back(4-1);
        face_list[1].push_back(2-1);

        face_element_type[2] = BAR_2;
        face_list[2].push_back(2-1);
        face_list[2].push_back(5-1);

        face_element_type[3] = BAR_2;
        face_list[3].push_back(5-1);
        face_list[3].push_back(3-1);

        face_element_type[4] = BAR_2;
        face_list[4].push_back(3-1);
        face_list[4].push_back(6-1);

        face_element_type[5] = BAR_2;
        face_list[5].push_back(6-1);
        face_list[5].push_back(1-1);


        //! Face 6
        //! 1 4
        //! 4 2
        //! 2 5
        //! 5 3
        //! 3 6
        //! 6 1
        face_element_type.resize(6);
        face_list.resize(6);
        face_element_type[0] = BAR_2;
        face_list[0].push_back(1-1);
        face_list[0].push_back(4-1);

        face_element_type[1] = BAR_2;
        face_list[1].push_back(4-1);
        face_list[1].push_back(2-1);

        face_element_type[2] = BAR_2;
        face_list[2].push_back(2-1);
        face_list[2].push_back(5-1);

        face_element_type[3] = BAR_2;
        face_list[3].push_back(5-1);
        face_list[3].push_back(3-1);

        face_element_type[4] = BAR_2;
        face_list[4].push_back(3-1);
        face_list[4].push_back(6-1);

        face_element_type[5] = BAR_2;
        face_list[5].push_back(6-1);
        face_list[5].push_back(1-1);

        //! TRI_6
        //! 1.........4..........2
        //! 2.........5..........3
        //! 3.........6..........1

        mp_struct[3].push_back(0);
        mp_struct[3].push_back(1);

        mp_struct[4].push_back(1);
        mp_struct[4].push_back(2);

        mp_struct[5].push_back(2);
        mp_struct[5].push_back(0);

        //! subcell
        //! TRI_3 1 4 6
        //! TRI_3 4 2 5
        //! TRI_3 4 5 6
        //! TRI_3 6 5 3 

        child_element_type.resize(4);
        child_element_index.resize(4);

        child_element_type[0] = TRI_3;
        child_element_index[0].push_back(1-1);
        child_element_index[0].push_back(4-1);
        child_element_index[0].push_back(6-1);

        child_element_type[1] = TRI_3;
        child_element_index[1].push_back(4-1);
        child_element_index[1].push_back(2-1);
        child_element_index[1].push_back(5-1);

        child_element_type[2] = TRI_3;
        child_element_index[2].push_back(4-1);
        child_element_index[2].push_back(5-1);
        child_element_index[2].push_back(6-1);

        child_element_type[3] = TRI_3;
        child_element_index[3].push_back(6-1);
        child_element_index[3].push_back(5-1);
        child_element_index[3].push_back(3-1);
    }
    else if (elem_type == QUAD_4)
    {
        //! QUAD_4
        //! 1....2....3 ...4....1

        //! Face 4
        //! 1 2
        //! 2 3
        //! 3 4
        //! 4 1
        face_element_type.resize(4);
        face_list.resize(4);

        face_element_type[0] = BAR_2;
        face_list[0].push_back(1-1);
        face_list[0].push_back(2-1);

        face_element_type[1] = BAR_2;
        face_list[1].push_back(2-1);
        face_list[1].push_back(3-1);

        face_element_type[2] = BAR_2;
        face_list[2].push_back(3-1);
        face_list[2].push_back(4-1);

        face_element_type[3] = BAR_2;
        face_list[3].push_back(4-1);
        face_list[3].push_back(1-1);
    }
    else if (elem_type == QUAD_6)
    {
        //! QUAD_6
        //! Composite Face 4
        composite_face_element_type.resize(2);
        composite_face_list.resize(2);

        composite_face_element_type[ 0 ] = BAR_3;
        composite_face_list[ 0 ].push_back(1 - 1);
        composite_face_list[ 0 ].push_back(2 - 1);
        composite_face_list[ 0 ].push_back(5 - 1);

        composite_face_element_type[ 1 ] = BAR_3;
        composite_face_list[ 1 ].push_back(3 - 1);
        composite_face_list[ 1 ].push_back(4 - 1);
        composite_face_list[ 1 ].push_back(6 - 1);

        //! Face
        face_element_type.resize(6);
        face_list.resize(6);

        face_element_type[ 0 ] = BAR_2;
        face_list[ 0 ].push_back(1 - 1);
        face_list[ 0 ].push_back(5 - 1);

        face_element_type[ 1 ] = BAR_2;
        face_list[ 1 ].push_back(5 - 1);
        face_list[ 1 ].push_back(2 - 1);

        face_element_type[ 2 ] = BAR_2;
        face_list[ 2 ].push_back(2 - 1);
        face_list[ 2 ].push_back(3 - 1);

        face_element_type[ 3 ] = BAR_2;
        face_list[ 3 ].push_back(3 - 1);
        face_list[ 3 ].push_back(6 - 1);

        face_element_type[ 4 ] = BAR_2;
        face_list[ 4 ].push_back(6 - 1);
        face_list[ 4 ].push_back(4 - 1);

        face_element_type[ 5 ] = BAR_2;
        face_list[ 5 ].push_back(4 - 1);
        face_list[ 5 ].push_back(1 - 1);


        //! QUAD_6
        mp_struct[ 4 ].push_back(0);
        mp_struct[ 4 ].push_back(1);

        mp_struct[ 5 ].push_back(2);
        mp_struct[ 5 ].push_back(3);

        //! subcell
        child_element_type.resize(2);
        child_element_index.resize(2);

        child_element_type[ 0 ] = QUAD_4;
        child_element_index[ 0 ].push_back(1 - 1);
        child_element_index[ 0 ].push_back(5 - 1);
        child_element_index[ 0 ].push_back(6 - 1);
        child_element_index[ 0 ].push_back(4 - 1);

        child_element_type[ 1 ] = QUAD_4;
        child_element_index[ 1 ].push_back(5 - 1);
        child_element_index[ 1 ].push_back(2 - 1);
        child_element_index[ 1 ].push_back(3 - 1);
        child_element_index[ 1 ].push_back(6 - 1);
    }
    else if (elem_type == QUAD_8)
    {
        //! QUAD_8
        //! 1..5.2..6..3..7..4..8..1

        //! Composite Face 4
        //! BAR_3 1 2 5
        //! BAR_3 2 3 6
        //! BAR_3 3 4 7
        //! BAR_3 4 1 8
        composite_face_element_type.resize(4);
        composite_face_list.resize(4);

        composite_face_element_type[0] = BAR_3;
        composite_face_list[0].push_back(1-1);
        composite_face_list[0].push_back(2-1);
        composite_face_list[0].push_back(5-1);

        composite_face_element_type[1] = BAR_3;
        composite_face_list[1].push_back(2-1);
        composite_face_list[1].push_back(3-1);
        composite_face_list[1].push_back(6-1);

        composite_face_element_type[2] = BAR_3;
        composite_face_list[2].push_back(3-1);
        composite_face_list[2].push_back(4-1);
        composite_face_list[2].push_back(7-1);

        composite_face_element_type[3] = BAR_3;
        composite_face_list[3].push_back(4-1);
        composite_face_list[3].push_back(1-1);
        composite_face_list[3].push_back(8-1);

        //! Face 8
        //! 1 5
        //! 5 2
        //! 2 6
        //! 6 3

        //! 3 7
        //! 7 4
        //! 4 8
        //! 8 1

        face_element_type.resize(8);
        face_list.resize(8);

        face_element_type[0] = BAR_2;
        face_list[0].push_back(1-1);
        face_list[0].push_back(5-1);

        face_element_type[1] = BAR_2;
        face_list[1].push_back(5-1);
        face_list[1].push_back(2-1);

        face_element_type[2] = BAR_2;
        face_list[2].push_back(2-1);
        face_list[2].push_back(6-1);

        face_element_type[3] = BAR_2;
        face_list[3].push_back(6-1);
        face_list[3].push_back(3-1);

        face_element_type[4] = BAR_2;
        face_list[4].push_back(3-1);
        face_list[4].push_back(7-1);

        face_element_type[5] = BAR_2;
        face_list[5].push_back(7-1);
        face_list[5].push_back(4-1);

        face_element_type[6] = BAR_2;
        face_list[6].push_back(4-1);
        face_list[6].push_back(8-1);

        face_element_type[7] = BAR_2;
        face_list[7].push_back(8-1);
        face_list[7].push_back(1-1);

        //! QUAD_8
        //! 1.........5..........2
        //! 2.........6..........3
        //! 3.........7..........4
        //! 4.........8..........1

        mp_struct[4].push_back(0);
        mp_struct[4].push_back(1);

        mp_struct[5].push_back(1);
        mp_struct[5].push_back(2);

        mp_struct[6].push_back(2);
        mp_struct[6].push_back(3);

        mp_struct[7].push_back(3);
        mp_struct[7].push_back(0);
    }
    else if (elem_type == QUAD_9)
    {
        //! QUAD_9
        //! 1..5.2..6..3..7..4..8..1

        //! Composite Face 4
        //! BAR_3 1 2 5
        //! BAR_3 2 3 6
        //! BAR_3 3 4 7
        //! BAR_3 4 1 8
        composite_face_element_type.resize(4);
        composite_face_list.resize(4);

        composite_face_element_type[0] = BAR_3;
        composite_face_list[0].push_back(1-1);
        composite_face_list[0].push_back(2-1);
        composite_face_list[0].push_back(5-1);

        composite_face_element_type[1] = BAR_3;
        composite_face_list[1].push_back(2-1);
        composite_face_list[1].push_back(3-1);
        composite_face_list[1].push_back(6-1);

        composite_face_element_type[2] = BAR_3;
        composite_face_list[2].push_back(3-1);
        composite_face_list[2].push_back(4-1);
        composite_face_list[2].push_back(7-1);

        composite_face_element_type[3] = BAR_3;
        composite_face_list[3].push_back(4-1);
        composite_face_list[3].push_back(1-1);
        composite_face_list[3].push_back(8-1);

        //! Face 8
        //! 1 5
        //! 5 2
        //! 2 6
        //! 6 3

        //! 3 7
        //! 7 4
        //! 4 8
        //! 8 1
        face_element_type.resize(8);
        face_list.resize(8);

        face_element_type[0] = BAR_2;
        face_list[0].push_back(1-1);
        face_list[0].push_back(5-1);

        face_element_type[1] = BAR_2;
        face_list[1].push_back(5-1);
        face_list[1].push_back(2-1);

        face_element_type[2] = BAR_2;
        face_list[2].push_back(2-1);
        face_list[2].push_back(6-1);

        face_element_type[3] = BAR_2;
        face_list[3].push_back(6-1);
        face_list[3].push_back(3-1);

        face_element_type[4] = BAR_2;
        face_list[4].push_back(3-1);
        face_list[4].push_back(7-1);

        face_element_type[5] = BAR_2;
        face_list[5].push_back(7-1);
        face_list[5].push_back(4-1);

        face_element_type[6] = BAR_2;
        face_list[6].push_back(4-1);
        face_list[6].push_back(8-1);

        face_element_type[7] = BAR_2;
        face_list[7].push_back(8-1);
        face_list[7].push_back(1-1);

        //!QUAD_9
        //! 1.........5..........2
        //! 2.........6..........3
        //! 3.........7..........4
        //! 4.........8..........1
        //! 1 2.........9..........3 4

        mp_struct[4].push_back(0);
        mp_struct[4].push_back(1);

        mp_struct[5].push_back(1);
        mp_struct[5].push_back(2);

        mp_struct[6].push_back(2);
        mp_struct[6].push_back(3);

        mp_struct[7].push_back(3);
        mp_struct[7].push_back(0);

        mp_struct[8].push_back(0);
        mp_struct[8].push_back(1);
        mp_struct[8].push_back(2);
        mp_struct[8].push_back(3);

        //! subcell
        //! QUAD_4 : 1 5 9 8
        //! QUAD_4 : 5 2 6 9
        //! QUAD_4 : 9 6 3 7
        //! QUAD_4 : 8 9 7 4

        child_element_type.resize(4);
        child_element_index.resize(4);

        child_element_type[0] = QUAD_4;
        child_element_index[0].push_back(1-1);
        child_element_index[0].push_back(5-1);
        child_element_index[0].push_back(9-1);
        child_element_index[0].push_back(8-1);

        child_element_type[1] = QUAD_4;
        child_element_index[1].push_back(5-1);
        child_element_index[1].push_back(2-1);
        child_element_index[1].push_back(6-1);
        child_element_index[1].push_back(9-1);

        child_element_type[2] = QUAD_4;
        child_element_index[2].push_back(9-1);
        child_element_index[2].push_back(6-1);
        child_element_index[2].push_back(3-1);
        child_element_index[2].push_back(7-1);

        child_element_type[3] = QUAD_4;
        child_element_index[3].push_back(8-1);
        child_element_index[3].push_back(9-1);
        child_element_index[3].push_back(7-1);
        child_element_index[3].push_back(4-1);
    }
    else if (elem_type == TETRA_4)
    {
        //! TETRA_4
        //! 1..2...3...4

        //! Face 4
        //! 1 3 2
        //! 1 2 4
        //! 2 3 4
        //! 3 1 4

        face_element_type.resize(4);
        face_list.resize(4);

        face_element_type[0] = TRI_3;
        face_list[0].push_back(1-1);
        face_list[0].push_back(3-1);
        face_list[0].push_back(2-1);

        face_element_type[1] = TRI_3;
        face_list[1].push_back(1-1);
        face_list[1].push_back(2-1);
        face_list[1].push_back(4-1);

        face_element_type[2] = TRI_3;
        face_list[2].push_back(2-1);
        face_list[2].push_back(3-1);
        face_list[2].push_back(4-1);

        face_element_type[3] = TRI_3;
        face_list[3].push_back(3-1);
        face_list[3].push_back(1-1);
        face_list[3].push_back(4-1);
    }
    else if (elem_type == TETRA_10)
    {
        //! TETRA_10
        //! 1..7..3..6..2..5..1

        //! Composite Face 4
        //! TRI_6 1  3  2  7  6  5
        //! TRI_6 2  3  4  6  10 9
        //! TRI_6 3  1  4  7  8  10
        //! TRI_6 1  2  4  5  9  8
        composite_face_element_type.resize(4);
        composite_face_list.resize(4);

        composite_face_element_type[0] = TRI_6;
        composite_face_list[0].push_back(1-1);
        composite_face_list[0].push_back(3-1);
        composite_face_list[0].push_back(2-1);
        composite_face_list[0].push_back(7-1);
        composite_face_list[0].push_back(6-1);
        composite_face_list[0].push_back(5-1);

        composite_face_element_type[1] = TRI_6;
        composite_face_list[1].push_back(2-1);
        composite_face_list[1].push_back(3-1);
        composite_face_list[1].push_back(4-1);
        composite_face_list[1].push_back(6-1);
        composite_face_list[1].push_back(10-1);
        composite_face_list[1].push_back(9-1);

        composite_face_element_type[2] = TRI_6;
        composite_face_list[2].push_back(3-1);
        composite_face_list[2].push_back(1-1);
        composite_face_list[2].push_back(4-1);
        composite_face_list[2].push_back(7-1);
        composite_face_list[2].push_back(8-1);
        composite_face_list[2].push_back(10-1);

        composite_face_element_type[3] = TRI_6;
        composite_face_list[3].push_back(1-1);
        composite_face_list[3].push_back(2-1);
        composite_face_list[3].push_back(4-1);
        composite_face_list[3].push_back(5-1);
        composite_face_list[3].push_back(9-1);
        composite_face_list[3].push_back(8-1);

        //! Face 16
        //! 1  7  5
        //! 7  6  5
        //! 5  6  2
        //! 6  7  3

        face_element_type.resize(16);
        face_list.resize(16);

        face_element_type[0]= TRI_3;
        face_list[0].push_back(1-1);
        face_list[0].push_back(7-1);
        face_list[0].push_back(5-1);

        face_element_type[1]= TRI_3;
        face_list[1].push_back(7-1);
        face_list[1].push_back(6-1);
        face_list[1].push_back(5-1);

        face_element_type[2]= TRI_3;
        face_list[2].push_back(5-1);
        face_list[2].push_back(6-1);
        face_list[2].push_back(2-1);

        face_element_type[3]= TRI_3;
        face_list[3].push_back(6-1);
        face_list[3].push_back(7-1);
        face_list[3].push_back(3-1);

        //! 2..6..3..10..4..9..2
        //! 2  6  9
        //! 9  6  10
        //! 10 6  3
        //! 10 4  9

        face_element_type[4]= TRI_3;
        face_list[4].push_back(2-1);
        face_list[4].push_back(6-1);
        face_list[4].push_back(9-1);

        face_element_type[5]= TRI_3;
        face_list[5].push_back(9 -1);
        face_list[5].push_back(6 -1);
        face_list[5].push_back(10-1);

        face_element_type[6]= TRI_3;
        face_list[6].push_back(10-1);
        face_list[6].push_back(6 -1);
        face_list[6].push_back(3 -1);

        face_element_type[7]= TRI_3;
        face_list[7].push_back(10-1);
        face_list[7].push_back(4 -1);
        face_list[7].push_back(9 -1);

        //! 1..8..4..10..3..7..1
        //! 1  8  7
        //! 8  10 7
        //! 10 3  7
        //! 10 8  4

        face_element_type[8]= TRI_3;
        face_list[8].push_back(1-1);
        face_list[8].push_back(8-1);
        face_list[8].push_back(7-1);

        face_element_type[9]= TRI_3;
        face_list[9].push_back(8 -1);
        face_list[9].push_back(10-1);
        face_list[9].push_back(7 -1);

        face_element_type[10]= TRI_3;
        face_list[10].push_back(10-1);
        face_list[10].push_back(3 -1);
        face_list[10].push_back(7 -1);

        face_element_type[11]= TRI_3;
        face_list[11].push_back(10-1);
        face_list[11].push_back(8 -1);
        face_list[11].push_back(4 -1);

        //! 1..5..2..9..4..8..1
        //! 1  5  8
        //! 8  5  9
        //! 9  5  2
        //! 9  4  8
        face_element_type[12]= TRI_3;
        face_list[12].push_back(1-1);
        face_list[12].push_back(5-1);
        face_list[12].push_back(8-1);

        face_element_type[13]= TRI_3;
        face_list[13].push_back(8-1);
        face_list[13].push_back(5-1);
        face_list[13].push_back(9-1);

        face_element_type[14]= TRI_3;
        face_list[14].push_back(9-1);
        face_list[14].push_back(5-1);
        face_list[14].push_back(2-1);

        face_element_type[15]= TRI_3;
        face_list[15].push_back(9-1);
        face_list[15].push_back(4-1);
        face_list[15].push_back(8-1);

        //! TETRA_10
        //! 1.........5..........2
        //! 2.........6..........3
        //! 3.........7..........1

        //! 1.........8..........4
        //! 2.........9..........4
        //! 3.........10.........4

        mp_struct[4].push_back(1-1);
        mp_struct[4].push_back(2-1);

        mp_struct[5].push_back(2-1);
        mp_struct[5].push_back(3-1);

        mp_struct[6].push_back(3-1);
        mp_struct[6].push_back(1-1);

        mp_struct[7].push_back(1-1);
        mp_struct[7].push_back(4-1);

        mp_struct[8].push_back(2-1);
        mp_struct[8].push_back(4-1);

        mp_struct[9].push_back(3-1);
        mp_struct[9].push_back(4-1);

        //! subcell
        //! TETRA_4 : 8  5  1  7
        //! TETRA_4 : 8  9  5  7
        //! TETRA_4 : 9  10 6  7
        //! TETRA_4 : 10 3  6  7
        //! TETRA_4 : 2  6  5  9
        //! TETRA_4 : 5  6  7  9
        //! TETRA_4 : 4  8  10 9
        //! TETRA_4 : 10 8  7  9

        child_element_type.resize(8);
        child_element_index.resize(8);

        child_element_type[0] = TETRA_4;
        child_element_index[0].push_back(8-1);
        child_element_index[0].push_back(5-1);
        child_element_index[0].push_back(1-1);
        child_element_index[0].push_back(7-1);

        child_element_type[1] = TETRA_4;
        child_element_index[1].push_back(8-1);
        child_element_index[1].push_back(9-1);
        child_element_index[1].push_back(5-1);
        child_element_index[1].push_back(7-1);

        child_element_type[2] = TETRA_4;
        child_element_index[2].push_back(9-1);
        child_element_index[2].push_back(10-1);
        child_element_index[2].push_back(6-1);
        child_element_index[2].push_back(7-1);

        child_element_type[3] = TETRA_4;
        child_element_index[3].push_back(10-1);
        child_element_index[3].push_back(3-1);
        child_element_index[3].push_back(6-1);
        child_element_index[3].push_back(7-1);

        child_element_type[4] = TETRA_4;
        child_element_index[4].push_back(2-1);
        child_element_index[4].push_back(6-1);
        child_element_index[4].push_back(5-1);
        child_element_index[4].push_back(9-1);

        child_element_type[5] = TETRA_4;
        child_element_index[5].push_back(5-1);
        child_element_index[5].push_back(6-1);
        child_element_index[5].push_back(7-1);
        child_element_index[5].push_back(9-1);

        child_element_type[6] = TETRA_4;
        child_element_index[6].push_back(4-1);
        child_element_index[6].push_back(8-1);
        child_element_index[6].push_back(10-1);
        child_element_index[6].push_back(9-1);

        child_element_type[7] = TETRA_4;
        child_element_index[7].push_back(10-1);
        child_element_index[7].push_back(8-1);
        child_element_index[7].push_back(7-1);
        child_element_index[7].push_back(9-1);

        //! Case2: A second quad-dividing method.
        child_element_index_case2.resize(8);

        child_element_index_case2[ 0 ].push_back(8 - 1);    //! Store the corresponding composite elements index number according to the local index numbering sequence of child elements!!! 
        child_element_index_case2[ 0 ].push_back(5 - 1);
        child_element_index_case2[ 0 ].push_back(1 - 1);
        child_element_index_case2[ 0 ].push_back(7 - 1);////////////

        child_element_index_case2[ 1 ].push_back(9 - 1);
        child_element_index_case2[ 1 ].push_back(8 - 1);
        child_element_index_case2[ 1 ].push_back(10 - 1);
        child_element_index_case2[ 1 ].push_back(5 - 1);

        child_element_index_case2[ 2 ].push_back(10 - 1);
        child_element_index_case2[ 2 ].push_back(5 - 1);
        child_element_index_case2[ 2 ].push_back(8 - 1);
        child_element_index_case2[ 2 ].push_back(7 - 1);

        child_element_index_case2[ 3 ].push_back(10 - 1);
        child_element_index_case2[ 3 ].push_back(3 - 1);
        child_element_index_case2[ 3 ].push_back(6 - 1);
        child_element_index_case2[ 3 ].push_back(7 - 1);////////////

        child_element_index_case2[ 4 ].push_back(2 - 1);
        child_element_index_case2[ 4 ].push_back(6 - 1);
        child_element_index_case2[ 4 ].push_back(5 - 1);
        child_element_index_case2[ 4 ].push_back(9 - 1);////////////

        child_element_index_case2[ 5 ].push_back(9 - 1);
        child_element_index_case2[ 5 ].push_back(6 - 1);
        child_element_index_case2[ 5 ].push_back(5 - 1);
        child_element_index_case2[ 5 ].push_back(10 - 1);

        child_element_index_case2[ 6 ].push_back(4 - 1);
        child_element_index_case2[ 6 ].push_back(8 - 1);
        child_element_index_case2[ 6 ].push_back(10 - 1);
        child_element_index_case2[ 6 ].push_back(9 - 1);////////////

        child_element_index_case2[ 7 ].push_back(10 - 1);
        child_element_index_case2[ 7 ].push_back(5 - 1);
        child_element_index_case2[ 7 ].push_back(7 - 1);
        child_element_index_case2[ 7 ].push_back(6 - 1);
    }
    else if (elem_type == PYRA_5)
    {
        //! PYRA_5
        //! 1..2..3..4..5

        //! Face 5
        //! 1  4  3  2
        //! 1  5  4
        //! 1  2  5
        //! 2  3  5
        //! 3  4  5

        face_element_type.resize(5);
        face_list.resize(5);

        face_element_type[0] = QUAD_4;
        face_list[0].push_back(1-1);
        face_list[0].push_back(4-1);
        face_list[0].push_back(3-1);
        face_list[0].push_back(2-1);

        face_element_type[1] = TRI_3;
        face_list[1].push_back(1-1);
        face_list[1].push_back(5-1);
        face_list[1].push_back(4-1);

        face_element_type[2] = TRI_3;
        face_list[2].push_back(1-1);
        face_list[2].push_back(2-1);
        face_list[2].push_back(5-1);

        face_element_type[3] = TRI_3;
        face_list[3].push_back(2-1);
        face_list[3].push_back(3-1);
        face_list[3].push_back(5-1);

        face_element_type[4] = TRI_3;
        face_list[4].push_back(3-1);
        face_list[4].push_back(4-1);
        face_list[4].push_back(5-1);

    }
    else if (elem_type == PYRA_14)
    {
        //! PYRA_14
        //! Face 20

        //! Composite Face 5
        //! QUAD_9 : 1  4  3  2  9   8   7  6  14
        //! TRI_6  : 1  2  5  6  11  10
        //! TRI_6  : 2  3  5  7  12  11
        //! TRI_6  : 3  4  5  8  13  12
        //! TRI_6  : 4  1  5  9  10  13
        composite_face_element_type.resize(5);
        composite_face_list.resize(5);

        composite_face_element_type[0] = QUAD_9;
        composite_face_list[0].push_back(1-1);
        composite_face_list[0].push_back(4-1);
        composite_face_list[0].push_back(3-1);
        composite_face_list[0].push_back(2-1);
        composite_face_list[0].push_back(9-1);
        composite_face_list[0].push_back(8-1);
        composite_face_list[0].push_back(7-1);
        composite_face_list[0].push_back(6-1);
        composite_face_list[0].push_back(14-1);

        composite_face_element_type[1] = TRI_6;
        composite_face_list[ 1 ].push_back(1 - 1);
        composite_face_list[ 1 ].push_back(2 - 1);
        composite_face_list[ 1 ].push_back(5 - 1);
        composite_face_list[ 1 ].push_back(6 - 1);
        composite_face_list[ 1 ].push_back(11 - 1);
        composite_face_list[ 1 ].push_back(10 - 1);

        composite_face_element_type[2] = TRI_6;
        composite_face_list[2].push_back(2-1);
        composite_face_list[2].push_back(3-1);
        composite_face_list[2].push_back(5-1);
        composite_face_list[2].push_back(7-1);
        composite_face_list[2].push_back(12-1);
        composite_face_list[2].push_back(11-1);

        composite_face_element_type[3] = TRI_6;
        composite_face_list[3].push_back(3-1);
        composite_face_list[3].push_back(4-1);
        composite_face_list[3].push_back(5-1);
        composite_face_list[3].push_back(8-1);
        composite_face_list[3].push_back(13-1);
        composite_face_list[3].push_back(12-1);

        composite_face_element_type[4] = TRI_6;
        composite_face_list[4].push_back(4-1);
        composite_face_list[4].push_back(1-1);
        composite_face_list[4].push_back(5-1);
        composite_face_list[4].push_back(9-1);
        composite_face_list[4].push_back(10-1);
        composite_face_list[4].push_back(13-1);

        //! 1  9  4  8  3  7  2  6  1   middle 14
        //! Face 4
        //! 1  9  14 6
        //! 9  4  8  14
        //! 6  14 7  2
        //! 14 8  3  7

        face_element_type.resize(20);
        face_list.resize(20);

        face_element_type[0] = QUAD_4;
        face_list[0].push_back(1-1);
        face_list[0].push_back(9-1);
        face_list[0].push_back(14-1);
        face_list[0].push_back(6-1);

        face_element_type[1] = QUAD_4;
        face_list[1].push_back(9-1);
        face_list[1].push_back(4-1);
        face_list[1].push_back(8-1);
        face_list[1].push_back(14-1);

        face_element_type[2] = QUAD_4;
        face_list[2].push_back(6-1);
        face_list[2].push_back(14-1);
        face_list[2].push_back(7-1);
        face_list[2].push_back(2-1);

        face_element_type[3] = QUAD_4;
        face_list[3].push_back(14-1);
        face_list[3].push_back(8-1);
        face_list[3].push_back(3-1);
        face_list[3].push_back(7-1);

        //! 1  10  5  13  4  9  1
        //! Face 4
        //! 1  10  9
        //! 9  10  13
        //! 4  9   13
        //! 10 5   13

        face_element_type[4] = TRI_3;
        face_list[4].push_back(1-1);
        face_list[4].push_back(10-1);
        face_list[4].push_back(9-1);

        face_element_type[5] = TRI_3;
        face_list[5].push_back(9-1);
        face_list[5].push_back(10-1);
        face_list[5].push_back(13-1);

        face_element_type[6] = TRI_3;
        face_list[6].push_back(4-1);
        face_list[6].push_back(9-1);
        face_list[6].push_back(13-1);

        face_element_type[7] = TRI_3;
        face_list[7].push_back(10-1);
        face_list[7].push_back(5-1);
        face_list[7].push_back(13-1);

        //! 1  6  2  11  5  10  1
        //! Face 4
        //! 1  6  10
        //! 6  11 10
        //! 6  2  11
        //! 10 11 5

        face_element_type[8] = TRI_3;
        face_list[8].push_back(1-1);
        face_list[8].push_back(6-1);
        face_list[8].push_back(10-1);

        face_element_type[9] = TRI_3;
        face_list[9].push_back(6-1);
        face_list[9].push_back(11-1);
        face_list[9].push_back(10-1);

        face_element_type[10] = TRI_3;
        face_list[10].push_back(6-1);
        face_list[10].push_back(2-1);
        face_list[10].push_back(11-1);

        face_element_type[11] = TRI_3;
        face_list[11].push_back(10-1);
        face_list[11].push_back(11-1);
        face_list[11].push_back(5-1);

        //! 2  7   3  12  5  11  2
        //! Face 4
        //! 2  7  11
        //! 7  12 11
        //! 7  3  12
        //! 11 12 5

        face_element_type[12] = TRI_3;
        face_list[12].push_back(2-1);
        face_list[12].push_back(7-1);
        face_list[12].push_back(11-1);

        face_element_type[13] = TRI_3;
        face_list[13].push_back(7-1);
        face_list[13].push_back(12-1);
        face_list[13].push_back(11-1);

        face_element_type[14] = TRI_3;
        face_list[14].push_back(7-1);
        face_list[14].push_back(3-1);
        face_list[14].push_back(12-1);

        face_element_type[15] = TRI_3;
        face_list[15].push_back(11-1);
        face_list[15].push_back(12-1);
        face_list[15].push_back(5-1);

        //! 3  8  4  13  5  12  3
        //! Face 4
        //! 3  8  12
        //! 8  13 12
        //! 8  4  13
        //! 12 13 5

        face_element_type[16] = TRI_3;
        face_list[16].push_back(3-1);
        face_list[16].push_back(8-1);
        face_list[16].push_back(12-1);

        face_element_type[17] = TRI_3;
        face_list[17].push_back(8-1);
        face_list[17].push_back(13-1);
        face_list[17].push_back(12-1);

        face_element_type[18] = TRI_3;
        face_list[18].push_back(8-1);
        face_list[18].push_back(4-1);
        face_list[18].push_back(13-1);

        face_element_type[19] = TRI_3;
        face_list[19].push_back(12-1);
        face_list[19].push_back(13-1);
        face_list[19].push_back(5-1);

        //! PYRA_14
        //! 1.........6..........2
        //! 2.........7..........3
        //! 3.........8..........4
        //! 4.........9..........1

        //! 1.........10.........5
        //! 2.........11.........5
        //! 3.........12.........5
        //! 4.........13.........5

        //! 1 2.......14.........3 4

        mp_struct[5].push_back(0);
        mp_struct[5].push_back(1);

        mp_struct[6].push_back(1);
        mp_struct[6].push_back(2);

        mp_struct[7].push_back(2);
        mp_struct[7].push_back(3);

        mp_struct[8].push_back(3);
        mp_struct[8].push_back(0);

        mp_struct[9].push_back(0);
        mp_struct[9].push_back(4);

        mp_struct[10].push_back(1);
        mp_struct[10].push_back(4);

        mp_struct[11].push_back(2);
        mp_struct[11].push_back(4);

        mp_struct[12].push_back(3);
        mp_struct[12].push_back(4);

        mp_struct[13].push_back(0);
        mp_struct[13].push_back(3);
        mp_struct[13].push_back(2);
        mp_struct[13].push_back(1);

        //! subcell
        //! PYRA_5  : 10 11 12 13 5
        //! PYRA_5  : 1  6  14 9  10
        //! PYRA_5  : 6  2  7  14 11
        //! PYRA_5  : 4  9  14 8  13
        //! PYRA_5  : 14 7  3  8  12
        //! PYRA_5  : 10 13 12 11 14
        //! TETRA_4 : 6  10 11 14
        //! TETRA_4 : 8  12 13 14
        //! TETRA_4 : 9  13 10 14
        //! TETRA_4 : 7  11 12 14

        child_element_type.resize(10);
        child_element_index.resize(10);

        child_element_type[0] = PYRA_5;
        child_element_index[0].push_back(10-1);
        child_element_index[0].push_back(11-1);
        child_element_index[0].push_back(12-1);
        child_element_index[0].push_back(13-1);
        child_element_index[0].push_back(5-1);

        child_element_type[1] = PYRA_5;
        child_element_index[1].push_back(1-1);
        child_element_index[1].push_back(6-1);
        child_element_index[1].push_back(14-1);
        child_element_index[1].push_back(9-1);
        child_element_index[1].push_back(10-1);

        child_element_type[2] = PYRA_5;
        child_element_index[2].push_back(6-1);
        child_element_index[2].push_back(2-1);
        child_element_index[2].push_back(7-1);
        child_element_index[2].push_back(14-1);
        child_element_index[2].push_back(11-1);

        child_element_type[3] = PYRA_5;
        child_element_index[3].push_back(4-1);
        child_element_index[3].push_back(9-1);
        child_element_index[3].push_back(14-1);
        child_element_index[3].push_back(8-1);
        child_element_index[3].push_back(13-1);

        child_element_type[4] = PYRA_5;
        child_element_index[4].push_back(14-1);
        child_element_index[4].push_back(7-1);
        child_element_index[4].push_back(3-1);
        child_element_index[4].push_back(8-1);
        child_element_index[4].push_back(12-1);

        child_element_type[5] = PYRA_5;
        child_element_index[5].push_back(10-1);
        child_element_index[5].push_back(13-1);
        child_element_index[5].push_back(12-1);
        child_element_index[5].push_back(11-1);
        child_element_index[5].push_back(14-1);

        child_element_type[6] = TETRA_4;
        child_element_index[6].push_back(6-1);
        child_element_index[6].push_back(10-1);
        child_element_index[6].push_back(11-1);
        child_element_index[6].push_back(14-1);

        child_element_type[7] = TETRA_4;
        child_element_index[7].push_back(8-1);
        child_element_index[7].push_back(12-1);
        child_element_index[7].push_back(13-1);
        child_element_index[7].push_back(14-1);

        child_element_type[8] = TETRA_4;
        child_element_index[8].push_back(9-1);
        child_element_index[8].push_back(13-1);
        child_element_index[8].push_back(10-1);
        child_element_index[8].push_back(14-1);

        child_element_type[9] = TETRA_4;
        child_element_index[9].push_back(7-1);
        child_element_index[9].push_back(11-1);
        child_element_index[9].push_back(12-1);
        child_element_index[9].push_back(14-1);
    }
    else if (elem_type == PENTA_6)
    {
        //! PENTA_6

        //! 1  2  3  4  5  6
        //! Face 5

        //! 1  3  2
        //! 4  5  6
        //! 1  2  5  4
        //! 2  3  6  5
        //! 3  1  4  6

        face_element_type.resize(5);
        face_list.resize(5);

        face_element_type[0] = TRI_3;
        face_list[0].push_back(1-1);
        face_list[0].push_back(3-1);
        face_list[0].push_back(2-1);

        face_element_type[1] = TRI_3;
        face_list[1].push_back(4-1);
        face_list[1].push_back(5-1);
        face_list[1].push_back(6-1);

        face_element_type[2] = QUAD_4;
        face_list[2].push_back(1-1);
        face_list[2].push_back(2-1);
        face_list[2].push_back(5-1);
        face_list[2].push_back(4-1);

        face_element_type[3] = QUAD_4;
        face_list[3].push_back(2-1);
        face_list[3].push_back(3-1);
        face_list[3].push_back(6-1);
        face_list[3].push_back(5-1);

        face_element_type[4] = QUAD_4;
        face_list[4].push_back(3-1);
        face_list[4].push_back(1-1);
        face_list[4].push_back(4-1);
        face_list[4].push_back(6-1);
    }
    else if (elem_type == PHSPACE::PENTA_12)
    {
        //! PENTA_12

        //******************************
        //! Composite Face 5
        //! QUAD_8 : 1  2  5  4  7  10
        //! QUAD_8 : 2  3  6  5  8  11
        //! QUAD_8 : 3  1  4  6  9  12
        //! TRI_6  : 1  3  2  9  8  7
        //! TRI_6  : 4  5  6  10 11 12

        composite_face_element_type.resize(5);
        composite_face_list.resize(5);

        composite_face_element_type[ 0 ] = QUAD_6;
        composite_face_list[ 0 ].push_back(1 - 1);
        composite_face_list[ 0 ].push_back(2 - 1);
        composite_face_list[ 0 ].push_back(5 - 1);
        composite_face_list[ 0 ].push_back(4 - 1);
        composite_face_list[ 0 ].push_back(7 - 1);
        composite_face_list[ 0 ].push_back(10 - 1);

        composite_face_element_type[ 1 ] = QUAD_6;
        composite_face_list[ 1 ].push_back(2 - 1);
        composite_face_list[ 1 ].push_back(3 - 1);
        composite_face_list[ 1 ].push_back(6 - 1);
        composite_face_list[ 1 ].push_back(5 - 1);
        composite_face_list[ 1 ].push_back(8 - 1);
        composite_face_list[ 1 ].push_back(11 - 1);

        composite_face_element_type[ 2 ] = QUAD_6;
        composite_face_list[ 2 ].push_back(3 - 1);
        composite_face_list[ 2 ].push_back(1 - 1);
        composite_face_list[ 2 ].push_back(4 - 1);
        composite_face_list[ 2 ].push_back(6 - 1);
        composite_face_list[ 2 ].push_back(9 - 1);
        composite_face_list[ 2 ].push_back(12 - 1);

        composite_face_element_type[ 3 ] = TRI_6;
        composite_face_list[ 3 ].push_back(1 - 1);
        composite_face_list[ 3 ].push_back(3 - 1);
        composite_face_list[ 3 ].push_back(2 - 1);
        composite_face_list[ 3 ].push_back(9 - 1);
        composite_face_list[ 3 ].push_back(8 - 1);//
        composite_face_list[ 3 ].push_back(7 - 1);

        composite_face_element_type[ 4 ] = TRI_6;
        composite_face_list[ 4 ].push_back(4 - 1);
        composite_face_list[ 4 ].push_back(5 - 1);
        composite_face_list[ 4 ].push_back(6 - 1);
        composite_face_list[ 4 ].push_back(10 - 1);
        composite_face_list[ 4 ].push_back(11 - 1);
        composite_face_list[ 4 ].push_back(12 - 1);

        //! Face 20

        //! 1  9  3  8  2  7  1
        //! Face 4
        //! 1  9  7
        //! 9  8  7
        //! 7  8  2
        //! 9  3  8

        face_element_type.resize(14);
        face_list.resize(14);

        face_element_type[ 0 ] = TRI_3;
        face_list[ 0 ].push_back(1 - 1);
        face_list[ 0 ].push_back(9 - 1);
        face_list[ 0 ].push_back(7 - 1);

        face_element_type[ 1 ] = TRI_3;
        face_list[ 1 ].push_back(9 - 1);
        face_list[ 1 ].push_back(8 - 1);
        face_list[ 1 ].push_back(7 - 1);

        face_element_type[ 2 ] = TRI_3;
        face_list[ 2 ].push_back(7 - 1);
        face_list[ 2 ].push_back(8 - 1);
        face_list[ 2 ].push_back(2 - 1);

        face_element_type[ 3 ] = TRI_3;
        face_list[ 3 ].push_back(9 - 1);
        face_list[ 3 ].push_back(3 - 1);
        face_list[ 3 ].push_back(8 - 1);

        //! 4  13  5  14  6  15  4
        //! Face 4
        //! 4  13 15
        //! 13 14 15
        //! 13 5  14
        //! 15 14 6

        face_element_type[ 4 ] = TRI_3;
        face_list[ 4 ].push_back(4 - 1);
        face_list[ 4 ].push_back(10 - 1);
        face_list[ 4 ].push_back(12 - 1);

        face_element_type[ 5 ] = TRI_3;
        face_list[ 5 ].push_back(10 - 1);
        face_list[ 5 ].push_back(11 - 1);
        face_list[ 5 ].push_back(12 - 1);

        face_element_type[ 6 ] = TRI_3;
        face_list[ 6 ].push_back(10 - 1);//
        face_list[ 6 ].push_back(5 - 1);
        face_list[ 6 ].push_back(11 - 1);

        face_element_type[ 7 ] = TRI_3;
        face_list[ 7 ].push_back(12 - 1);
        face_list[ 7 ].push_back(11 - 1);
        face_list[ 7 ].push_back(6 - 1);

        //! 1  7  2  11  5  13  4  10  1  middle16
        //! Face 4
        //! 1  7  16 10
        //! 7  2  11 16
        //! 10 16 13 4
        //! 16 11 5  13
        face_element_type[ 8 ] = QUAD_4;
        face_list[ 8 ].push_back(1 - 1);
        face_list[ 8 ].push_back(7 - 1);
        face_list[ 8 ].push_back(10 - 1);
        face_list[ 8 ].push_back(4 - 1);

        face_element_type[ 9 ] = QUAD_4;
        face_list[ 9 ].push_back(7 - 1);
        face_list[ 9 ].push_back(2 - 1);
        face_list[ 9 ].push_back(5 - 1);
        face_list[ 9 ].push_back(10 - 1);

        //! 2  8  3  12  6  14  5  11  2  middle17
        //! Face 4
        //! 2  8  17 11
        //! 8  3  12 17
        //! 11 17 14 5
        //! 17 12 6  14

        face_element_type[ 10 ] = QUAD_4;
        face_list[ 10 ].push_back(2 - 1);
        face_list[ 10 ].push_back(8 - 1);
        face_list[ 10 ].push_back(11 - 1);
        face_list[ 10 ].push_back(5 - 1);

        face_element_type[ 11 ] = QUAD_4;
        face_list[ 11 ].push_back(8 - 1);
        face_list[ 11 ].push_back(3 - 1);
        face_list[ 11 ].push_back(6 - 1);
        face_list[ 11 ].push_back(11 - 1);

        //! 1  10  4  15  6  12  3  9  1  middle18
        //! Face 4
        //! 1  10 18 9
        //! 9  18 12 3
        //! 10 4  15 18
        //! 18 15 6  12

        face_element_type[ 12 ] = QUAD_4;
        face_list[ 12 ].push_back(3 - 1);
        face_list[ 12 ].push_back(9 - 1);
        face_list[ 12 ].push_back(12 - 1);
        face_list[ 12 ].push_back(6 - 1);

        face_element_type[ 13 ] = QUAD_4;
        face_list[ 13 ].push_back(9 - 1);
        face_list[ 13 ].push_back(1 - 1);
        face_list[ 13 ].push_back(4 - 1);
        face_list[ 13 ].push_back(12 - 1);
        //******************************

        //! PENTA_12

        mp_struct[ 6 ].push_back(0);
        mp_struct[ 6 ].push_back(1);

        mp_struct[ 7 ].push_back(1);
        mp_struct[ 7 ].push_back(2);

        mp_struct[ 8 ].push_back(2);
        mp_struct[ 8 ].push_back(0);//

        mp_struct[ 9 ].push_back(3);
        mp_struct[ 9 ].push_back(4);

        mp_struct[ 10 ].push_back(4);
        mp_struct[ 10 ].push_back(5);

        mp_struct[ 11 ].push_back(5);
        mp_struct[ 11 ].push_back(3);


        //! subcell

        child_element_type.resize(4);
        child_element_index.resize(4);

        child_element_type[ 0 ] = PENTA_6;
        child_element_index[ 0 ].push_back(1 - 1);
        child_element_index[ 0 ].push_back(7 - 1);
        child_element_index[ 0 ].push_back(9 - 1);
        child_element_index[ 0 ].push_back(4 - 1);
        child_element_index[ 0 ].push_back(10 - 1);
        child_element_index[ 0 ].push_back(12 - 1);

        child_element_type[ 1 ] = PENTA_6;
        child_element_index[ 1 ].push_back(7 - 1);
        child_element_index[ 1 ].push_back(8 - 1);
        child_element_index[ 1 ].push_back(9 - 1);
        child_element_index[ 1 ].push_back(10 - 1);
        child_element_index[ 1 ].push_back(11 - 1);
        child_element_index[ 1 ].push_back(12 - 1);

        child_element_type[ 2 ] = PENTA_6;
        child_element_index[ 2 ].push_back(7 - 1);
        child_element_index[ 2 ].push_back(2 - 1);
        child_element_index[ 2 ].push_back(8 - 1);
        child_element_index[ 2 ].push_back(10 - 1);
        child_element_index[ 2 ].push_back(5 - 1);
        child_element_index[ 2 ].push_back(11 - 1);


        child_element_type[ 3 ] = PENTA_6;
        child_element_index[ 3 ].push_back(9 - 1);
        child_element_index[ 3 ].push_back(8 - 1);
        child_element_index[ 3 ].push_back(3 - 1);
        child_element_index[ 3 ].push_back(12 - 1);
        child_element_index[ 3 ].push_back(11 - 1);
        child_element_index[ 3 ].push_back(6 - 1);
    }
    else if (elem_type == PENTA_15)
    {
        //! PENTA_15
        //! 1.........7..........2
        //! 2.........8..........3
        //! 3.........9..........1

        //! 1.........10.........4
        //! 2.........11.........5
        //! 3.........12.........6

        //! 4.........13..........5
        //! 5.........14..........6
        //! 6.........15..........4

        mp_struct[6].push_back(0);
        mp_struct[6].push_back(1);

        mp_struct[7].push_back(1);
        mp_struct[7].push_back(2);

        mp_struct[8].push_back(2);
        mp_struct[8].push_back(3);

        mp_struct[9].push_back(0);
        mp_struct[9].push_back(3);

        mp_struct[10].push_back(1);
        mp_struct[10].push_back(4);

        mp_struct[11].push_back(2);
        mp_struct[11].push_back(5);

        mp_struct[12].push_back(3);
        mp_struct[12].push_back(4);

        mp_struct[13].push_back(4);
        mp_struct[13].push_back(5);

        mp_struct[14].push_back(5);
        mp_struct[14].push_back(3);
    }
    else if (elem_type == PENTA_18)
    {
        //! PENTA_18

        //! Composite Face 5
        //! QUAD_9 : 1  2  5  4  7  11  13  10  16
        //! QUAD_9 : 2  3  6  5  8  12  14  11  17
        //! QUAD_9 : 3  1  4  6  9  10  15  12  18
        //! TRI_6  : 1  3  2  9  8  7
        //! TRI_6  : 4  5  6  13 14 15

        composite_face_element_type.resize(5);
        composite_face_list.resize(5);

        composite_face_element_type[0] = QUAD_9;
        composite_face_list[0].push_back(1-1);
        composite_face_list[0].push_back(2-1);
        composite_face_list[0].push_back(5-1);
        composite_face_list[0].push_back(4-1);
        composite_face_list[0].push_back(7-1);
        composite_face_list[0].push_back(11-1);
        composite_face_list[0].push_back(13-1);
        composite_face_list[0].push_back(10-1);
        composite_face_list[0].push_back(16-1);

        composite_face_element_type[1] = QUAD_9;
        composite_face_list[1].push_back(2-1);
        composite_face_list[1].push_back(3-1);
        composite_face_list[1].push_back(6-1);
        composite_face_list[1].push_back(5-1);
        composite_face_list[1].push_back(8-1);
        composite_face_list[1].push_back(12-1);
        composite_face_list[1].push_back(14-1);
        composite_face_list[1].push_back(11-1);
        composite_face_list[1].push_back(17-1);

        composite_face_element_type[2] = QUAD_9;
        composite_face_list[2].push_back(3-1);
        composite_face_list[2].push_back(1-1);
        composite_face_list[2].push_back(4-1);
        composite_face_list[2].push_back(6-1);
        composite_face_list[2].push_back(9-1);
        composite_face_list[2].push_back(10-1);
        composite_face_list[2].push_back(15-1);
        composite_face_list[2].push_back(12-1);
        composite_face_list[2].push_back(18-1);

        composite_face_element_type[3] = TRI_6;
        composite_face_list[3].push_back(1-1);
        composite_face_list[3].push_back(3-1);
        composite_face_list[3].push_back(2-1);
        composite_face_list[3].push_back(9-1);
        composite_face_list[3].push_back(8-1);
        composite_face_list[3].push_back(7-1);

        composite_face_element_type[4] = TRI_6;
        composite_face_list[4].push_back(4-1);
        composite_face_list[4].push_back(5-1);
        composite_face_list[4].push_back(6-1);
        composite_face_list[4].push_back(13-1);
        composite_face_list[4].push_back(14-1);
        composite_face_list[4].push_back(15-1);

        //! Face 20

        //! 1  9  3  8  2  7  1
        //! Face 4
        //! 1  9  7
        //! 9  8  7
        //! 7  8  2
        //! 9  3  8

        face_element_type.resize(20);
        face_list.resize(20);

        face_element_type[0] = TRI_3;
        face_list[0].push_back(1-1);
        face_list[0].push_back(9-1);
        face_list[0].push_back(7-1);

        face_element_type[1] = TRI_3;
        face_list[1].push_back(9-1);
        face_list[1].push_back(8-1);
        face_list[1].push_back(7-1);

        face_element_type[2] = TRI_3;
        face_list[2].push_back(7-1);
        face_list[2].push_back(8-1);
        face_list[2].push_back(2-1);

        face_element_type[3] = TRI_3;
        face_list[3].push_back(9-1);
        face_list[3].push_back(3-1);
        face_list[3].push_back(8-1);

        //! 4  13  5  14  6  15  4
        //! Face 4
        //! 4  13 15
        //! 13 14 15
        //! 12 5  14
        //! 15 14 6

        face_element_type[4] = TRI_3;
        face_list[4].push_back(4-1);
        face_list[4].push_back(13-1);
        face_list[4].push_back(15-1);

        face_element_type[5] = TRI_3;
        face_list[5].push_back(13-1);
        face_list[5].push_back(14-1);
        face_list[5].push_back(15-1);

        face_element_type[6] = TRI_3;
        face_list[6].push_back(13-1);
        face_list[6].push_back(5-1);
        face_list[6].push_back(14-1);

        face_element_type[7] = TRI_3;
        face_list[7].push_back(15-1);
        face_list[7].push_back(14-1);
        face_list[7].push_back(6-1);

        //! 1  7  2  11  5  13  4  10  1  middle16
        //! Face 4
        //! 1  7  16 10
        //! 7  2  11 16
        //! 10 16 13 4
        //! 16 11 5  13

        face_element_type[8] = QUAD_4;
        face_list[8].push_back(1-1);
        face_list[8].push_back(7-1);
        face_list[8].push_back(16-1);
        face_list[8].push_back(10-1);

        face_element_type[9] = QUAD_4;
        face_list[9].push_back(7-1);
        face_list[9].push_back(2-1);
        face_list[9].push_back(11-1);
        face_list[9].push_back(16-1);

        face_element_type[10] = QUAD_4;
        face_list[10].push_back(10-1);
        face_list[10].push_back(16-1);
        face_list[10].push_back(13-1);
        face_list[10].push_back(4-1);

        face_element_type[11] = QUAD_4;
        face_list[11].push_back(16-1);
        face_list[11].push_back(11-1);
        face_list[11].push_back(5-1);
        face_list[11].push_back(13-1);

        //! 2  8  3  12  6  14  5  11  2  middle17
        //! Face 4
        //! 2  8  17 11
        //! 8  3  12 17
        //! 11 17 14 5
        //! 17 12 6  14

        face_element_type[12] = QUAD_4;
        face_list[12].push_back(2-1);
        face_list[12].push_back(8-1);
        face_list[12].push_back(17-1);
        face_list[12].push_back(11-1);

        face_element_type[13] = QUAD_4;
        face_list[13].push_back(8-1);
        face_list[13].push_back(3-1);
        face_list[13].push_back(12-1);
        face_list[13].push_back(17-1);

        face_element_type[14] = QUAD_4;
        face_list[14].push_back(11-1);
        face_list[14].push_back(17-1);
        face_list[14].push_back(14-1);
        face_list[14].push_back(5-1);

        face_element_type[15] = QUAD_4;
        face_list[15].push_back(17-1);
        face_list[15].push_back(12-1);
        face_list[15].push_back(6-1);
        face_list[15].push_back(14-1);

        //! 1  10  4  15  6  12  3  9  1  middle18
        //! Face 4
        //! 1  10 18 9
        //! 9  18 12 3
        //! 10 4  15 18
        //! 18 15 6  12

        face_element_type[16] = QUAD_4;
        face_list[16].push_back(1-1);
        face_list[16].push_back(10-1);
        face_list[16].push_back(18-1);
        face_list[16].push_back(9-1);

        face_element_type[17] = QUAD_4;
        face_list[17].push_back(9-1);
        face_list[17].push_back(18-1);
        face_list[17].push_back(12-1);
        face_list[17].push_back(3-1);

        face_element_type[18] = QUAD_4;
        face_list[18].push_back(10-1);
        face_list[18].push_back(4-1);
        face_list[18].push_back(15-1);
        face_list[18].push_back(18-1);

        face_element_type[19] = QUAD_4;
        face_list[19].push_back(18-1);
        face_list[19].push_back(15-1);
        face_list[19].push_back(6-1);
        face_list[19].push_back(12-1);

        //! PENTA_18
        //! 1.........7..........2
        //! 2.........8..........3
        //! 3.........9..........1

        //! 1.........10.........4
        //! 2.........11.........5
        //! 3.........12.........6

        //! 4.........13.........5
        //! 5.........14.........6
        //! 6.........15.........4

        //! 10........16.........11
        //! 11........17.........12
        //! 12........18.........10

        mp_struct[6].push_back(0);
        mp_struct[6].push_back(1);

        mp_struct[7].push_back(1);
        mp_struct[7].push_back(2);

        mp_struct[8].push_back(2);
        mp_struct[8].push_back(0);

        mp_struct[9].push_back(0);
        mp_struct[9].push_back(3);

        mp_struct[10].push_back(1);
        mp_struct[10].push_back(4);

        mp_struct[11].push_back(2);
        mp_struct[11].push_back(5);

        mp_struct[12].push_back(3);
        mp_struct[12].push_back(4);

        mp_struct[13].push_back(4);
        mp_struct[13].push_back(5);

        mp_struct[14].push_back(5);
        mp_struct[14].push_back(3);

        mp_struct[ 15 ].push_back(0);//
        mp_struct[ 15 ].push_back(1);//
        mp_struct[ 15 ].push_back(4);//
        mp_struct[ 15 ].push_back(3);//

        mp_struct[ 16 ].push_back(1);//
        mp_struct[ 16 ].push_back(2);//
        mp_struct[ 16 ].push_back(5);//
        mp_struct[ 16 ].push_back(4);//

        mp_struct[ 17 ].push_back(2);//
        mp_struct[ 17 ].push_back(0);//
        mp_struct[ 17 ].push_back(3);//
        mp_struct[ 17 ].push_back(5);//

        //! subcell
        //! PENTA_6  : 1  7  9  10 16 18
        //! PENTA_6  : 7  8  9  16 17 18
        //! PENTA_6  : 7  2  8  16 11 17
        //! PENTA_6  : 9  8  3  18 17 12
        //! PENTA_6  : 10 16 18 4  13 15
        //! PENTA_6  : 16 17 18 13 14 15
        //! PENTA_6  : 16 11 17 12 5  14
        //! PENTA_6  : 18 17 12 15 14 6

        child_element_type.resize(8);
        child_element_index.resize(8);

        child_element_type[0] = PENTA_6;
        child_element_index[0].push_back(1-1);
        child_element_index[0].push_back(7-1);
        child_element_index[0].push_back(9-1);
        child_element_index[0].push_back(10-1);
        child_element_index[0].push_back(16-1);
        child_element_index[0].push_back(18-1);

        child_element_type[1] = PENTA_6;
        child_element_index[1].push_back(7-1);
        child_element_index[1].push_back(8-1);
        child_element_index[1].push_back(9-1);
        child_element_index[1].push_back(16-1);
        child_element_index[1].push_back(17-1);
        child_element_index[1].push_back(18-1);

        child_element_type[2] = PENTA_6;
        child_element_index[2].push_back(7-1);
        child_element_index[2].push_back(2-1);
        child_element_index[2].push_back(8-1);
        child_element_index[2].push_back(16-1);
        child_element_index[2].push_back(11-1);
        child_element_index[2].push_back(17-1);


        child_element_type[3] = PENTA_6;
        child_element_index[3].push_back(9-1);
        child_element_index[3].push_back(8-1);
        child_element_index[3].push_back(3-1);
        child_element_index[3].push_back(18-1);
        child_element_index[3].push_back(17-1);
        child_element_index[3].push_back(12-1);

        child_element_type[4] = PENTA_6;
        child_element_index[4].push_back(10-1);
        child_element_index[4].push_back(16-1);
        child_element_index[4].push_back(18-1);
        child_element_index[4].push_back(4-1);
        child_element_index[4].push_back(13-1);
        child_element_index[4].push_back(15-1);

        child_element_type[5] = PENTA_6;
        child_element_index[5].push_back(16-1);
        child_element_index[5].push_back(17-1);
        child_element_index[5].push_back(18-1);
        child_element_index[5].push_back(13-1);
        child_element_index[5].push_back(14-1);
        child_element_index[5].push_back(15-1);

        child_element_type[6] = PENTA_6;
        child_element_index[6].push_back(16-1);
        child_element_index[6].push_back(11-1);
        child_element_index[6].push_back(17-1);
        child_element_index[6].push_back(13-1);
        child_element_index[6].push_back(5-1);
        child_element_index[6].push_back(14-1);

        child_element_type[7] = PENTA_6;
        child_element_index[7].push_back(18-1);
        child_element_index[7].push_back(17-1);
        child_element_index[7].push_back(12-1);
        child_element_index[7].push_back(15-1);
        child_element_index[7].push_back(14-1);
        child_element_index[7].push_back(6-1);
    }
    else if (elem_type == HEXA_8)
    {
        //! HEXA_8

        //! 1  2 3 4 5 6 7 8
        //! Face 6
        //! 1  4  3  2
        //! 5  6  7  8
        //! 1  5  8  4
        //! 2  3  7  6
        //! 4  3  8  7
        //! 1  2  6  5

        face_element_type.resize(6);
        face_list.resize(6);

        face_element_type[0] = QUAD_4;
        face_list[0].push_back(1-1);
        face_list[0].push_back(4-1);
        face_list[0].push_back(3-1);
        face_list[0].push_back(2-1);

        face_element_type[1] = QUAD_4;
        face_list[1].push_back(5-1);
        face_list[1].push_back(6-1);
        face_list[1].push_back(7-1);
        face_list[1].push_back(8-1);

        face_element_type[2] = QUAD_4;
        face_list[2].push_back(1-1);
        face_list[2].push_back(5-1);
        face_list[2].push_back(8-1);
        face_list[2].push_back(4-1);

        face_element_type[3] = QUAD_4;
        face_list[3].push_back(2-1);
        face_list[3].push_back(3-1);
        face_list[3].push_back(7-1);
        face_list[3].push_back(6-1);

        face_element_type[4] = QUAD_4;
        face_list[4].push_back(3-1);
        face_list[4].push_back(4-1);
        face_list[4].push_back(8-1);
        face_list[4].push_back(7-1);

        face_element_type[5] = QUAD_4;
        face_list[5].push_back(1-1);
        face_list[5].push_back(2-1);
        face_list[5].push_back(6-1);
        face_list[5].push_back(5-1);
    }
    else if (elem_type == HEXA_20)
    {
        //! HEXA_20
        //! 1.........9. ........2
        //! 2.........10.........3
        //! 3.........11.........4
        //! 4.........12.........1

        //! 1.........13.........5
        //! 2.........14.........6
        //! 3.........15.........7
        //! 4.........16.........8

        //! 5.........17.........6
        //! 6.........18.........7
        //! 7.........19.........8
        //! 8.........20.........5

        mp_struct[8].push_back(0);
        mp_struct[8].push_back(1);

        mp_struct[9].push_back(1);
        mp_struct[9].push_back(2);

        mp_struct[10].push_back(2);
        mp_struct[10].push_back(3);

        mp_struct[11].push_back(3);
        mp_struct[11].push_back(0);

        mp_struct[12].push_back(0);
        mp_struct[12].push_back(4);

        mp_struct[13].push_back(1);
        mp_struct[13].push_back(5);

        mp_struct[14].push_back(2);
        mp_struct[14].push_back(6);

        mp_struct[15].push_back(3);
        mp_struct[15].push_back(7);

        mp_struct[16].push_back(4);
        mp_struct[16].push_back(5);

        mp_struct[17].push_back(5);
        mp_struct[17].push_back(6);

        mp_struct[18].push_back(6);
        mp_struct[18].push_back(7);

        mp_struct[19].push_back(7);
        mp_struct[19].push_back(4);
    }
    else if (elem_type == HEXA_27)
    {
        //! HEXA_27

        //! Composite Face 6
        //! QUAD_9 : 1  4  3  2  12  11  10  9   21
        //! QUAD_9 : 5  6  7  8  17  18  19  20  26
        //! QUAD_9 : 1  2  6  5  9   14  17  13  22
        //! QUAD_9 : 2  3  7  6  10  15  18  14  23
        //! QUAD_9 : 3  4  8  7  11  16  19  15  24
        //! QUAD_9 : 1  4  8  5  12  16  20  10  25

        composite_face_element_type.resize(6);
        composite_face_list.resize(6);

        composite_face_element_type[0] = QUAD_9;
        composite_face_list[0].push_back(1-1);
        composite_face_list[0].push_back(4-1);
        composite_face_list[0].push_back(3-1);
        composite_face_list[0].push_back(2-1);
        composite_face_list[0].push_back(12-1);
        composite_face_list[0].push_back(11-1);
        composite_face_list[0].push_back(10-1);
        composite_face_list[0].push_back(9-1);
        composite_face_list[0].push_back(21-1);

        composite_face_element_type[1] = QUAD_9;
        composite_face_list[1].push_back(5-1);
        composite_face_list[1].push_back(6-1);
        composite_face_list[1].push_back(7-1);
        composite_face_list[1].push_back(8-1);
        composite_face_list[1].push_back(17-1);
        composite_face_list[1].push_back(18-1);
        composite_face_list[1].push_back(19-1);
        composite_face_list[1].push_back(20-1);
        composite_face_list[1].push_back(26-1);

        composite_face_element_type[2] = QUAD_9;
        composite_face_list[2].push_back(1-1);
        composite_face_list[2].push_back(2-1);
        composite_face_list[2].push_back(6-1);
        composite_face_list[2].push_back(5-1);
        composite_face_list[2].push_back(9-1);
        composite_face_list[2].push_back(14-1);
        composite_face_list[2].push_back(17-1);
        composite_face_list[2].push_back(13-1);
        composite_face_list[2].push_back(22-1);

        composite_face_element_type[3] = QUAD_9;
        composite_face_list[3].push_back(2-1);
        composite_face_list[3].push_back(3-1);
        composite_face_list[3].push_back(7-1);
        composite_face_list[3].push_back(6-1);
        composite_face_list[3].push_back(10-1);
        composite_face_list[3].push_back(15-1);
        composite_face_list[3].push_back(18-1);
        composite_face_list[3].push_back(14-1);
        composite_face_list[3].push_back(23-1);

        composite_face_element_type[4] = QUAD_9;
        composite_face_list[4].push_back(3-1);
        composite_face_list[4].push_back(4-1);
        composite_face_list[4].push_back(8-1);
        composite_face_list[4].push_back(7-1);
        composite_face_list[4].push_back(11-1);
        composite_face_list[4].push_back(16-1);
        composite_face_list[4].push_back(19-1);
        composite_face_list[4].push_back(15-1);
        composite_face_list[4].push_back(24-1);

        composite_face_element_type[5] = QUAD_9;
        composite_face_list[5].push_back(1-1);
        composite_face_list[5].push_back(5-1);
        composite_face_list[5].push_back(8-1);
        composite_face_list[5].push_back(4-1);
        composite_face_list[5].push_back(13-1);
        composite_face_list[5].push_back(20-1);
        composite_face_list[5].push_back(16-1);
        composite_face_list[5].push_back(12-1);
        composite_face_list[5].push_back(25-1);

        //! Face 24

        //! 1 12 4 11 3  10  2  9  1  middle 21
        //! Face 4
        //! 1  12 21 9
        //! 12 4  11 21
        //! 9  21 10 2
        //! 21 11 3 10

        face_element_type.resize(24);
        face_list.resize(24);

        face_element_type[0] = QUAD_4;
        face_list[0].push_back(1-1);
        face_list[0].push_back(12-1);
        face_list[0].push_back(21-1);
        face_list[0].push_back(9-1);

        face_element_type[1] = QUAD_4;
        face_list[1].push_back(12-1);
        face_list[1].push_back(4-1);
        face_list[1].push_back(11-1);
        face_list[1].push_back(21-1);

        face_element_type[2] = QUAD_4;
        face_list[2].push_back(9-1);
        face_list[2].push_back(21-1);
        face_list[2].push_back(10-1);
        face_list[2].push_back(2-1);

        face_element_type[3] = QUAD_4;
        face_list[3].push_back(21-1);
        face_list[3].push_back(11-1);
        face_list[3].push_back(3-1);
        face_list[3].push_back(10-1);

        //! 5 17 6  18 7 19 8 20 5  middle 26
        //! Face 4
        //! 5  17 26 20
        //! 20 26 19 8
        //! 17 6  18 26
        //! 26 18 7  19

        face_element_type[4] = QUAD_4;
        face_list[4].push_back(5-1);
        face_list[4].push_back(17-1);
        face_list[4].push_back(26-1);
        face_list[4].push_back(20-1);

        face_element_type[5] = QUAD_4;
        face_list[5].push_back(20-1);
        face_list[5].push_back(26-1);
        face_list[5].push_back(19-1);
        face_list[5].push_back(8-1);

        face_element_type[6] = QUAD_4;
        face_list[6].push_back(17-1);
        face_list[6].push_back(6-1);
        face_list[6].push_back(18-1);
        face_list[6].push_back(26-1);

        face_element_type[7] = QUAD_4;
        face_list[7].push_back(26-1);
        face_list[7].push_back(18-1);
        face_list[7].push_back(7-1);
        face_list[7].push_back(19-1);

        //! 1  13  5  20  8  16 4 12 1  middle 25
        //! Face 4
        //! 1  13 25 12
        //! 12 25 16  4
        //! 13 5  20 25
        //! 25 20 8  16

        face_element_type[8] = QUAD_4;
        face_list[8].push_back(1-1);
        face_list[8].push_back(13-1);
        face_list[8].push_back(25-1);
        face_list[8].push_back(12-1);

        face_element_type[9] = QUAD_4;
        face_list[9].push_back(12-1);
        face_list[9].push_back(25-1);
        face_list[9].push_back(16-1);
        face_list[9].push_back(4-1);

        face_element_type[10] = QUAD_4;
        face_list[10].push_back(13-1);
        face_list[10].push_back(5-1);
        face_list[10].push_back(20-1);
        face_list[10].push_back(25-1);

        face_element_type[11] = QUAD_4;
        face_list[11].push_back(25-1);
        face_list[11].push_back(20-1);
        face_list[11].push_back(8-1);
        face_list[11].push_back(16-1);

        //! 2  10  3  15  7 18 6 14 2  middle 23
        //! Face 4
        //! 2  10 23 14
        //! 10 3  15 23
        //! 14 23 18 6
        //! 23 15 7 18

        face_element_type[12] = QUAD_4;
        face_list[12].push_back(2-1);
        face_list[12].push_back(10-1);
        face_list[12].push_back(23-1);
        face_list[12].push_back(14-1);

        face_element_type[13] = QUAD_4;
        face_list[13].push_back(10-1);
        face_list[13].push_back(3-1);
        face_list[13].push_back(15-1);
        face_list[13].push_back(23-1);

        face_element_type[14] = QUAD_4;
        face_list[14].push_back(14-1);
        face_list[14].push_back(23-1);
        face_list[14].push_back(18-1);
        face_list[14].push_back(6-1);

        face_element_type[15] = QUAD_4;
        face_list[15].push_back(23-1);
        face_list[15].push_back(15-1);
        face_list[15].push_back(7-1);
        face_list[15].push_back(18-1);

        //! 4 16 8 19  7 15 3 11 4  middle 24
        //! Face 4
        //! 4 16 24 11
        //! 11 24 15 3
        //! 16 8 19 24
        //! 24 19 7 15

        face_element_type[16] = QUAD_4;
        face_list[16].push_back(4-1);
        face_list[16].push_back(16-1);
        face_list[16].push_back(24-1);
        face_list[16].push_back(11-1);

        face_element_type[17] = QUAD_4;
        face_list[17].push_back(11-1);
        face_list[17].push_back(24-1);
        face_list[17].push_back(15-1);
        face_list[17].push_back(3-1);

        face_element_type[18] = QUAD_4;
        face_list[18].push_back(16-1);
        face_list[18].push_back(8-1);
        face_list[18].push_back(19-1);
        face_list[18].push_back(24-1);

        face_element_type[19] = QUAD_4;
        face_list[19].push_back(24-1);
        face_list[19].push_back(19-1);
        face_list[19].push_back(7-1);
        face_list[19].push_back(15-1);

        //! 1 9 2 14 6 17 5 13 1  middle 22
        //! Face 4
        //! 1  9  22 10
        //! 9  2  14 22
        //! 10 22 17 5
        //! 22 14 6 17

        face_element_type[20] = QUAD_4;
        face_list[20].push_back(1-1);
        face_list[20].push_back(9-1);
        face_list[20].push_back(22-1);
        face_list[20].push_back(13-1);

        face_element_type[21] = QUAD_4;
        face_list[21].push_back(9-1);
        face_list[21].push_back(2-1);
        face_list[21].push_back(14-1);
        face_list[21].push_back(22-1);

        face_element_type[22] = QUAD_4;
        face_list[22].push_back(13-1);
        face_list[22].push_back(22-1);
        face_list[22].push_back(17-1);
        face_list[22].push_back(5-1);

        face_element_type[23] = QUAD_4;
        face_list[23].push_back(22-1);
        face_list[23].push_back(14-1);
        face_list[23].push_back(6-1);
        face_list[23].push_back(17-1);

        //! HEXA_27
        //! 1.........9. ........2
        //! 2.........10.........3
        //! 3.........11.........4
        //! 4.........12.........1

        //! 1.........13.........5
        //! 2.........14.........6
        //! 3.........15.........7
        //! 4.........16.........8

        //! 5.........17.........6
        //! 6.........18.........7
        //! 7.........19.........8
        //! 8.........20.........5

        //! 1 2.......21.........3 4
        //! 1 2.......22.........6 5
        //! 2 6.......23.........7 3
        //! 3 7.......24.........8 4
        //! 1 5.......25.........8 4
        //! 5 6.......26.........7 8
        //! 1 2 3 4...27.........5 6 7 8

        mp_struct[8].push_back(0);
        mp_struct[8].push_back(1);

        mp_struct[9].push_back(1);
        mp_struct[9].push_back(2);

        mp_struct[10].push_back(2);
        mp_struct[10].push_back(3);

        mp_struct[11].push_back(3);
        mp_struct[11].push_back(0);

        mp_struct[12].push_back(0);
        mp_struct[12].push_back(4);

        mp_struct[13].push_back(1);
        mp_struct[13].push_back(5);

        mp_struct[14].push_back(2);
        mp_struct[14].push_back(6);

        mp_struct[15].push_back(3);
        mp_struct[15].push_back(7);

        mp_struct[16].push_back(4);
        mp_struct[16].push_back(5);

        mp_struct[17].push_back(5);
        mp_struct[17].push_back(6);

        mp_struct[18].push_back(6);
        mp_struct[18].push_back(7);

        mp_struct[19].push_back(7);
        mp_struct[19].push_back(4);

        //! 1 2.......21.........3 4
        //! 1 2.......22.........6 5
        //! 2 6.......23.........7 3
        //! 3 7.......24.........8 4
        //! 1 5.......25.........8 4
        //! 5 6.......26.........7 8
        //! 1 2 3 4...27.........5 6 7 8

        mp_struct[20].push_back(0);
        mp_struct[20].push_back(3);
        mp_struct[20].push_back(2);
        mp_struct[20].push_back(1);

        mp_struct[21].push_back(0);
        mp_struct[21].push_back(1);
        mp_struct[21].push_back(5);
        mp_struct[21].push_back(4);

        mp_struct[22].push_back(1);
        mp_struct[22].push_back(2);
        mp_struct[22].push_back(6);
        mp_struct[22].push_back(5);

        mp_struct[23].push_back(2);
        mp_struct[23].push_back(3);
        mp_struct[23].push_back(7);
        mp_struct[23].push_back(6);

        mp_struct[24].push_back(0);
        mp_struct[24].push_back(4);
        mp_struct[24].push_back(7);
        mp_struct[24].push_back(3);

        mp_struct[25].push_back(4);
        mp_struct[25].push_back(5);
        mp_struct[25].push_back(6);
        mp_struct[25].push_back(7);

        mp_struct[26].push_back(0);
        mp_struct[26].push_back(1);
        mp_struct[26].push_back(2);
        mp_struct[26].push_back(3);
        mp_struct[26].push_back(4);
        mp_struct[26].push_back(5);
        mp_struct[26].push_back(6);
        mp_struct[26].push_back(7);

        //! subcell
        //! HEXA_8  : 4  12 21 11 16 25 27 24
        //! HEXA_8  : 12 1  9  21 25 13 22 27
        //! HEXA_8  : 11 21 10 3  24 27 23 15
        //! HEXA_8  : 21 9  2  10 27 22 14 23
        //! HEXA_8  : 16 25 27 24 8  20 26 19
        //! HEXA_8  : 25 13 22 27 20 5  17 26
        //! HEXA_8  : 24 27 23 15 19 26 18 7
        //! HEXA_8  : 27 22 14 23 26 17 6  18

        child_element_type.resize(8);
        child_element_index.resize(8);

        child_element_type[0] = HEXA_8;
        child_element_index[0].push_back(4-1);
        child_element_index[0].push_back(12-1);
        child_element_index[0].push_back(21-1);
        child_element_index[0].push_back(11-1);
        child_element_index[0].push_back(16-1);
        child_element_index[0].push_back(25-1);
        child_element_index[0].push_back(27-1);
        child_element_index[0].push_back(24-1);

        child_element_type[1] = HEXA_8;
        child_element_index[1].push_back(12-1);
        child_element_index[1].push_back(1-1);
        child_element_index[1].push_back(9-1);
        child_element_index[1].push_back(21-1);
        child_element_index[1].push_back(25-1);
        child_element_index[1].push_back(13-1);
        child_element_index[1].push_back(22-1);
        child_element_index[1].push_back(27-1);

        child_element_type[2] = HEXA_8;
        child_element_index[2].push_back(11-1);
        child_element_index[2].push_back(21-1);
        child_element_index[2].push_back(10-1);
        child_element_index[2].push_back(3-1);
        child_element_index[2].push_back(24-1);
        child_element_index[2].push_back(27-1);
        child_element_index[2].push_back(23-1);
        child_element_index[2].push_back(15-1);

        child_element_type[3] = HEXA_8;
        child_element_index[3].push_back(21-1);
        child_element_index[3].push_back(9-1);
        child_element_index[3].push_back(2-1);
        child_element_index[3].push_back(10-1);
        child_element_index[3].push_back(27-1);
        child_element_index[3].push_back(22-1);
        child_element_index[3].push_back(14-1);
        child_element_index[3].push_back(23-1);

        child_element_type[4] = HEXA_8;
        child_element_index[4].push_back(16-1);
        child_element_index[4].push_back(25-1);
        child_element_index[4].push_back(27-1);
        child_element_index[4].push_back(24-1);
        child_element_index[4].push_back(8-1);
        child_element_index[4].push_back(20-1);
        child_element_index[4].push_back(26-1);
        child_element_index[4].push_back(19-1);

        child_element_type[5] = HEXA_8;
        child_element_index[5].push_back(25-1);
        child_element_index[5].push_back(13-1);
        child_element_index[5].push_back(22-1);
        child_element_index[5].push_back(27-1);
        child_element_index[5].push_back(20-1);
        child_element_index[5].push_back(5-1);
        child_element_index[5].push_back(17-1);
        child_element_index[5].push_back(26-1);

        child_element_type[6] = HEXA_8;
        child_element_index[6].push_back(24-1);
        child_element_index[6].push_back(27-1);
        child_element_index[6].push_back(23-1);
        child_element_index[6].push_back(15-1);
        child_element_index[6].push_back(19-1);
        child_element_index[6].push_back(26-1);
        child_element_index[6].push_back(18-1);
        child_element_index[6].push_back(7-1);

        child_element_type[7] = HEXA_8;
        child_element_index[7].push_back(27-1);
        child_element_index[7].push_back(22-1);
        child_element_index[7].push_back(14-1);
        child_element_index[7].push_back(23-1);
        child_element_index[7].push_back(26-1);
        child_element_index[7].push_back(17-1);
        child_element_index[7].push_back(6-1);
        child_element_index[7].push_back(18-1);
    }
}

bool BasicElement::ISBasicElement(int elem_type)
{
    bool results = false;
    if (elem_type == NODE)
    {
    }
    else if (elem_type == BAR_2)
    {
        results = true;
    }
    else if (elem_type == BAR_3)
    {
    }
    else if (elem_type == TRI_3)
    {
        if (GetDim() == TWO_D)
        {
            results = true;
        }
    }
    else if (elem_type == TRI_6)
    {
    }
    else if (elem_type == QUAD_4)
    {
        if (GetDim() == TWO_D)
        {
            results = true;
        }
    }
    else if (elem_type == QUAD_8)
    {
    }
    else if (elem_type == QUAD_9)
    {
    }
    else if (elem_type == TETRA_4)
    {
        results = true;
    }
    else if (elem_type == TETRA_10)
    {
    }
    else if (elem_type == PYRA_5)
    {
        results = true;
    }
    else if (elem_type == PYRA_14)
    {
    }
    else if (elem_type == PENTA_6)
    {
        results = true;
    }
    else if (elem_type == PENTA_15)
    {
    }
    else if (elem_type == PENTA_18)
    {
    }
    else if (elem_type == HEXA_8)
    {
        results = true;
    }
    else if (elem_type == HEXA_20)
    {
    }
    else if (elem_type == HEXA_27)
    {
    }
    return results;
}

bool BasicElement::IsFaceElementType(int elementType)
{
    bool results = false;

    if (elementType == BAR_2)
    {
        if (GetDim() == TWO_D)
        {
            results = true;
        }
        else
        {
            results = false;
        }
    }
    else if (elementType == TRI_3)
    {
        if (GetDim() == TWO_D)
        {
            results = false;
        }
        else
        {
            results = true;
        }
    }
    else if (elementType == QUAD_4)
    {
        if (GetDim() == TWO_D)
        {
            results = false;
        }
        else
        {
            results = true;
        }
    }
    return results;
}


bool BasicElement::IsFaceElementTypeAtLeast(int elementType)
{
    bool results = ISBasicElement(elementType);

    if (elementType == BAR_2)
    {
        if (GetDim() == TWO_D)
        {
            results = true;
        }
        else
        {
            results = false;
        }
    }
    return results;
}


bool BasicElement::IsBasicVolumeElementType(int elementType)
{
    bool results = false;
    if (elementType == TRI_3)
    {
        if (GetDim() == TWO_D)
        {
            results = true;
        }
    }
    else if (elementType == QUAD_4)
    {
        if (GetDim() == TWO_D)
        {
            results = true;
        }
    }
    else if (elementType == TETRA_4)
    {
        results = true;
    }
    else if (elementType == PYRA_5)
    {
        results = true;
    }
    else if (elementType == PENTA_6)
    {
        results = true;
    }
    else if (elementType == HEXA_8)
    {
        results = true;
    }
    return results;
}

int BasicElement::GetRefineElemType(int elem_type, const bool & anisotropicAdapt)
{
    int next_elem_type = -1;
    if (elem_type == NODE)
    {
    }
    else if (elem_type == BAR_2)
    {
        next_elem_type = BAR_3;
    }
    else if (elem_type == BAR_3)
    {
    }
    else if (elem_type == TRI_3)
    {
        next_elem_type = TRI_6;
    }
    else if (elem_type == TRI_6)
    {
    }
    else if (elem_type == QUAD_4)
    {
        if (anisotropicAdapt)
        {
            next_elem_type = QUAD_6;
        }
        else
        {
            next_elem_type = QUAD_9;
        }
    }
    else if (elem_type == QUAD_8)
    {
    }
    else if (elem_type == QUAD_9)
    {
    }
    else if (elem_type == TETRA_4)
    {
        next_elem_type = TETRA_10;
    }
    else if (elem_type == TETRA_10)
    {
    }
    else if (elem_type == PYRA_5)
    {
        next_elem_type = PYRA_14;
    }
    else if (elem_type == PYRA_14)
    {
    }
    else if (elem_type == PENTA_6)
    {
        if (anisotropicAdapt)
        {
            next_elem_type = PENTA_12;
        }
        else
        {
            next_elem_type = PENTA_18;
        }
    }
    else if (elem_type == PENTA_15)
    {
    }
    else if (elem_type == PENTA_18)
    {
    }
    else if (elem_type == HEXA_8)
    {
        next_elem_type = HEXA_27;
    }
    else if (elem_type == HEXA_20)
    {
    }
    else if (elem_type == HEXA_27)
    {
    }
    return next_elem_type;
}

int BasicElement::GetSimpleElemType(int elem_type)
{
    int prev_elem_type = -1;
    if (elem_type == NODE)
    {
    }
    else if (elem_type == BAR_2)
    {
    }
    else if (elem_type == BAR_3)
    {
        prev_elem_type = BAR_2;
    }
    else if (elem_type == TRI_3)
    {
    }
    else if (elem_type == TRI_6)
    {
        prev_elem_type = TRI_3;
    }
    else if (elem_type == QUAD_4)
    {
    }
    else if (elem_type == QUAD_6)
    {
        prev_elem_type = QUAD_4;
    }
    else if (elem_type == QUAD_8)
    {
    }
    else if (elem_type == QUAD_9)
    {
        prev_elem_type = QUAD_4;
    }
    else if (elem_type == TETRA_4)
    {
    }
    else if (elem_type == TETRA_10)
    {
        prev_elem_type = TETRA_4;
    }
    else if (elem_type == PYRA_5)
    {
    }
    else if (elem_type == PYRA_14)
    {
        prev_elem_type = PYRA_5;
    }
    else if (elem_type == PENTA_6)
    {
    }
    else if (elem_type == PENTA_12)
    {
        prev_elem_type = PENTA_6;
    }
    else if (elem_type == PENTA_15)
    {
    }
    else if (elem_type == PENTA_18)
    {
        prev_elem_type = PENTA_6;
    }
    else if (elem_type == HEXA_8)
    {
    }
    else if (elem_type == HEXA_20)
    {
    }
    else if (elem_type == HEXA_27)
    {
        prev_elem_type = HEXA_8;
    }
    return prev_elem_type;
}


ElementProxy::ElementProxy()
{
    basic_element_number = GetNumberOfBasicElement();
    basic_element = new BasicElement * [ basic_element_number ];
    for (int i = 0; i < basic_element_number; ++ i)
    {
        basic_element[i] = new BasicElement();
        basic_element[i]->Init(i);
    }
}

ElementProxy::~ElementProxy()
{
    for (int i = 0; i < basic_element_number; ++ i)
    {
        delete basic_element[i];
    }
    delete [] basic_element;
}

//! The following sections is copy from HyperFLOW 3.0.
//! Copy start.
ElementProxy * elementProxy = 0;
ElementInitializingClass::ElementInitializingClass()
{
    elementProxy = new ElementProxy();
}

ElementInitializingClass::~ElementInitializingClass()
{
    delete elementProxy;
}

ElementInitializingClass element_initializing;

int  GetRefineElemType(int elementType, const bool & anisotropicAdapt)
{ 
    return GetBasicElement(elementType)->GetRefineElemType(elementType, anisotropicAdapt);
}

int  GetSimpleElementType(int elementType)
{ 
    return GetBasicElement(elementType)->GetSimpleElemType(elementType);
}

bool ISBasicElement(int elementType)
{ 
    return GetBasicElement(elementType)->ISBasicElement(elementType); 
}

bool IsFaceElementType(int elementType)
{
    return GetBasicElement(elementType)->IsFaceElementType(elementType); 
}

bool IsFaceElementTypeAtLeast(int elementType)
{
    return GetBasicElement(elementType)->IsFaceElementTypeAtLeast(elementType); 
}

BasicElement * GetBasicElement(int elementType)
{ 
    return elementProxy->GetBasicElement(elementType);
}

uint_t GetChildElementNumbers(int elementType)
{ 
    return elementProxy->GetBasicElement(elementType)->GetChildElementNumber();
}

bool IsBasicVolumeElementType(int elementType)
{
    BasicElement * element = GetBasicElement(elementType);
    return element->ISBasicElement(elementType);
}

uint_t GetChildElementNumber(int elementType)
{ 
    return elementProxy->GetBasicElement(elementType)->GetChildElementNumber();
}

int GetElementNodeNumbers(int elementType) 
{ 
    return elementProxy->GetBasicElement(elementType)->GetElementNodeNumber(elementType);
}

vector<int> & GetRelatedPointListForMiddlePointComputation(int elem_type, int i)
{ 
    return elementProxy->GetBasicElement(elem_type)->GetMPlist(i);
}

bool IsBoundaryFace(const int & rightCellIndex)
{
    return (rightCellIndex < 0);
}

bool IsElementNotExist(const int & rightCellIndex)
{
    return (rightCellIndex < 0);
}
//! Copy end.

int GetNumberOfBasicElement()
{
    //const int PENTA_12 = NofValidElementTypes;
    //const int QUAD_6   = NofValidElementTypes+1;
    return (NofValidElementTypes + 2);
}

}
