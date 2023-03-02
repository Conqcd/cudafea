//
//  IntegrationMatrix.cpp
//
//

#include "IntegrationMatrix.h"

// not a petsc matrix !

IntegrationMatrix::IntegrationMatrix(int _type = 0): type(_type) // a(1), b(1) c(1)
{
    create();
}

IntegrationMatrix::~IntegrationMatrix()
{
	
}

void IntegrationMatrix::create()
{
    // Create the Integration Matrix
    //matrix.resize(SAMPLES_PER_ELEMENT);
    switch (type)
	{
    case 0:
	    {
	        matrix[0][0] = -1; matrix[0][1] = -1; matrix[0][2] = -1;
	        matrix[1][0] = +1; matrix[1][1] = -1; matrix[1][2] = -1;
	        matrix[2][0] = -1; matrix[2][1] = +1; matrix[2][2] = -1;
	        matrix[3][0] = +1; matrix[3][1] = +1; matrix[3][2] = -1;
	        matrix[4][0] = -1; matrix[4][1] = -1; matrix[4][2] = +1;
	        matrix[5][0] = +1; matrix[5][1] = -1; matrix[5][2] = +1;
	        matrix[6][0] = -1; matrix[6][1] = +1; matrix[6][2] = +1;
	        matrix[7][0] = +1; matrix[7][1] = +1; matrix[7][2] = +1;

	        for (int j = 0; j < DOF_3D; j++)
	        {
	            for (int i = 0; i < SAMPLES_PER_ELEMENT; i++)
	            {
	                // 8 integration points
	                matrix[i][j] = matrix[i][j] / (2 * sqrt(3)); // R 
	            }
	        }
    		break;
	    }
	case 1:
	{
		matrix[0][0] = -1; matrix[0][1] = -1; matrix[0][2] = -1;
		matrix[1][0] = +1; matrix[1][1] = -1; matrix[1][2] = -1;
		matrix[2][0] = -1; matrix[2][1] = +1; matrix[2][2] = -1;
		matrix[3][0] = +1; matrix[3][1] = +1; matrix[3][2] = -1;
		matrix[4][0] = -1; matrix[4][1] = -1; matrix[4][2] = +1;
		matrix[5][0] = +1; matrix[5][1] = -1; matrix[5][2] = +1;
		matrix[6][0] = -1; matrix[6][1] = +1; matrix[6][2] = +1;
		matrix[7][0] = +1; matrix[7][1] = +1; matrix[7][2] = +1;

		for (int j = 0; j < DOF_3D; j++)
		{
			for (int i = 0; i < SAMPLES_PER_ELEMENT; i++)
			{
				// 8 integration points
				matrix[i][j] = matrix[i][j] / (2 * sqrt(3)); // R 
			}
		}
		break;
	}
    default:
		break;
	}
}

Diff_Ni_j::Diff_Ni_j()
{
    create();
}

Diff_Ni_j::~Diff_Ni_j()
{
	
}

void Diff_Ni_j::create()
{
    IntegrationMatrix temp(0);

    //matrix.resize(SAMPLES_PER_ELEMENT);
	
    for (int k = 0; k < SAMPLES_PER_ELEMENT; k++)
    {
        for (int j = 0; j < NODES_PER_ELEMENT; j++)
        {
            for (int i = 0; i < DOF_3D; i++)
            {
                int id1, id2;
                id1 = (i + 1) % DOF_3D;
                id2 = (id1 + 1) % DOF_3D;
                matrix[k][i][j] = SIGN(temp(j, i)) * (0.5 + SIGN(temp(j, id1)) * temp(k, id1)) * (0.5 + SIGN(temp(j, id2)) * temp(k, id2)) / 2;
            }
        }
    }
}