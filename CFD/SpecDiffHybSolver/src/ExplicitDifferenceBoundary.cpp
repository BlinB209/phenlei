#include "GlobalDataBase.h"
#include "TK_Exit.h"
#include "ExplicitDifferenceBoundary.h"
#include "PHMpi.h"
#include <iostream>
#include <sstream>
#include <cstdlib>
#include "PHHeader.h"

using namespace std;

namespace PHSPACE
{
	ExplicitDifferenceBoundary::ExplicitDifferenceBoundary()
	{
	}

	void ExplicitDifferenceBoundary::GetSchemeCoefDiffBoundaryData()
	{
		double rb0, rb1, rb2, rb3, rb4;

		nBoundary_Neumann = 3;

		Range I(0, nBoundary_Neumann + 1);
		boundary_Neumann_0 = new RDouble1D(I, fortranArray);

		Range J(-(nBoundary_Neumann + 1), 0);
		boundary_Neumann_N = new RDouble1D(J, fortranArray);

		int ibcscheme = GlobalDataBase::GetIntParaFromDB("ibcscheme");
		if (ibcscheme == 3)
		{
			rb0 = -11.0/6.0 ; rb1 = 3.0; rb2 = -3.0/2.0; rb3 = 1.0/3.0; rb4 = 0.0     ;
		}
		else if (ibcscheme == 4)
		{
			rb0 = -25.0/12.0; rb1 = 4.0; rb2 = -3.0    ; rb3 = 4.0/3.0; rb4 = -1.0/4.0;
		}
		else
		{
			TK_Exit::UnexpectedVarValue( "ibcscheme", ibcscheme );
		}

		(*boundary_Neumann_0)(0) = rb0;
		(*boundary_Neumann_0)(1) = rb1;
		(*boundary_Neumann_0)(2) = rb2;
		(*boundary_Neumann_0)(3) = rb3;
		(*boundary_Neumann_0)(4) = rb4;

		(*boundary_Neumann_N)( 0) = -rb0;
		(*boundary_Neumann_N)(-1) = -rb1;
		(*boundary_Neumann_N)(-2) = -rb2;
		(*boundary_Neumann_N)(-3) = -rb3;
		(*boundary_Neumann_N)(-4) = -rb4;
	}

	ExplicitDifferenceBoundary::~ExplicitDifferenceBoundary()
	{
		delete boundary_Neumann_0; boundary_Neumann_0 = NULL;
        delete boundary_Neumann_N; boundary_Neumann_N = NULL;
	}
}

