/*
 * Copyright (c) 2019
 *	Side Effects Software Inc.  All rights reserved.
 *
 * Redistribution and use of Houdini Development Kit samples in source and
 * binary forms, with or without modification, are permitted provided that the
 * following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. The name of Side Effects Software may not be used to endorse or
 *    promote products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SIDE EFFECTS SOFTWARE `AS IS' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN
 * NO EVENT SHALL SIDE EFFECTS SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *----------------------------------------------------------------------------
 * The Subtribasic SOP.  This is provided as athe simplest possible example of the graded subdivision tringle mesh library.
 */

#include "SOP_SubtriBasic.h"

#include <GU/GU_Detail.h>
#include <GA/GA_Iterator.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_Vector3.h>
#include <SYS/SYS_Math.h>
#include <stddef.h>
#include <PRM/PRM_DialogScript.h>
#include <OP/OP_Caller.h>
#include <GU/GU_PrimPoly.h>


using namespace IDW_Tools;
using namespace subTri;

void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
        "subtribasic",
        "SubtriBasic",
        SOP_SubtriBasic::myConstructor,
        SOP_SubtriBasic::myTemplateList,
        0,
        0,
        0,
		OP_FLAG_GENERATOR));
}

static PRM_Name names[] = {
  PRM_Name("max_iter",	"Max Iter"),
  PRM_Name("angle",	"main angle"),
  PRM_Name("radius",	"sphere radius"),
  PRM_Name("resmult",      "resmult"),
};

static PRM_Default      maxIterDefault(10);
static PRM_Default      angleDefault(1);
static PRM_Default      radiusDefault(1.0);
static PRM_Default      resmultDefault(.01);

PRM_Template
SOP_SubtriBasic::myTemplateList[] = {
  PRM_Template(PRM_INT_J, 1, &names[0], &maxIterDefault),
  PRM_Template(PRM_FLT_J, 1, &names[1], &angleDefault),
  PRM_Template(PRM_FLT_J, 1, &names[2], &radiusDefault),
  PRM_Template(PRM_FLT_J, 1, &names[3], &resmultDefault),
  PRM_Template(),
};


OP_Node *
SOP_SubtriBasic::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
  cout<<"myConstructor"<<endl;
    return new SOP_SubtriBasic(net, name, op);
}

SOP_SubtriBasic::SOP_SubtriBasic(OP_Network *net, const char *name, OP_Operator *op)
    : SOP_Node(net, name, op)
{
    mySopFlags.setManagesDataIDs(true);
}

SOP_SubtriBasic::~SOP_SubtriBasic() {}

OP_ERROR
SOP_SubtriBasic::cookMySop(OP_Context &context)
{
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT)
        return error();
	// standard SOP setup
    fpreal now = context.getTime();
    OP_Node::flags().setTimeDep(true);
    setVariableOrder(3, 2, 0, 1);
    setCurGdh(0, myGdpHandle);
    setupLocalVars();
	gdp->clearAndDestroy();

	//get UI values
	int maxIter = MAXITER(now);
	float triradius = ANGLE(now);
	float sphRad = RADIUS(now);
	float resmult = RESMULT(now);

	//-----------------------------------------------------------------------------------
	// -------------- THIS IS WHERE WE CALL THE ACTUAL GRADED SUBDIV CODE ---------------
	subdividedMesh<float> sphere;
	sphere.initTet(sphRad);
	sphere.exampleSubdivide(maxIter, triradius, sphRad, resmult);
	//-----------------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------

	// the rest of the code just iterates over the mesh, building data in Houdini.
	// add point attribs to GDP
	GA_RWHandleIA pidxAttr = GA_RWHandleIA(gdp->addIntArray(GA_ATTRIB_POINT, "primindices"));

	// iterate over the subdiv mesh points, creating a list of Houdini points with attribs
	vector<smPoint<float>>::iterator sphPt;
	for (sphPt = sphere.pointList.begin(); sphPt != sphere.pointList.end(); ++sphPt){
	  UT_Vector3 newLoc;
	  newLoc.assign(sphPt->loc[0], sphPt->loc[1], sphPt->loc[2]); 
	  GA_Offset ptoff = gdp->appendPointOffset();
	  gdp->setPos3(ptoff, newLoc);
	  UT_Int32Array pidxArr;
	  pidxArr.setSize(exint(sphPt->trisIndxs.size()));
	  for(int i = 0; i < sphPt->trisIndxs.size(); ++i){
		pidxArr[i] =  sphPt->trisIndxs[i];
	  }
	  pidxAttr.set(ptoff, pidxArr);
	}

	// add prim attribs to GDP
	GA_RWHandleV3 neighborsAttr = GA_RWHandleV3(gdp->addFloatTuple(GA_ATTRIB_PRIMITIVE, "neighbors", 3));
	GA_RWHandleI indexAttr = GA_RWHandleI(gdp->addIntTuple(GA_ATTRIB_PRIMITIVE, "index", 1));
	GA_RWHandleI splitAttr = GA_RWHandleI(gdp->addIntTuple(GA_ATTRIB_PRIMITIVE, "split", 1));
	GA_RWHandleV3 childrenAttr = GA_RWHandleV3(gdp->addFloatTuple(GA_ATTRIB_PRIMITIVE, "children", 3));
	GA_RWHandleI cenaddrAttr = GA_RWHandleI(gdp->addIntTuple(GA_ATTRIB_PRIMITIVE, "cenAddr", 1));
	GA_RWHandleV3 vertsAttr = GA_RWHandleV3(gdp->addFloatTuple(GA_ATTRIB_PRIMITIVE, "verts", 3));
	GA_RWHandleV3 flippedAttr = GA_RWHandleV3(gdp->addFloatTuple(GA_ATTRIB_PRIMITIVE, "flipped", 3));
	
	// iterate over the subdiv mesh triangles, creating a list of Houdini prims with attribs
	vector<subtri<float>>::iterator tri;
	for (tri = sphere.subtriList.begin(); tri != sphere.subtriList.end(); ++tri){
	  if(tri->splitState != FINISHED){
		GU_PrimPoly *poly = GU_PrimPoly::build(gdp, 3,0,0);
		poly->setVertexPoint(0, tri->vertAddr(0));
		poly->setVertexPoint(1, tri->vertAddr(1));
		poly->setVertexPoint(2, tri->vertAddr(2));
		UT_Vector3 neighbors(tri->neighbors[0].get(),tri->neighbors[1].get(),tri->neighbors[2].get());
		neighborsAttr.set(poly->getMapOffset(), neighbors);
		indexAttr.set(poly->getMapOffset(), tri->index);
		splitAttr.set(poly->getMapOffset(), tri->splitState);
		UT_Vector3 children(tri->children[0].get(),tri->children[1].get(),tri->children[2].get());
		childrenAttr.set(poly->getMapOffset(), children);
		cenaddrAttr.set(poly->getMapOffset(), tri->cenAddr.get());
		UT_Vector3 verts(tri->verts[0].get(),tri->verts[1].get(),tri->verts[2].get());
		vertsAttr.set(poly->getMapOffset(), verts);
		UT_Vector3 flipped(tri->flipped[0],tri->flipped[1],tri->flipped[2]);
		flippedAttr.set(poly->getMapOffset(), flipped);
	  }
	}
    resetLocalVarRefs();
    return error();
}

const char *
SOP_SubtriBasic::inputLabel(unsigned) const
{
    return "no inputs to this SOP";
}
