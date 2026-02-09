/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */


// $Revision: 1.00 $
// $Date: 2026/01/20 $

#ifndef embeddedShell_h
#define embeddedShell_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <vector>

class SectionForceDeformation;
class ASDShellQ4Transformation;
class ASDShellQ4LocalCoordinateSystem;


class embeddedShell : public Element {
public:
	// life cycle
	embeddedShell();
	embeddedShell(
		int tag,
		int node1,
		int node2,
		int node3,
		int node4,
		SectionForceDeformation* section,
		bool corotational
	);
	~embeddedShell();

		const char* getClassType(void) const { return "embeddedShell"; };

		// domain methods
		void setDomain(Domain* theDomain);

		// print
		void Print(OPS_Stream& s, int flag);

		// methods dealing with nodes and number of external dof
		int getNumExternalNodes() const;
		const ID& getExternalNodes();
		Node** getNodePtrs(void);
		int getNumDOF(void);

		// methods dealing with committed state and update
		int commitState(void);
		int revertToLastCommit(void);
		int revertToStart(void);
		int update(void);

		// methods to return the current linearized stiffness and mass matrices
		const Matrix& getTangentStiff(void);
		const Matrix& getInitialStiff(void);
		const Matrix& getMass(void);

		// methods to addInertiaLoadToUnbalance -> do we need this? Dangerous
		int addInertiaLoadToUnbalance(const Vector& accel);

		const Vector& getResistingForce(void);
		const Vector& getResistingForceIncInertia(void);

		// public methods for element output
		int sendSelf(int commitTag, Channel& theChannel);
		int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

		Response* setResponse(const char** argv, int argc, OPS_Stream& output);
		int getResponse(int responseID, Information& eleInformation);
		
		// calculate characteristic length
		// double getCharacteristicLenght(void);

		// display
		int displaySelf(Renderer& theViewer, int displayMode, float fact, const char** displayModes = 0, int numModes = 0);


private:

	// internal methods to compute and store matrices and vectors
	int calculateAll(Matrix& LHS, Vector& RHS, int options);

private:

	// cross section
	SectionForceDeformation** theSections;

	// nodals ids
	ID m_node_ids = ID(4);
	Node* nodePointers[4] = { nullptr, nullptr, nullptr, nullptr };

	// coordinate transformation
	ASDShellQ4Transformation* m_transformation = nullptr;

	// section orientation with respect to local coordinate system
	Vector* m_local_x = nullptr;
	double m_angle;
	// initialization flag
	bool m_initialized = false;

	// input natural coordinates 
	double eta1, xi1;
	double eta2, xi2;


};



#endif
