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


// Written: Stefano Ercolessi
// Created: 22/01/2026

// Purpose: This file contains the class definition for 
// embedded Timoshenko beam element in a MITC4 shell with Drilling DOF


#include <embeddedShell.h>
#include <ASDShellQ4CorotationalTransformation.h>

#include <SectionForceDeformation.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <ElementResponse.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <Renderer.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sstream>



namespace {

	// calculation options
	constexpr int OPT_NONE = 0;
	constexpr int OPT_UPDATE = (1 << 0);
	constexpr int OPT_LHS = (1 << 1);
	constexpr int OPT_RHS = (1 << 2);
	constexpr int OPT_LHS_IS_INITIAL = (1 << 3);
	
	// this should be an input !!!!!!!!!
	inline double computeBeamLength(const double xi1, const double eta1, const double xi2, const double eta2) {	
		// length in the natural coordinate system
		const double dxi = xi2 - xi1;
		const double deta = eta2 - eta1;
		//return std::sqrt(dxi * dxi + deta * deta);
		return 1.0;
	}

	inline double  computeBeamLenghtFromLocalCS(const ASDShellQ4LocalCoordinateSystem& local_cs,
		const double xi1, const double eta1, const double xi2, const double eta2) {
		// get the coordinates of the two beam nodes in the local CS
		auto ShapeQ4 = [](double xi, double eta, double N[4]) {
			N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
			N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
			N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
			N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
			};
		double N1[4], N2[4];
		ShapeQ4(xi1, eta1, N1);
		ShapeQ4(xi2, eta2, N2);

		const auto& P = local_cs.Nodes();

		double P1[3] = { 0.0, 0.0, 0.0 };
		double P2[3] = { 0.0, 0.0, 0.0 };
		for (int i = 0; i < 4; i++) {
			P1[0] += N1[i] * P[i](0);
			P1[1] += N1[i] * P[i](1);
			P1[2] += N1[i] * P[i](2);
			P2[0] += N2[i] * P[i](0);
			P2[1] += N2[i] * P[i](1);
			P2[2] += N2[i] * P[i](2);
		}

		const double dx = P2[0] - P1[0];
		const double dy = P2[1] - P1[1];
		const double dz = P2[2] - P1[2];

		return std::sqrt(dx * dx + dy * dy + dz * dz);
	}

	// Compute Timoshenko beam shape functions
	inline Matrix computeBeamShapeFunctions(const double xi) {
		// shape function matrix
		Matrix N(6, 12);
		N.Zero();
		// first node
		N(0, 0) = N(1,1) = N(2,2) = N(3,3) = N(4,4) = N(5,5) = 0.5 * (1.0 - xi);
		N(6, 0) = N(7, 1) = N(8, 2) = N(9, 3) = N(10, 4) = N(11, 5) = 0.5 * (1.0 + xi);
	}

	inline Matrix computeBeamShapeFunctionDerivatives(const double xi) {
		// shape function derivative matrix
		Matrix dN(6, 12);
		dN.Zero();
		// first node
		dN(0, 0) = dN(1,1) = dN(2,2) = dN(3,3) = dN(4,4) = dN(5,5) = -0.5;
		dN(6, 0) = dN(7, 1) = dN(8, 2) = dN(9, 3) = dN(10, 4) = dN(11, 5) = 0.5;
	}

	inline Matrix computeStrainMatrix(const double xi, const double L) {
		// strain displacment matrix
		Matrix B(6, 12);
		B.Zero();
		// they are the same do not call the functions
		double N1 = 0.5;
		double N2 = 0.5;
		double dN1 = -1.0/L;
		double dN2 = 1.0/L;
		// axial strain
		B(0, 0) = dN1; B(0, 6) = dN2;
		// torsional curvature
		B(1, 3) = dN1; B(1, 9) = dN2;
		// bending curvature kappa_y
		B(2, 4) = dN2; B(2, 10) = dN1;
		// bending curvature kappa_z
		B(3, 5) = dN1; B(3, 11) = dN2;
		// shear strain gamma_y
		B(4, 1) = dN1; B(4, 5) = -N1; B(4, 7) = dN2; B(4, 11) = -N2;
		// shear strain gamma_z
		B(5, 2) = dN1; B(5, 4) = N1;  B(5, 8) = dN2; B(5, 10) = N2;

		return B;
	}

	// Compute the rotation matrix from local shell to local beam
	inline Matrix computeShellToBeamRotationMatrix(const double alpha) {
		// compute cos and sin
		const double c = std::cos(alpha);
		const double s = std::sin(alpha);
		// build the rotation matrix
		Matrix R(12, 12);
		R.Zero();
		// first Node
		R(0, 0) = R(1, 1) = c;
		R(0, 1) = s; R(1, 0) = -s;
		R(2, 2) = 1.0;
		R(3, 3) = R(4, 4) = R(5, 5) = 1.0;
		// second Node
		R(6, 6) = R(7, 7) = c;
		R(6, 7) = s; R(7, 6) = -s;
		R(8, 8) = 1.0;	
		R(9, 9) = R(10, 10) = R(11, 11) = 1.0;
		// set the rotation matrix
		return R;
	};

	// Compute the transformation matrix from local shell to beam element
	// MPC constraints Transformation matrix
	inline Matrix computeShellToBeamTransformationMatrix(const double xi1, const double eta1,
		const double xi2, const double eta2) {
		// build the transformation matrix
		Matrix T(12, 24);
		T.Zero();
		// Compute the shape functions at the two beam nodes
		auto shapeQ4 = [](double xi, double eta, double N[4]) {
			N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
			N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
			N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
			N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
			};
		double N1[4], N2[4];
		shapeQ4(xi1, eta1, N1);
		shapeQ4(xi2, eta2, N2);
		// first beam node

		T(0, 0) = T(1, 1) = T(2, 2) = T(3, 3) = T(4, 4) = T(5, 5) = N1[0];
		T(0, 6) = T(1, 7) = T(2, 8) = T(3, 9) = T(4, 10) = T(5, 11) = N1[1];
		T(0, 12) = T(1, 13) = T(2, 14) = T(3, 15) = T(4, 16) = T(5, 17) = N1[2];
		T(0, 18) = T(1, 19) = T(2, 20) = T(3, 21) = T(4, 22) = T(5, 23) = N1[3];

		T(6, 0) = T(7, 1) = T(8, 2) = T(9, 3) = T(10, 4) = T(11, 5) = N2[0];
		T(6, 6) = T(7, 7) = T(8, 8) = T(9, 9) = T(10, 10) = T(11, 11) = N2[1];
		T(6, 12) = T(7, 13) = T(8, 14) = T(9, 15) = T(10, 16) = T(11, 17) = N2[2];
		T(6, 18) = T(7, 19) = T(8, 20) = T(9, 21) = T(10, 22) = T(11, 23) = N2[3];

		// return the transformation matrix
		return T;
	};

	// Some usefull functions
	inline bool computeProjectedBeamAngle(Node* nodePointers[4], double xi1, double eta1,
		double xi2, double eta2, double& outAngle, double outDir[3]) {
		constexpr double small = 1.0e-12;

		// Q4 shape functions
		auto shapeQ4 = [](double xi, double eta, double N[4]) {
			N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
			N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
			N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
			N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
			};

		// map natural to global
		auto map2Global = [&](double xi, double eta, double P[3]) -> bool {
			double N[4];
			shapeQ4(xi, eta, N);
			P[0] = P[1] = P[2] = 0.0;
			for (int i = 0; i < 4; i++) {
				if (nodePointers[i] == nullptr) return false;
				const Vector& coord = nodePointers[i]->getCrds();
				P[0] += N[i] * coord(0);
				P[1] += N[i] * coord(1);
				P[2] += N[i] * coord(2);
			}
			return true;
			};

		double P1[3], P2[3];
		if (!map2Global(xi1, eta1, P1) || !map2Global(xi2, eta2, P2))
			return false;

		// element node coordinates
		double PN[4][3];
		for (int i = 0; i < 4; i++) {
			const Vector& crd = nodePointers[i]->getCrds();
			PN[i][0] = crd(0);
			PN[i][1] = crd(1);
			PN[i][2] = crd(2);
		}

		// reference X (same default used in ASDShellQ4)
		double x_ref[3] = {
			0.5 * (PN[1][0] + PN[2][0]) - 0.5 * (PN[0][0] + PN[3][0]),
			0.5 * (PN[1][1] + PN[2][1]) - 0.5 * (PN[0][1] + PN[3][1]),
			0.5 * (PN[1][2] + PN[2][2]) - 0.5 * (PN[0][2] + PN[3][2])
		};
		double xlen = std::sqrt(x_ref[0]*x_ref[0] + x_ref[1]*x_ref[1] + x_ref[2]*x_ref[2]);
		if (xlen < small) return false;
		x_ref[0]/=xlen; x_ref[1]/=xlen; x_ref[2]/=xlen;

		// approximate shell normal (use two edges)
		double vA[3] = { PN[1][0] - PN[0][0], PN[1][1] - PN[0][1], PN[1][2] - PN[0][2] };
		double vB[3] = { PN[3][0] - PN[0][0], PN[3][1] - PN[0][1], PN[3][2] - PN[0][2] };
		double n_shell[3] = {
			vA[1] * vB[2] - vA[2] * vB[1],
			vA[2] * vB[0] - vA[0] * vB[2],
			vA[0] * vB[1] - vA[1] * vB[0]
		};
		double nlen = std::sqrt(n_shell[0]*n_shell[0] + n_shell[1]*n_shell[1] + n_shell[2]*n_shell[2]);
		if (nlen < small) return false;
		n_shell[0]/=nlen; n_shell[1]/=nlen; n_shell[2]/=nlen;

		// beam vector and projection on shell plane
		double v[3] = { P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2] };
		double vdotn = v[0]*n_shell[0] + v[1]*n_shell[1] + v[2]*n_shell[2];
		double v_proj[3] = { v[0] - vdotn*n_shell[0], v[1] - vdotn*n_shell[1], v[2] - vdotn*n_shell[2] };
		double vplen = std::sqrt(v_proj[0]*v_proj[0] + v_proj[1]*v_proj[1] + v_proj[2]*v_proj[2]);

		if (vplen < small) {
			// degenerate: fallback to reference axis
			outDir[0] = x_ref[0];
			outDir[1] = x_ref[1];
			outDir[2] = x_ref[2];
			outAngle = 0.0;
			return true;
		}

		v_proj[0] /= vplen; v_proj[1] /= vplen; v_proj[2] /= vplen;

		// signed angle between reference axis and beam projection
		double dot_ev = x_ref[0]*v_proj[0] + x_ref[1]*v_proj[1] + x_ref[2]*v_proj[2];
		dot_ev = std::max(-1.0, std::min(1.0, dot_ev));
		double cross_ev[3] = {
			x_ref[1]*v_proj[2] - x_ref[2]*v_proj[1],
			x_ref[2]*v_proj[0] - x_ref[0]*v_proj[2],
			x_ref[0]*v_proj[1] - x_ref[1]*v_proj[0]
		};
		double cross_dot_n = cross_ev[0]*n_shell[0] + cross_ev[1]*n_shell[1] + cross_ev[2]*n_shell[2];
		outAngle = std::atan2(cross_dot_n, dot_ev);

		outDir[0] = v_proj[0]; outDir[1] = v_proj[1]; outDir[2] = v_proj[2];

		return true;
	}

	inline bool computeBeamLocalAngle(const double xi1, const double etai1, const double xi2, const double etai2,
		double& outAngle, double outDir[3]) {
		// compute m_angle which is the angle between beam and the local x of the shell
		outAngle = std::atan2( etai2 - etai1, xi2 - xi1 );
		outDir[0] = std::cos(outAngle);
		outDir[1] = std::sin(outAngle);
		outDir[2] = 0.;

		return true;
	}

	

	class Globals {
	private:
		Globals() = default;
		Globals(const Globals&) = delete;
		Globals& operator = (const Globals&) = delete;

	public:
		static Globals& instance() {
			static Globals instance;
			return instance;
		}

	public:
		// beam lenght
		double L = 0.0;
		// global displacements
		Vector UG = Vector(24);
		// local displacments
		Vector UL = Vector(24);
		// Vector beam local displacments
		Vector dBloc = Vector(12);

		// B matrix
		Matrix B = Matrix(6, 12);
		// N shape function matrix timoshenko beam
		Matrix N = Matrix(6, 12);
		// Material local matrix
		Matrix D = Matrix(6, 6);

		// Strain local vector
		Vector Eps = Vector(6);
		// Stress local vector
		Vector Sig = Vector(6);
		// Vector local RHS -> return from the beam section 
		Vector Rlocal = Vector(6);

		// T' = R T
		Matrix Tprime = Matrix(12, 24);

		// beam local stiffness matrix
		Matrix Kbeam = Matrix(12, 12);

		// LHS matrix tangent stiffness matrix
		Matrix LHS = Matrix(24, 24);
		// LHS initial tangent stiffness matrix
		Matrix LHS_initial = Matrix(24, 24);
		// LHS mass matrix
		Matrix LHS_mass = Matrix(24, 24);
		// RHS rsidual vector
		Vector RHS = Vector(24);
		// RHS residual vector with inertia terms
		Vector RHS_winertia = Vector(24);

	};
}

static void debugPrintVector(const Vector& v, const char* name, int maxDump = 8) {
	if (&v == nullptr) return;
	opserr << name << " size=" << v.Size() << " norm=" << v.Norm() << "\n";
	// check finite and print first entries
	for (int i = 0; i < v.Size() && i < maxDump; ++i) {
		double val = v(i);
		if (!std::isfinite(val)) opserr << "  " << name << "(" << i << ") = " << val << " [!finite]\n";
		else opserr << "  " << name << "(" << i << ") = " << val << "\n";
	}
	if (v.Size() > maxDump) opserr << "  ...\n";
}
static void debugPrintMatrix(const Matrix& M, const char* name, int maxDump = 6) {
	opserr << name << " size=" << M.noRows() << "x" << M.noCols() << "\n";
	// Frobenius-like check
	double maxv = 0.0, minv = 0.0;
	for (int i = 0; i < M.noRows(); ++i) for (int j = 0; j < M.noCols(); ++j) {
		double val = M(i, j);
		if (!std::isfinite(val)) { opserr << "  " << name << "(" << i << "," << j << ") = " << val << " [!finite]\n"; return; }
		maxv = std::max(maxv, fabs(val));
		minv = std::min(minv, -fabs(val));
	}
	opserr << "  max abs entry ~ " << maxv << "\n";
	// print small block and diag
	int r = std::min(maxDump, M.noRows());
	int c = std::min(maxDump, M.noCols());
	for (int i = 0; i < r; ++i) {
		std::ostringstream line;
		for (int j = 0; j < c; ++j) line << M(i, j) << " ";
		opserr << "  " << name << "[" << i << ",0.." << (c - 1) << "] = " << line.str().c_str() << "\n";
	}
	if (M.noRows() > maxDump || M.noCols() > maxDump) opserr << "  ...\n";
}

// define the OPS_ function
void*
OPS_embeddedShell(void)
{
	static bool first_done = false;
	if (!first_done) {
		opserr << "Using embeddedShell - Developed by: ASDEA\n";
		first_done = true;
	}

	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 6) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: element embeddedShell eleTag node1 node2 node3 node4 sectionTag\n";
		return 0;
	}

	// get the input data INT
	int iData[6];
	int numData = 6;
	if (OPS_GetInt(&numData, iData) != 0) {
		opserr << "WARNING invalid integer tag: element embeddedShell \n";
		return 0;
	}

	bool corotational = false;

	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* type = OPS_GetString();
		if (strcmp(type, "-corotational") == 0 || (strcmp(type, "-Corotational") == 0)) {
			corotational = true;
		}
		else {
			OPS_ResetCurrentInputArg(-1);
			break;
		}
	}

	SectionForceDeformation* section = OPS_getSectionForceDeformation(iData[5]);

	if (section == 0) {
		opserr << "ERROR element embeddedShell " << iData[0] << " section " << iData[5] << " not found\n";
		return 0;
	}

	return new embeddedShell(
		iData[0],
		iData[1],
		iData[2],
		iData[3],
		iData[4],
		section,
		corotational
	);


}

// define empty constructor
embeddedShell::embeddedShell()
	: Element(0, ELE_TAG_embeddedShell)
{
}

// define the full constructor
embeddedShell::embeddedShell(
	int tag,
	int node1,
	int node2,
	int node3,
	int node4,
	SectionForceDeformation* section,
	bool corotational
)
	: Element(tag, ELE_TAG_embeddedShell),
	m_transformation(corotational ? new ASDShellQ4CorotationalTransformation() : new ASDShellQ4Transformation())
{
	m_node_ids(4);
	m_node_ids(0) = node1;
	m_node_ids(1) = node2;
	m_node_ids(2) = node3;
	m_node_ids(3) = node4;

	theSections = new SectionForceDeformation*[1];
	theSections[0] = section;

	// initialize other member variables
	for (int i = 0; i < 4; i++) {
		nodePointers[i] = nullptr;
	}
	

	xi1 = -1.0;
	xi2 = 1.0;
	eta1 = 0.0;
	eta2 = 0.0;

}

embeddedShell::~embeddedShell()
{
	// clean up section
	delete theSections[0];
	delete[] theSections;
}

// set the domain
void
embeddedShell::setDomain(Domain* theDomain)
{
	// if domain is null
	if (theDomain == nullptr) {
		for (int i = 0; i < 4; i++)
			nodePointers[i] = nullptr;

	    // set transformation ????????????????????
		m_transformation->setDomain(theDomain, m_node_ids, m_initialized);
		//call base class implementation
		DomainComponent::setDomain(theDomain);
	}

	// node pointers
	for (int i = 0; i < 4; i++) {
		nodePointers[i] = theDomain->getNode(m_node_ids(i));
	}

	// set domain on transformation
	m_transformation->setDomain(theDomain, m_node_ids, m_initialized);

	// if is not initialized, initialize the local coordinate system
	if (!m_initialized) {
        // compute projected beam angle and store direction in m_local_x
        double beamDir[3] = {0.0, 0.0, 0.0};
        if (!computeBeamLocalAngle(xi1, eta1, xi2, eta2, m_angle, beamDir)) {
            opserr << "embeddedShell::setDomain - computeProjectedBeamAngle failed for element " << this->getTag() << endln;
            // fallback: set default local x to element reference X (use node coords)
            const Vector& c0 = nodePointers[0]->getCrds();
            const Vector& c1 = nodePointers[1]->getCrds();
            beamDir[0] = c1(0) - c0(0);
            beamDir[1] = c1(1) - c0(1);
            beamDir[2] = c1(2) - c0(2);
            double blen = std::sqrt(beamDir[0]*beamDir[0] + beamDir[1]*beamDir[1] + beamDir[2]*beamDir[2]);
            if (blen > 0.0) { beamDir[0]/=blen; beamDir[1]/=blen; beamDir[2]/=blen; }
            m_angle = 0.0;
        }

        // store in m_local_x (allocate if needed) — keep same semantics as ASDShellQ4 expects (normalized)
        if (m_local_x == nullptr)
            m_local_x = new Vector(3);
        (*m_local_x)(0) = beamDir[0];
        (*m_local_x)(1) = beamDir[1];
        (*m_local_x)(2) = beamDir[2];
		//opserr << *m_local_x;
		
        m_initialized = true;
    }

    // call base class implementation
    DomainComponent::setDomain(theDomain);
}

int embeddedShell::getNumExternalNodes() const
{
	return 4;
}

const ID& embeddedShell::getExternalNodes()
{
	return m_node_ids;
}

int embeddedShell::getNumDOF()
{
	return 24;
}

Node** embeddedShell::getNodePtrs()
{
	return m_transformation->getNodes().data();
}

// to be filled 
int embeddedShell::sendSelf(int commitTag, Channel& theChannel)
{
	return 0;
}

int embeddedShell::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	return 0;
}	

void embeddedShell::Print(OPS_Stream& s, int flag)
{
	
}

Response*embeddedShell::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	return nullptr;
}

int embeddedShell::getResponse(int responseID, Information& eleInfo)
{
	return 0;
}

int embeddedShell::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** displayModes, int numModes)
{
	return 0;
}

// methods dealing with the commit state and update
int embeddedShell::commitState() {
	// commitState of the material
	int retVal = 0;
	retVal += theSections[0]->commitState();
	return retVal;
}

int embeddedShell::revertToLastCommit() {
	// revertToLastCommit of the material
	int retVal = 0;
	retVal += theSections[0]->revertToLastCommit();
	return retVal;
}

int embeddedShell::revertToStart() {
	// revertToStart of the material
	int retVal = 0;
	retVal += theSections[0]->revertToStart();
	return retVal;
}

int embeddedShell::update() {
	// compute
 	auto& LHS = Globals::instance().LHS;
	auto& RHS = Globals::instance().RHS;
	return calculateAll(LHS, RHS, OPT_UPDATE);
}

const Matrix& embeddedShell::getTangentStiff()
{
	// compute
	auto& LHS = Globals::instance().LHS;
	auto& RHS = Globals::instance().RHS;
	calculateAll(LHS, RHS, (OPT_LHS));
	return LHS;
}

const Matrix& embeddedShell::getInitialStiff()
{
	// compute
	auto& LHS = Globals::instance().LHS;
	auto& RHS = Globals::instance().RHS;
	calculateAll(LHS, RHS, (OPT_LHS | OPT_LHS_IS_INITIAL));
	return LHS;
}

const Matrix& embeddedShell::getMass() {
	// Output matrix
	auto& LHS = Globals::instance().LHS_mass;
	auto& beamLength = Globals::instance().L;
	LHS.Zero();

	// Jacobian 
	double xi = 0.0;
	double w = 4.0;
	double jacobian = beamLength / 2.0;

	// we need to get the rho and A of the section from the section fiber
	double rho = theSections[0]->getRho();

	// compute local mass matrix
	double m = 0.5 * rho * jacobian * w;
	// consistent mass matrix for the beam element
	Matrix Mloc(12, 12);
	Mloc.Zero();
	// axial
	Mloc(0, 0) = Mloc(6, 6) = m;
	Mloc(0, 6) = Mloc(6, 0) = m / 2.0;
	// torsion
	Mloc(3, 3) = Mloc(9, 9) = m;
	Mloc(3, 9) = Mloc(9, 3) = m / 2.0;
	// bending y
	Mloc(4, 4) = Mloc(10, 10) = m;
	Mloc(4, 10) = Mloc(10, 4) = m / 2.0;
	// bending z
	Mloc(5, 5) = Mloc(11, 11) = m;
	Mloc(5, 11) = Mloc(11, 5) = m / 2.0;
	// shear y
	Mloc(1, 1) = Mloc(7, 7) = m;
	Mloc(1, 7) = Mloc(7, 1) = m / 2.0;
	// shear z
	Mloc(2, 2) = Mloc(8, 8) = m;
	Mloc(2, 8) = Mloc(8, 2) = m / 2.0;

	// transform to shell local system
	auto& Tprime = Globals::instance().Tprime;
	Tprime.addMatrixProduct(0.0, computeShellToBeamRotationMatrix(m_angle), computeShellToBeamTransformationMatrix(xi1, eta1, xi2, eta2), 1.0);
	// LHS mass
	LHS.addMatrixTripleProduct(0.0, Tprime, Mloc, 1.0);
	
	return LHS;
	
}

int embeddedShell::calculateAll(Matrix& LHS, Vector& RHS, int options)
{
	// Check th options
	if (!m_transformation->isLinear()) {
		// corotational alculatioin of the tangent LHS requires the RHS
		if(options & OPT_LHS)
			options |= OPT_RHS;
	}

	// initialize result
	int result = 0;
	// Zero the output
	if (options & OPT_LHS)
		LHS.Zero();
	if (options & OPT_RHS)
		RHS.Zero();

	// Global displacements
	auto& UG = Globals::instance().UG;
	auto& beamLength = Globals::instance().L;
	m_transformation->computeGlobalDisplacements(UG);

	// update transformation
	if(options & OPT_UPDATE)
		m_transformation->update(UG);

	// compute the reference coordinate system if the shell
	ASDShellQ4LocalCoordinateSystem cs = m_transformation->createReferenceCoordinateSystem();

	// compute the local coordinate system of the shell
	ASDShellQ4LocalCoordinateSystem local_cs = m_transformation->createLocalCoordinateSystem(UG);

	// set up all needed parameters fot the embedded shell

	// local displacmeents
	auto& UL = Globals::instance().UL;
	m_transformation->calculateLocalDisplacements(local_cs, UG, UL);

	// Doing the integation stuff
	// Reduced integration for shear locking
	// 1 Gauss point at the center

	double xi = 0.0;
	double w = 2.0;
	beamLength = computeBeamLenghtFromLocalCS(local_cs, xi1, eta1, xi2, eta2);
	double jacobian = beamLength / 2.0;

	// compute and integrate at the gauss point
	// Compute the strain matrix
	auto& B = Globals::instance().B;
	auto& D = Globals::instance().D;
	auto& dbeam = Globals::instance().dBloc;
	auto& strain = Globals::instance().Eps;
	auto& sigma = Globals::instance().Sig;
	auto& Tprime = Globals::instance().Tprime;
	auto& RHS_local = Globals::instance().Rlocal;

	B = computeStrainMatrix(xi, beamLength);

#if 0
	std::cout << "Print B matrix: for the Beam: \n";
	debugPrintMatrix(B, "B", 12);
#endif

	// Update Strain
	if (options & OPT_UPDATE) {
		// section deformation
		// compute the deformation vector to send to the fiber section
		// define the local beam displacements from the R and T
		// dbeam = R T dshell
		Tprime.addMatrixProduct(0.0, computeShellToBeamRotationMatrix(m_angle), computeShellToBeamTransformationMatrix(xi1, eta1, xi2, eta2), 1.0);

		dbeam.addMatrixVector(0.0, Tprime, UL, 1.0);		

		strain.addMatrixVector(0.0, B, dbeam, 1.0);
		
		// Update results in the section 
#if 0
		std::cout << "Print Rotation Matrix: \n";
		debugPrintMatrix(computeShellToBeamRotationMatrix(m_angle), "R", 12);
		std::cout << "Print Transformation Matrix: \n";
		debugPrintMatrix(computeShellToBeamTransformationMatrix(xi1, eta1, xi2, eta2), "T", 24);
		std::cout << "Print RT: \n";
		debugPrintMatrix(Tprime, "RT", 12);
		debugPrintVector(dbeam, "dbeam", 12);
		debugPrintVector(strain, "strain", 12);
#endif

		// send strain to material
		result += theSections[0]->setTrialSectionDeformation(strain);
	}

	// Integrate RHS
	if(options & OPT_RHS){

		// get the stress from the section		
		sigma = theSections[0]->getStressResultant();
		// transform to shell local system
		Tprime.addMatrixProduct(0.0, computeShellToBeamRotationMatrix(m_angle), computeShellToBeamTransformationMatrix(xi1, eta1, xi2, eta2), 1.0);
		// RHS local
		RHS_local.addMatrixTransposeVector(0.0, B, sigma, 1.0);
		
	}
	
	// Integrate RHS
	RHS.addMatrixTransposeVector(0.0, Tprime, RHS_local, w * jacobian);
	//opserr << "Section: \n" << theSections[0]->getSectionTangent() ;


	if (options & OPT_RHS || options & OPT_LHS) {

		// Section tangent stiffness or Initial?
		const Matrix & Csec = (options & OPT_LHS_IS_INITIAL) ?
			theSections[0]->getInitialTangent() :
			theSections[0]->getSectionTangent();

		// Obtain section type ID array so we know which entry corresponds to which
		const ID& secCode = theSections[0]->getType();
		int secOrder = secCode.Size();

		// helper to map SECTION_RESPONSE to a beam index
		auto secResponseToBeamIndex = [](int secResp) -> int {
			switch (secResp) {
			case SECTION_RESPONSE_P: return 0; // axial
			case SECTION_RESPONSE_T: return 1; // torsion
			case SECTION_RESPONSE_MY: return 2; // bending y curvature about y
			case SECTION_RESPONSE_MZ: return 3; // bending z curvature about z
			case SECTION_RESPONSE_VY: return 4; // shear force in y
			case SECTION_RESPONSE_VZ: return 5; // shear force in z
			default: return -1;
			}
		};

		// Im thinking that we need to force the user to use all the cases for the section
		// check 
		int nsec = Csec.noRows();
		if (nsec != Csec.noCols()) {
			opserr << "embeddedShell::calculaAll - section tangent not square for element ";
		}
		// build the D matrix

		D.Zero();

		for (int i = 0; i < secOrder; ++i) {
			int secRespI = secCode(i);
			int bi = secResponseToBeamIndex(secRespI);
			if (bi < 0) continue; // unknown quantity for beam

			for (int j = 0; j < secOrder; ++j) {
				int secRespJ = secCode(j);
				int bj = secResponseToBeamIndex(secRespJ);
				if (bj < 0) continue;

				// protect against sections that return fewer rows/cols than getOrder()
				if (i < Csec.noRows() && j < Csec.noCols())
					D(bi, bj) = Csec(i, j);
			}
		}

#if 0
		std::cout << "Matrix D: \n";
		debugPrintMatrix(D, "D", 8);
#endif

	}

	if (options & OPT_LHS) {
		// compute the beams tiffness matrix
		Matrix Klocal(12, 12);
		Klocal.Zero();
		Klocal.addMatrixTripleProduct(0.0, B, D, B, w * jacobian);

#if 0
		//Klocal.addMatrixTripleProduct(0.0, B, D, w * jacobian);
		std::cout << "klocal: \n";
		debugPrintMatrix(Klocal, "klocal", 12);
		//debugPrintMatrix(Tprime, "RT", 12);
#endif

		LHS.addMatrixTripleProduct(0.0, Tprime, Klocal, Tprime, 1.0);

#if 0
		std::cout << "Print Tprime matrix: for the Beam: \n";
		debugPrintMatrix(Tprime, "Tprime", 24);
		std::cout << "Print Kglobal matrix: for the Beam: \n";
		debugPrintMatrix(LHS, "Kglobal", 12);
#endif
	}

	m_transformation->transformToGlobal(local_cs, UG, UL, LHS, RHS, (options & OPT_LHS));

#if 0
	std::cout << "Print global displacements vector: for the Beam: \n";
	debugPrintVector(UG, "UG", 24);
	std::cout << "Print local displacements vector: for the Beam: \n";
	debugPrintVector(UL, "UL", 24);

	std::cout << "Print LHS matrix: for the Beam: \n";
	debugPrintMatrix(LHS, "LHS", 24);
	std::cout << "Print RHS vector: for the Beam: \n";
	debugPrintVector(RHS, "RHS", 24);

#endif
	return 0;
};

// methods for obtaining resisting force (it does not have elemental loads)
const Vector& embeddedShell::getResistingForce()
{
	// compute
	auto& LHS = Globals::instance().LHS;
	auto& RHS = Globals::instance().RHS;
	calculateAll(LHS, RHS, (OPT_RHS));

#if 0
	std::cout << "Print RHS vector: for the Beam: \n";
	debugPrintVector(RHS, "RHS", 24);
#endif

	// set the output
	return RHS;
}

const Vector& embeddedShell::getResistingForceIncInertia()
{
	// compute
	auto& LHS = Globals::instance().LHS;
	auto& RHS = Globals::instance().RHS_winertia;
	auto& Tprime = Globals::instance().Tprime;
	calculateAll(LHS, RHS, (OPT_RHS));
	// add inertia terms
	Vector accel(24);
	accel.Zero();
	// Global
	// Compute mass
	const auto& M = this->getMass();
	// Compose global acceleration vector
	Vector AG(24), AL(12);
	AG.Zero();
	AL.Zero();
	// Compute global accelerations
	for (int i = 0; i < 4; i++) {
		const Vector& aNode = nodePointers[i]->getTrialAccel();
		for (int j = 0; j < 6; j++) {
			AG(i * 6 + j) = aNode(j);
		}
	}
	// compute local accelerations 
	Tprime.addMatrixProduct(0.0, computeShellToBeamRotationMatrix(m_angle), computeShellToBeamTransformationMatrix(xi1, eta1, xi2, eta2), 1.0);
	AL.addMatrixVector(0.0, Tprime, AG, 1.0);
	// compute local inertia forces
	Vector R_inertia_local(12);
	R_inertia_local.addMatrixVector(0.0, M, AG, 1.0);
	// transform to global and return
	RHS.addMatrixTransposeVector(1.0, Tprime, R_inertia_local, 1.0);
	
	return RHS;
}

int embeddedShell::addInertiaLoadToUnbalance(const Vector& accel)
{
	// does nothing
	// remeber to add inertia in getResistingForceIncInertia
	// is used for the earthquake excitations
	return 0;
}


