/* Test for quad element with embedded bar|truss */

#ifndef quadEmbedded_h
#define quadEmbedded_h

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Damping.h>
#include <cmath>

class Node;
class UniaxialMaterial;
class Response;

class quadEmbedded : public Element
{
public:
    quadEmbedded(int tag, int nd1, int nd2, int nd3, int nd4,
        UniaxialMaterial &m, double xa, double ya, double xb, double yb, double phi,
        double rho = 0.0);
    quadEmbedded();
    ~quadEmbedded();

    const char* getClassType(void) const { return "quadEmbedded"; };

    int getNumExternalNodes(void) const;
    const ID& getExternalNodes(void);
    Node** getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain* theDomain);
    int setDamping(Domain* theDomain, Damping* theDamping);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix& getTangentStiff(void);
    const Matrix& getInitialStiff(void);
    const Matrix& getMass(void);

    void zeroLoad();
    int addLoad(ElementalLoad* theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector& accel);

    const Vector& getResistingForce(void);
    const Vector& getResistingForceIncInertia(void);

    // public methods for element output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker
        & theBroker);

    int displaySelf(Renderer&, int mode, float fact, const char** displayModes = 0, int numModes = 0);
    void Print(OPS_Stream& s, int flag = 0);

    Response* setResponse(const char** argv, int argc,
        OPS_Stream& s);

    int getResponse(int responseID, Information& eleInformation);

    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int parameterID, Information& info);

protected:

private:
    // private attributes - a copy for each object of the class

    
    UniaxialMaterial* theMaterial;

    ID connectedExternalNodes; // Tags of quad nodes

    Node* theNodes[4];

    static double matrixData[64];  // array data for matrix
	static Vector u0;		// previous step displacements

    // Note: positive for outward normal
 
    void shapeFunction(const double xi, const double eta, 
        double& N1, double& N2, double& N3, double& N4); // Return shape funcito of the element  
	// Initial Stiffness matrix
    Matrix* Ki;
    
    // Cooridnate of the node A
	double xa;
    double ya;
    // Coordinate of the node B
    double xb;
	double yb;
    // Bar Diameter -> should area for generic section?
    double phi;
    // Bar rho
    double rho;
    // Bar Area
    double A;

    // Define the rotaiton matrix
    Matrix rotationMatrix2D();
    // Define the transformation matrix to handle quand restrain
    Matrix transformationMatrix2D();
    // Compute the stiffness matrix
    const Matrix& computeStiffness();
    // Compute the mass matrix
    const Matrix& computeMass();
    // Compute the B matrix
    const Matrix& computeB();

 


};

#endif // !quadEmbedded_h
