#include <quadEmbedded.h>


#include <Node.h>
#include <NDMaterial.h>
#include <UniaxialMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <elementAPI.h>


// anonimus namespace for utilities

namespace {

    class Globals {
    private:
        Globals() = default;
		Globals(const Globals&) = delete;
		Globals& operator=(const Globals&) = delete;

    public:
        static Globals& instance() {
            static Globals instance;
            return instance;
        }
        
    public:

        // Stiffness matrix
        Matrix K = Matrix(8, 8);
        // Mass matrix
        Matrix M = Matrix(8, 8);
        // Strain-displacement matrix
		Matrix B = Matrix(1,8);
		// Element resisting force vector
        Vector P = Vector(8);
        
        
    };
}



// Define the OPS
void* OPS_quadEmbedded()
{
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();

    if (ndm != 2 || ndf != 2) {
        opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with quadEmbedded element\n";
        return 0;
    }

    // min args: eleTag nd1 nd2 nd3 nd4 matTag xa ya xb yb phi rho = opt => 6 int + 7 double = 13 args
    if (OPS_GetNumRemainingInputArgs() < 12) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: element quadEmbedded eleTag nd1 nd2 nd3 nd4 xa ya xb yb phi E rho <-damp dampTag>\n";
        return 0;
    }

    // eleTag, nd1..nd4
    int idata[5];
    int num = 5;
    if (OPS_GetIntInput(&num, idata) < 0) {
        opserr << "WARNING: invalid integer inputs (eleTag nd1 nd2 nd3 nd4)\n";
        return 0;
    }

    int matTag;
    num = 1;
    if (OPS_GetIntInput(&num, &matTag)) {
        opserr << "WARNING: invalid matTag\n";
        return 0;
    }

	UniaxialMaterial* mat = OPS_GetUniaxialMaterial(matTag);

    if (mat == 0) {
        opserr << "WARNING material not found \n";
        opserr << "Material: " << matTag;
        opserr << "\n quadEmbedded element " << idata[0] << endln;
        return 0;
    }
    // xa, ya, xb, yb, phi
    double ddata[5];
    num = 5;
    if (OPS_GetDoubleInput(&num, ddata) < 0) {
        opserr << "WARNING: invalid double inputs (xa ya xb yb phi)\n";
        return 0;
    }

    double rhoOpt;
    num = OPS_GetNumRemainingInputArgs();
    if (num == 0) {
        rhoOpt = 0.0;
    }
    else {
        if (OPS_GetDoubleInput(&num, &rhoOpt) < 0) {
            opserr << "WARNING: invalid rho input \n";
        }
    }

    return new quadEmbedded(idata[0], idata[1], idata[2], idata[3], idata[4], *mat,
         ddata[0], ddata[1], ddata[2], ddata[3], ddata[4], rhoOpt);
}


Vector quadEmbedded::u0(8);


quadEmbedded::quadEmbedded(int tag, int nd1, int nd2, int nd3, int nd4, 
    UniaxialMaterial &m, double xa, double ya, double xb, double yb, double phi, double rho):
    Element(tag, ELE_TAG_quadEmbedded),
    theMaterial(0), connectedExternalNodes(4), xa(xa), ya(ya), xb(xb), yb(yb),
    phi(phi), Ki(0)
{
    // Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
    connectedExternalNodes(2) = nd3;
    connectedExternalNodes(3) = nd4;

    //get Copy of the material
    theMaterial = m.getCopy();

    if (theMaterial == 0) {
        opserr << "quadEmbedded::quadEmbedded -- failed to get a copy of the material\n";
        exit(-1);
    }

    for (int i = 0; i < 4; i++)
        theNodes[i] = 0;

    // Compute the Area of the bar
    A = 3.14159265358979323846 * phi * phi * 0.25;

    // define rho
    if (rho != 0.0) {
        this->rho = rho;
    }
    else {
        if (theMaterial->getRho() != 0.0)
            this->rho = theMaterial->getRho();
        else {
            this->rho = 0.0;
        }
    }

}

quadEmbedded::quadEmbedded() :
    Element(0, ELE_TAG_quadEmbedded), connectedExternalNodes(4), rho(0.0)
{
    for (int i = 0; i < 4; i++)
        theNodes[i] = 0;
}

quadEmbedded::~quadEmbedded()
{
    // Delete the Material
    delete theMaterial;
    // Delete Iitial Stiffnes
    delete [] Ki;    
}

int
quadEmbedded::getNumExternalNodes() const
{
    return 4;
}

const ID&
quadEmbedded::getExternalNodes() {
	return connectedExternalNodes;
}

Node ** 
quadEmbedded::getNodePtrs(void) {
    return theNodes;
}

int
quadEmbedded::getNumDOF(void)
{
    return 8;
}

void
quadEmbedded::setDomain(Domain* theDomain)
{
    // if domain is null
    if (theDomain == 0) {
        theNodes[0] = nullptr;
		theNodes[1] = nullptr;
		theNodes[2] = nullptr;
        theNodes[3] = nullptr;

        DomainComponent::setDomain(theDomain);
		return;
    }

    // get nodes
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    int Nd3 = connectedExternalNodes(2);
    int Nd4 = connectedExternalNodes(3);

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);
    theNodes[2] = theDomain->getNode(Nd3);
    theNodes[3] = theDomain->getNode(Nd4);

    if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 || theNodes[3] == 0) {
        opserr << "FATAL ERROR: quadEmbedded empty nodes\n";
        return;
    }

    //get dofs numbers
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    int dofNd3 = theNodes[2]->getNumberDOF();
    int dofNd4 = theNodes[3]->getNumberDOF();

    if (dofNd1 != 2 || dofNd2 != 2 || dofNd3 != 2 || dofNd4 != 2) {
        opserr << "FATAL ERROR quadEmbedded nodes have different dofs\n";
        return;
    }

    //store initial displacements of the associated nodes
    u0.Zero();
    u0(0) = theNodes[0]->getTrialDisp()[0];
    u0(1) = theNodes[0]->getTrialDisp()[1];
    u0(2) = theNodes[1]->getTrialDisp()[0];
    u0(3) = theNodes[1]->getTrialDisp()[1];
    u0(4) = theNodes[2]->getTrialDisp()[0];
    u0(5) = theNodes[2]->getTrialDisp()[1];
    u0(6) = theNodes[3]->getTrialDisp()[0];
    u0(7) = theNodes[3]->getTrialDisp()[1];

    this->DomainComponent::setDomain(theDomain);

    // Compute consistent nodal loads due to pressure
    //this->setPressureLoadAtNodes();

}

int
quadEmbedded::setDamping(Domain* theDomain, Damping* damping)
{

    return 0;
}

int
quadEmbedded::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
        opserr << "quadEmbedded::commitState () - failed in base class";
    }

    retVal += theMaterial->commitState();
    
    return retVal;
}

int
quadEmbedded::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    retVal += theMaterial->revertToLastCommit();

    return retVal;
}

int
quadEmbedded::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    retVal += theMaterial->revertToStart();

    return retVal;
}

// modify this Interpolate strains at the integration points
int
quadEmbedded::update()
{   
    const Vector& disp1 = theNodes[0]->getTrialDisp();
    const Vector& disp2 = theNodes[1]->getTrialDisp();
    const Vector& disp3 = theNodes[2]->getTrialDisp();
    const Vector& disp4 = theNodes[3]->getTrialDisp();
        
    // Here we need to set the material strain
    Vector disp(8);
    disp(0) = disp1(0) - u0(0);
    disp(1) = disp1(1) - u0(1);
    
    disp(2) = disp2(0) - u0(2);
    disp(3) = disp2(1) - u0(3);

    disp(4) = disp3(0) - u0(4);
    disp(5) = disp3(1) - u0(5);

    disp(6) = disp4(0) - u0(6);
    disp(7) = disp4(1) - u0(7);

    // Compute the Strain

    Vector strain(1);
    Matrix B = computeB();
    strain.addMatrixVector(1.0, B, disp, 1.0);

    int ret;
    ret = theMaterial->setTrialStrain(strain(0));

    return ret;
}

// modify this Interpolate stresses at the integration points
const Matrix&
quadEmbedded::getTangentStiff()
{
	return computeStiffness();

}

// modify this Interpolate initial stiffness at the integration points
const Matrix&
quadEmbedded::getInitialStiff()
{
    if(Ki!=0)
		return *Ki;

    Matrix tmp = computeStiffness();
    Ki = new Matrix(tmp);
    return tmp;
}

// modify to account for the truss
const Matrix&
quadEmbedded::getMass()
{
    return computeMass();
}

void
quadEmbedded::zeroLoad(void)
{
	//auto& globals = Globals::instance();
	//auto& Q = globals.Q;

    //Q.Zero();
    return;
}

int
quadEmbedded::addLoad(ElementalLoad* theLoad, double loadFactor)
{
    // Added option for applying body forces in load pattern: C.McGann, U.Washington
    //int type;
    //const Vector& data = theLoad->getData(type, loadFactor);

    //if (type == LOAD_TAG_SelfWeight) {
    //    applyLoad = 1;
    //    appliedB[0] += loadFactor * data(0) * b[0];
    //    appliedB[1] += loadFactor * data(1) * b[1];
        return 0;
    //}
    //else {
    //    opserr << "FourNodeQuad::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
    //    return -1;
    //}

    //return -1;
}

// modify to account for the truss
int
quadEmbedded::addInertiaLoadToUnbalance(const Vector& accel)
{
    auto& globals = Globals::instance();
    auto& M = globals.M;
    auto& P = globals.P;


    if (rho == 0) {
        this->getResistingForce();
        return 0;
    }

    // Get R * accel from the nodes
    const Vector& Raccel1 = theNodes[0]->getRV(accel);
    const Vector& Raccel2 = theNodes[1]->getRV(accel);
    const Vector& Raccel3 = theNodes[2]->getRV(accel);
    const Vector& Raccel4 = theNodes[3]->getRV(accel);

    if (2 != Raccel1.Size() || 2 != Raccel2.Size() || 2 != Raccel3.Size() ||
        2 != Raccel4.Size()) {
        opserr << "FourNodeQuad::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
        return -1;
    }

    Vector ra(8);

    ra(0) = Raccel1(0);
    ra(1) = Raccel1(1);
    ra(2) = Raccel2(0);
    ra(3) = Raccel2(1);
    ra(4) = Raccel3(0);
    ra(5) = Raccel3(1);
    ra(6) = Raccel4(0);
    ra(7) = Raccel4(1);

    // Compute mass matrix
    this->getMass();

    // Want to add ( - fact * M R * accel ) to unbalance
    P.addMatrixVector(1.0, M, ra, 1.0);

    return 0;
}


// modify to account for the truss
const Vector&
quadEmbedded::getResistingForce()
{
	auto& globals = Globals::instance();
	auto& K = globals.K;
	auto& P = globals.P;

	// Retrieve nodal displacements
	const Vector& disp1 = theNodes[0]->getTrialDisp();
	const Vector& disp2 = theNodes[1]->getTrialDisp();
	const Vector& disp3 = theNodes[2]->getTrialDisp();
	const Vector& disp4 = theNodes[3]->getTrialDisp();

    Vector u(8);
	u(0) = disp1(0) - u0(0);
	u(1) = disp1(1) - u0(1);
	u(2) = disp2(0) - u0(2);
	u(3) = disp2(1) - u0(3);
	u(4) = disp3(0) - u0(4);
	u(5) = disp3(1) - u0(5);
	u(6) = disp4(0) - u0(6);
	u(7) = disp4(1) - u0(7);

    //implementa P = B^T EA Bemb ue
    
    Vector tmp(8); tmp.Zero();
    tmp.addMatrixVector(0.0, K, u, 1.0);

    P.Zero();
    P.addVector(1.0, tmp, 1.0);

    return P;
}

// get resisting force including inertia forces
const Vector&
quadEmbedded::getResistingForceIncInertia()
{
    auto& globals = Globals::instance();
    auto& M = globals.M;
    auto& P = globals.P;


    if (rho == 0) {
        this->getResistingForce();
        return P;
    }

    const Vector& accel1 = theNodes[0]->getTrialAccel();
    const Vector& accel2 = theNodes[1]->getTrialAccel();
    const Vector& accel3 = theNodes[2]->getTrialAccel();
    const Vector& accel4 = theNodes[3]->getTrialAccel();

   Vector a(8);

    a(0) = accel1(0) - u0(0);
    a(1) = accel1(1) - u0(0);
    a(2) = accel2(0) - u0(0);
    a(3) = accel2(1) - u0(0);
    a(4) = accel3(0) - u0(0);
    a(5) = accel3(1) - u0(0);
    a(6) = accel4(0) - u0(0);
    a(7) = accel4(1) - u0(0);

    // Compute the current resisting force
    this->getResistingForce();

    // Compute the mass matrix
    this->getMass();

    // Take advantage of lumped mass matrix
    P.addMatrixVector(1.0, M, a, 1.0);

    return P;
}

int
quadEmbedded::sendSelf(int commitTag, Channel& theChannel)
{
    int res = 0;

    // note: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
    int dataTag = this->getDbTag();

    // Quad packs its data into a Vector and sends this to theChannel
    // along with its dbTag and the commitTag passed in the arguments
    static Vector data(1);
    data(0) = this->getTag();

    res += theChannel.sendVector(dataTag, commitTag, data);
    if (res < 0) {
        opserr << "WARNING quadEmbedded::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return res;
    }


    // Now quad sends the ids of its materials
    int matDbTag;

    static ID idData(5);

    int i;
    for (i = 0; i < 1; i++) {
        idData(i) = theMaterial->getClassTag();
        matDbTag = theMaterial->getDbTag();
        // NOTE: we do have to ensure that the material has a database
        // tag if we are sending to a database channel.
        if (matDbTag == 0) {
            matDbTag = theChannel.getDbTag();
            if (matDbTag != 0)
                theMaterial->setDbTag(matDbTag);
        }
        idData(i + 4) = matDbTag;
    }

    idData(1) = connectedExternalNodes(0);
    idData(2) = connectedExternalNodes(1);
    idData(3) = connectedExternalNodes(2);
    idData(4) = connectedExternalNodes(3);

    res += theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING FourNodeQuad::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return res;
    }

    // Finally, quad asks its material objects to send themselves
    for (i = 0; i < 4; i++) {
        res += theMaterial->sendSelf(commitTag, theChannel);
        if (res < 0) {
            opserr << "WARNING FourNodeQuad::sendSelf() - " << this->getTag() << " failed to send its Material\n";
            return res;
        }
    }

    // Ask the Damping to send itself


    return res;
}

int
quadEmbedded::recvSelf(int commitTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    int res = 0;

    int dataTag = this->getDbTag();

    // Quad creates a Vector, receives the Vector and then sets the 
    // internal data with the data in the Vector
    static Vector data(11);
    res += theChannel.recvVector(dataTag, commitTag, data);
    if (res < 0) {
        opserr << "WARNING FourNodeQuad::recvSelf() - failed to receive Vector\n";
        return res;
    }

    this->setTag((int)data(0));

    static ID idData(12);
    // Quad now receives the tags of its four external nodes
    res += theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING FourNodeQuad::recvSelf() - " << this->getTag() << " failed to receive ID\n";
        return res;
    }

    connectedExternalNodes(0) = idData(1);
    connectedExternalNodes(1) = idData(2);
    connectedExternalNodes(2) = idData(3);
    connectedExternalNodes(3) = idData(4);


    if (theMaterial == 0) {
        // Allocate new material
        int matClassTag = idData(0);
        int matDbTag = idData(1);  // 

        // Allocate new material with the sent class tag
        theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
        if (theMaterial == 0) {
            opserr << "quadEmbedded::recvSelf() - Broker could not create UniaxialMaterial of class type " << matClassTag << endln;
            return -1;
        }

        // Now receive material into the newly allocated space
        theMaterial->setDbTag(matDbTag);
        res += theMaterial->recvSelf(commitTag, theChannel, theBroker);
        if (res < 0) {
            opserr << "quadEmbedded::recvSelf() - material failed to recv itself\n";
            return res;
        }
    }
    // material exists, ensure material of correct type and recvSelf on it
    else {
        int matClassTag = idData(0);
        int matDbTag = idData(1);

        // Check that material is of the right type; if not,
        // delete it and create a new one of the right type
        if (theMaterial->getClassTag() != matClassTag) {
            delete theMaterial;
            theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
            if (theMaterial == 0) {
                opserr << "quadEmbedded::recvSelf() - material failed to create\n";
                return -1;
            }
        }

        // Receive the material
        theMaterial->setDbTag(matDbTag);
        res += theMaterial->recvSelf(commitTag, theChannel, theBroker);
        if (res < 0) {
            opserr << "quadEmbedded::recvSelf() - material failed to recv itself\n";
            return res;
        }
    }
    return res;
}

void
quadEmbedded::Print(OPS_Stream& s, int flag)
{
    if (flag == 2) {

        s << "#FourNodeQuad\n";

        int i;
        const int numNodes = 4;
        const int nstress = 1;

        for (i = 0; i < numNodes; i++) {
            const Vector& nodeCrd = theNodes[i]->getCrds();
            // const Vector &nodeDisp = theNodes[i]->getDisp();
            s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << endln;
        }

        // spit out the section location & invoke print on the scetion
        const int numMaterials = 1;

        static Vector avgStress(nstress);
        static Vector avgStrain(nstress);
        avgStress.Zero();
        avgStrain.Zero();
        for (i = 0; i < numMaterials; i++) {
            avgStress += theMaterial->getStress();
            avgStrain += theMaterial->getStrain();
        }
        avgStress /= numMaterials;
        avgStrain /= numMaterials;

        s << "#AVERAGE_STRESS ";
        for (i = 0; i < nstress; i++)
            s << avgStress(i) << " ";
        s << endln;

        s << "#AVERAGE_STRAIN ";
        for (i = 0; i < nstress; i++)
            s << avgStrain(i) << " ";
        s << endln;
    }

    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "\nFourNodeQuad, element id:  " << this->getTag() << endln;
        s << "\tConnected external nodes:  " << connectedExternalNodes;
        //s << "\tthickness:  " << thickness << endln;
        s << "\tmass density:  " << rho << endln;
        theMaterial->Print(s, flag);
        s << "\tStress (xx yy xy)" << endln;
        for (int i = 0; i < 1; i++)
            s << "\t\tGauss point " << i + 1 << ": " << theMaterial->getStress();
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"FourNodeQuad\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
        s << connectedExternalNodes(1) << ", ";
        s << connectedExternalNodes(2) << ", ";
        s << connectedExternalNodes(3) << "], ";
        //s << "\"thickness\": " << thickness << ", ";
        s << "\"masspervolume\": " << rho << ", ";
        s << "\"material\": \"" << theMaterial->getTag() << "\"}";
    }
}

int
quadEmbedded::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
    // get the end point display coords
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    static Vector v4(3);
    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);
    theNodes[2]->getDisplayCrds(v3, fact, displayMode);
    theNodes[3]->getDisplayCrds(v4, fact, displayMode);

    // place values in coords matrix
    static Matrix coords(4, 3);
    for (int i = 0; i < 3; i++) {
        coords(0, i) = v1(i);
        coords(1, i) = v2(i);
        coords(2, i) = v3(i);
        coords(3, i) = v4(i);
    }

    // set the quantity to be displayed at the nodes;
    // if displayMode is 1 through 3 we will plot material stresses otherwise 0.0
    static Vector values(4);
    if (displayMode < 4 && displayMode > 0) {
        for (int i = 0; i < 1; i++) {
            const Vector& stress = theMaterial->getStress();
            values(i) = stress(displayMode - 1);
        }
    }
    else {
        for (int i = 0; i < 4; i++)
            values(i) = 0.0;
    }

    // draw the polygon
    return theViewer.drawPolygon(coords, values, this->getTag());
}

Response*
quadEmbedded::setResponse(const char** argv, int argc,
    OPS_Stream& output)
{
	Globals& globals = Globals::instance();
	auto& P = globals.P;

    Response* theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType", "quadEmbedded");
    output.attr("eleTag", this->getTag());
    output.attr("node1", connectedExternalNodes[0]);
    output.attr("node2", connectedExternalNodes[1]);
    output.attr("node3", connectedExternalNodes[2]);
    output.attr("node4", connectedExternalNodes[3]);

    char dataOut[10];
    if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {

        for (int i = 1; i <= 4; i++) {
            sprintf(dataOut, "P1_%d", i);
            output.tag("ResponseType", dataOut);
            sprintf(dataOut, "P2_%d", i);
            output.tag("ResponseType", dataOut);
        }

        theResponse = new ElementResponse(this, 1, P);
    }

    /*
    else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "integrPoint") == 0) {

        int pointNum = atoi(argv[1]);
        if (pointNum > 0 && pointNum <= 4) {

            output.tag("GaussPoint");
            output.attr("number", pointNum);
            output.attr("eta", pts[pointNum - 1][0]);
            output.attr("neta", pts[pointNum - 1][1]);

            theResponse = theMaterial[pointNum - 1]->setResponse(&argv[2], argc - 2, output);

            output.endTag();

        }
    }
    else if ((strcmp(argv[0], "stresses") == 0) || (strcmp(argv[0], "stress") == 0)) {
        for (int i = 0; i < 4; i++) {
            output.tag("GaussPoint");
            output.attr("number", i + 1);
            output.attr("eta", pts[i][0]);
            output.attr("neta", pts[i][1]);

            output.tag("NdMaterialOutput");
            output.attr("classType", theMaterial[i]->getClassTag());
            output.attr("tag", theMaterial[i]->getTag());

            output.tag("ResponseType", "sigma11");
            output.tag("ResponseType", "sigma22");
            output.tag("ResponseType", "sigma12");

            output.endTag(); // GaussPoint
            output.endTag(); // NdMaterialOutput
        }
        theResponse = new ElementResponse(this, 3, Vector(12));
    }

    else if ((strcmp(argv[0], "stressesAtNodes") == 0) || (strcmp(argv[0], "stressAtNodes") == 0)) {
        for (int i = 0; i < 4; i++) { // nnodes
            output.tag("NodalPoint");
            output.attr("number", i + 1);
            // output.attr("eta",pts[i][0]);
            // output.attr("neta",pts[i][1]);

            // output.tag("NdMaterialOutput");
            // output.attr("classType", theMaterial[i]->getClassTag());
            // output.attr("tag", theMaterial[i]->getTag());

            output.tag("ResponseType", "sigma11");
            output.tag("ResponseType", "sigma22");
            output.tag("ResponseType", "sigma12");

            output.endTag(); // GaussPoint
            // output.endTag(); // NdMaterialOutput
        }
        theResponse = new ElementResponse(this, 11, Vector(12)); // 3 * nnodes
    }

    else if ((strcmp(argv[0], "strain") == 0) || (strcmp(argv[0], "strains") == 0)) {
        for (int i = 0; i < 4; i++) {
            output.tag("GaussPoint");
            output.attr("number", i + 1);
            output.attr("eta", pts[i][0]);
            output.attr("neta", pts[i][1]);

            output.tag("NdMaterialOutput");
            output.attr("classType", theMaterial[i]->getClassTag());
            output.attr("tag", theMaterial[i]->getTag());

            output.tag("ResponseType", "eta11");
            output.tag("ResponseType", "eta22");
            output.tag("ResponseType", "eta12");

            output.endTag(); // GaussPoint
            output.endTag(); // NdMaterialOutput
        }
        theResponse = new ElementResponse(this, 4, Vector(12));
    }
    */
    output.endTag(); // ElementOutput

    return theResponse;
}

int
quadEmbedded::getResponse(int responseID, Information& eleInfo)
{
    return 0;
}

int
quadEmbedded::setParameter(const char** argv, int argc, Parameter& param)
{
    int res = -1;

    if (argc < 1)
        return -1;


    // specific material point
    if (strstr(argv[0], "material") != 0) {

        if (argc < 3)
            return -1;

        int pointNum = atoi(argv[1]);
        if (pointNum > 0 && pointNum <= 4)
            return theMaterial->setParameter(&argv[2], argc - 2, param);
        else
            return -1;
    }

    // all material points
    for (int i = 0; i < 4; i++) {
        int matRes = theMaterial->setParameter(argv, argc, param);
        if (matRes != -1)
            res = matRes;
    }
    return res;
}

int
quadEmbedded::updateParameter(int parameterID, Information& info)
{
    int res = -1;
    int matRes = res;
    switch (parameterID) {
    case -1:
        return -1;

    case 1:

        for (int i = 0; i < 4; i++) {
            matRes = theMaterial->updateParameter(parameterID, info);
        }
        if (matRes != -1) {
            res = matRes;
        }
        return res;

    case 2:
        //pressure = info.theDouble;
        //this->setPressureLoadAtNodes();	// update consistent nodal loads
        return 0;

    default:
        /*
        if (parameterID >= 100) { // material parameter
          int pointNum = parameterID/100;
          if (pointNum > 0 && pointNum <= 4)
            return theMaterial[pointNum-1]->updateParameter(parameterID-100*pointNum, info);
          else
            return -1;
        } else // unknown
        */
        return -1;
    }
}


void
quadEmbedded::shapeFunction(const double xi, const double eta,
        double& N1, double& N2, double& N3, double& N4) {
    // Return the shape funcitons 
	N1 = 0.25 * (1.0 - xi) * (1.0 - eta);
	N2 = 0.25 * (1.0 + xi) * (1.0 - eta);
	N3 = 0.25 * (1.0 + xi) * (1.0 + eta);
	N4 = 0.25 * (1.0 - xi) * (1.0 + eta);
}

Matrix quadEmbedded::rotationMatrix2D() {
    // Get the alpha angle wrt to the hor
    double alpha = std::atan2(yb - ya, xb - xa); // Angle in radians
    // Define the Rotation MAtrix    
    static Matrix R(2, 4);

    R(0, 0) = std::cos(alpha);
    R(0, 1) = std::sin(alpha);
    R(1, 2) = std::cos(alpha);
    R(1, 3) = std::sin(alpha);

    return R;
}

Matrix
quadEmbedded::transformationMatrix2D() {
    
    // Initialize the transformation matrix
    Matrix Tqt(4, 8);
    Tqt.Zero();
    double N1A, N2A, N3A, N4A;
    double N1B, N2B, N3B, N4B;
    // Compute Shape funcitons in the nodes
    shapeFunction(xa, ya, N1A, N2A, N3A, N4A);
    shapeFunction(xb, yb, N1B, N2B, N3B, N4B);

    // Define the transformation matrix
    Tqt.Zero();
    Tqt(0, 0) = N1A; Tqt(0, 2) = N2A; Tqt(0, 4) = N3A; Tqt(0, 6) = N4A;
    Tqt(1, 1) = N1A; Tqt(1, 3) = N2A; Tqt(1, 5) = N3A; Tqt(1, 7) = N4A;
    Tqt(2, 0) = N1B; Tqt(2, 2) = N2B; Tqt(2, 4) = N3B; Tqt(2, 6) = N4B;
    Tqt(3, 1) = N1B; Tqt(3, 3) = N2B; Tqt(3, 5) = N3B; Tqt(3, 7) = N4B;

    return Tqt;

}

const Matrix&
quadEmbedded::computeStiffness() {
    auto& globals = Globals::instance();
    auto& K = globals.K;
    // Inititialize Matrices
    K.Zero();
    // Rotation matrix
    Matrix R = rotationMatrix2D();
    // Transformation matrix
    Matrix Tqt = transformationMatrix2D();
    // Transformation to global coordinates
    Matrix T(2, 8); T.Zero();
    T = R * Tqt;
    Matrix Dtruss(2, 2); Dtruss.Zero();
    Matrix Kemb(8, 8); Kemb.Zero();
    double trussLen = std::sqrt((xb - xa) * (xb - xa) + (yb - ya) * (yb - ya));
    double E = theMaterial->getTangent();

    Dtruss(0, 0) = Dtruss(1, 1) = E * A / trussLen;
    Dtruss(0, 1) = Dtruss(1, 0) = -E * A / trussLen;
    
    K.addMatrixTripleProduct(0.0, T, Dtruss, T, 1.0);
    return K;
}

const Matrix&
quadEmbedded::computeMass() {
    auto& globals = Globals::instance();
    auto& M = globals.M;
    M.Zero();

    Matrix Mtruss(2, 2);
    double lenght = std::sqrt((xb - xa) * (xb - xa) + (yb - ya) * (yb - ya));
    Mtruss(0, 0) = rho * A * lenght / 3.;
    Mtruss(1, 1) = rho * A * lenght / 3.;
    Mtruss(0, 1) = rho * A * lenght / 6.;
    Mtruss(1, 0) = rho * A * lenght / 6.;

    Matrix R = rotationMatrix2D();
    Matrix Tqt = transformationMatrix2D();
    Matrix T(4, 8); T.Zero();
    T = R * Tqt;
    M.addMatrixTripleProduct(0.0, T, Mtruss, T, 1.0);

    return M;
}

const Matrix& 
quadEmbedded::computeB() {
    auto& globals = Globals::instance();
    auto& B = globals.B;
    B.Zero();

    Matrix R = rotationMatrix2D();
    Matrix Tqt = transformationMatrix2D();

    Matrix b(1, 2); b.Zero();

    double lenght = std::sqrt((xb - xa) * (xb - xa) + (yb - ya) * (yb - ya));

    b(0, 0) = -1.0 / lenght;
    b(0, 1) = 1.0 / lenght;

    // Compute B = b^T R Tqt
    Matrix dd = b * R;
    B = dd * Tqt;
    
    
    
    return B; 
}