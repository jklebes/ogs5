#include "rf_pcs.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <ctime>
#include <iostream>
#include "memory.h"
#include <string>
#include <files0.h>
#include "display.h"

namespace FiniteElement
{
	class CFiniteElementStd;
}
using FiniteElement::CFiniteElementStd;

struct biotimevalue
{
	double btime;
	double b_delta_metabolism;
};

class CBiologicalProperties
{
public:
	CFiniteElementStd * Fem_Ele_Std;
	string bio_species;//Chemical species dependency
	string bacteria_name; //Speices of bacteria
	double bio_H0;//Initial concentration of species
	double min_bio_H0;//Minimum value of bugs present
	double bio_initial_porosity; //Initial porosity
	double bio_initial_permeability; //Initial peremability
	double bio_initial_d10_1;//Initial d10 value
	double bio_initial_r1;//Initial smallest channel radius
	double bio_growth_rate;//Exponential rate of growth
	double bio_species_ic;//Initial concentration
	double bio_protein_production;//g protein per molmole H metabolised per litre
	double cell_mass;//dry mass g/protein
	double bio_bact_vol;//Volume of bacteria
    double min_bio_porosity; // Limit overall permeability//porosity reduction
	int biogrowth_model;
	double H_max;
	double lifespan;
	int number;
	bool biospecificelementoutput = false;
	bool biospecificnodeoutput = false;
	double n1;//Specific Bio Output
	double e1, e2, e3;//Specific Bio Output
	int nsteps, esteps;
	std::ios::pos_type Read(std::ifstream*);
	friend class FiniteElement::CFiniteElementStd;
	double CalPrimaryVariable(string pcs_name_vector, CRFProcess* m_pcs, long idx, int mode);

private:
	//friend class FiniteElement::CFiniteElementStd;
};
#define BIO_FILE_EXTENSION ".bp"
bool BPRead(const std::string given_file_base_name);
extern std::vector<CBiologicalProperties*> bp_vector;
void CalculateBiomassGrowth(CRFProcess*);
double Logistic_Growth_Model(double N0, double C1);
double Logistic_Growth_Model_Time_Integral(double N0, double C1);
double GetBioSourceTerm(long);
void Update_sourceterm_exchangeareas();
double Max_Feasible_Metabolic_Rate(double C1, double met_rate, double ne);
double GetBioSourceTerm(long node_no);
