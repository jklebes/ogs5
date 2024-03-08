/**
* \copyright
* Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

/**************************************************************************
FEMLib - Object: Biological Processes
Task:
Programing:
06/2020 CMCD
last modified
**************************************************************************/
#include "makros.h"
#include "rf_bio.h"
#include "fem_ele_std.h"
#include "rf_mmp_new.h"
#include "MshEditor.h" 
using namespace Display;
using namespace MeshLib;
// FileIO
#include "readNonBlankLineFromInputStream.h"

//std::ios::pos_type CBiologicalProperties::BPRead(std::ifstream* bp_file)
//bool MSPRead(const std::string& given_file_base_name);
std::vector<CBiologicalProperties*> bp_vector;
std::vector<biotimevalue*>bio_time_step;
std::vector<std::vector<biotimevalue*>>bio_history;
#ifdef BIO_OUT
ofstream outfile("biological_sourceterms_node_output.txt");
ofstream outfile_e("biological_elementvalues_output.txt");
#endif // BIO_OUT

bool BPRead(std::string base_file_name)
{
	//----------------------------------------------------------------------
	// OK  MMPDelete();
	//----------------------------------------------------------------------
	ScreenMessage("BPRead ... ");;
	CBiologicalProperties* m_mat_bp = NULL;
	char line[MAX_ZEILE];
	std::string sub_line;
	std::string line_string;
	std::ios::pos_type position;
	//========================================================================
	// file handling
	std::string bp_file_name;
	bp_file_name = base_file_name + BIO_FILE_EXTENSION;
	std::ifstream bp_file(bp_file_name.data(), std::ios::in);
	if (!bp_file.good())
	{
		ScreenMessage("No biological model property data \n");
		return false;
	}
	bp_file.seekg(0L, std::ios::beg);
	//========================================================================
	// keyword loop
	while (!bp_file.eof())
	{
		bp_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != string::npos)
		{
			ScreenMessage("done, read %d medium properties\n",
				bp_vector.size());

			return true;
		}
		//----------------------------------------------------------------------
		// keyword found
		if (line_string.find("#BIOLOGICAL_PROPERTIES") != string::npos)
		{
			m_mat_bp = new CBiologicalProperties();
			position = m_mat_bp->Read(&bp_file);
			m_mat_bp->number = (int)bp_vector.size();
			bp_vector.push_back(m_mat_bp);
			bp_file.seekg(position, std::ios::beg);
		}  // keyword found
	}      // eof
	return true;
	// Tests
}

std::ios::pos_type CBiologicalProperties::Read(std::ifstream* bp_file)
{
	int i, j, k = 0;
	std::string line_string;
	std::stringstream in;
	std::ios::pos_type position;
	std::string dollar("$");
	std::string hash("#");
	// WW bool new_subkeyword = false;
	bool new_keyword = false;
	std::string m_string;
	// WW
	std::stringstream buff;
	std::vector<string> tokens;
	char* pch;
	char seps[] = "+\n";
	char seps1[] = "*";
	double f_buff;

	while (!new_keyword)
	{
		// WW new_subkeyword = false;
		position = bp_file->tellg();
		line_string = GetLineFromFile1(bp_file);
		if (line_string.size() < 1)
			break;
		if (line_string.find(hash) != std::string::npos)
		{
			new_keyword = true;
			break;
		}
		if (line_string.find("$MAX_CONCENTRATION") != string::npos)
		{
			in.str(GetLineFromFile1(bp_file));
			// in >> idummy >> this->vol_bio; CB
			in >> H_max;
			in.clear();
			continue;
		}

		// subkeyword found
		if (line_string.find("$GROWTH_MODEL") != string::npos)
		{
			in.str(GetLineFromFile1(bp_file));
			// in >> idummy >> this->vol_bio; CB
			in >> biogrowth_model;
			switch (biogrowth_model)
			{
			case 1:  // 
				in >> bio_H0;//Initial concentration of species
				in >> bio_growth_rate;//Exponential rate of growth
				bio_species_ic = bio_H0;
				break;
			default:
				std::cout << "Error in MMPRead: no valid vol_bio_model"
					<< "\n";
				break;
			}
			in.clear();
			continue;
		}
		if (line_string.find("$MIN_CONCENTRATION_SPECIES") != string::npos)
		{
			in.str(GetLineFromFile1(bp_file));
			in >> min_bio_H0; //min concentration of bugs
			in.clear();
			continue;
		}
		if (line_string.find("$PROTEIN_PRODUCTIVITY") != string::npos)
		{
			in.str(GetLineFromFile1(bp_file));
			in >> bio_protein_production; //g protein produed per mole H metabolised
			in.clear();
			continue;
		}
		if (line_string.find("$LIFESPAN") != string::npos)
		{
			in.str(GetLineFromFile1(bp_file));
			in >> lifespan;;//How long they live
			in.clear();
			continue;
		}
		if (line_string.find("$CELL_VOLUME") != string::npos)
		{
			in.str(GetLineFromFile1(bp_file));
			in >> bio_bact_vol;;//Volume of individual bacteria
			in.clear();
			continue;
		}
		if (line_string.find("$CELL_MASS") != string::npos)
		{
			in.str(GetLineFromFile1(bp_file));
			in >> cell_mass;//dry cell mass g protein
			in.clear();
			continue;
		}
		if (line_string.find("$DEPENDENCY_NAME") != string::npos)
		{
			in.str(GetLineFromFile1(bp_file));
			in >> bio_species;
			in.clear();
			continue;
		}
		if (line_string.find("$BACTERIA_NAME") != string::npos)
		{
			in.str(GetLineFromFile1(bp_file));
			in >> bacteria_name;
			in.clear();
			continue;
		}
		if (line_string.find("$PERMEABILITY_POROSITY") != string::npos)
		{
			in.str(GetLineFromFile1(bp_file));
			in >> bio_species;
			in.clear();
			continue;
		}
        if (line_string.find("$MIN_BIO_POROSITY") != string::npos)
        {
            in.str(GetLineFromFile1(bp_file));
            in >> min_bio_porosity;
            in.clear();
            continue;
        }
		if (line_string.find("$BIO_ELEMENT_OUTPUT") != string::npos)
		{
			in.str(GetLineFromFile1(bp_file));
			in >> esteps >> e1 >> e2 >> e3;
			in.clear();
			biospecificelementoutput = true;
			continue;
		}
		if (line_string.find("$BIO_SOURCETERM_OUTPUT") != string::npos)
		{
			in.str(GetLineFromFile1(bp_file));
			in >> nsteps >> n1;
			in.clear();
			biospecificnodeoutput = true;
			continue;
		}
		
	}
	return position;
}
void CalculateBiomassGrowth(CRFProcess* m_pcs)
{
	//Call it during transport to calculate the source term to couple into the transport equations
	//Called after the calcualtions as well to set the new porosity and permeability

	double n1,n2; //Porosity
	double k1, k2; //Permeability
	double d10_1, d10_2; //KC 10% grain size value
	double r1,dr;
	double Bf_pv = 0.0;//Bacteria filling pore volume
	double geom = 3.0 / (4.0*PI);
	double invgeom = 1.0 / geom;
	CBiologicalProperties* bp;
	bp = bp_vector[0];
	long no_ele = m_pcs->m_msh->ele_vector.size();
	long idx;
	double C1, H0;
	double H;//Moles H used per litre
	double Hmax;
	double growth_rate = bp->bio_growth_rate;
	double biomass;
	double delta_biovolume;
	double biovolume;
	double Vp;
	double NP;
	double numerator, denominator;
	int j;
	double living_metabolites_theta1;
	double living_metabolites_theta0;
	double changelivingbiomass;
	double changedeadbiomass;
	double dead_metabilotes;
	double living_metabolites;
	double delta_living_metabolites;
	double living_and_dead_metabolites;
	double dead_biomass;
	double living_biomass;
	double forced_death;
	double net_death;
	double dying_cells;
	bool natural_death = false;
	bool starvation = false;


	CRFProcess* fluid_pcs = PCSGetFluxProcess();
	CRFProcess* transport_pcs =	PCSGet("MASS_TRANSPORT", bp->bio_species);

	int idx_k = fluid_pcs->GetElementValueIndex("PERMEABILITY");
	int idx_k_i = fluid_pcs->GetElementValueIndex("INITIAL_PERMEABILITY");
	int idx_ne = fluid_pcs->GetElementValueIndex("POROSITY");
	int idx_ne_i = fluid_pcs->GetElementValueIndex("INITIAL_POROSITY");
	int idx_b = m_pcs->GetElementValueIndex("Bacteria_volume");
	int idx_N = m_pcs->GetElementValueIndex("Bacteria_metabolism");
	int idx_N0 = m_pcs->GetElementValueIndex("Initial_living_metabolites");
	int idx_d = m_pcs->GetElementValueIndex("Dead_metabolites");
	int idx_C1 = transport_pcs->GetElementValueIndex(bp->bio_species);
	int Patch;
	biotimevalue* rec;

	double porosity, *permeability;
	if (aktueller_zeitschritt == 1) {
		for (idx = 0; idx < no_ele; idx++) {
			if ((idx == 0) || (idx == 46) ||(idx==232))
			{
				idx = idx;
			}
			Patch = m_pcs->m_msh->ele_vector[idx]->GetPatchIndex();
			porosity = mmp_vector[Patch]->porosity;
			bp->bio_initial_porosity = porosity;
			permeability = mmp_vector[Patch]->PermeabilityTensor(idx);
			//Set Element Variables
			fluid_pcs->SetElementValue(idx, idx_ne, porosity);//POROSITY
			fluid_pcs->SetElementValue(idx, idx_ne_i, porosity);//Initial POROSITY
			fluid_pcs->SetElementValue(idx, idx_k, permeability[0]);//PERMEABILITY
			fluid_pcs->SetElementValue(idx, idx_k_i, permeability[0]);//Initial PERMEABILITY
			m_pcs->SetElementValue(idx, idx_N0, bp->bio_H0);//Initial conditions of bacteria concentration. Set here right now to keep all biological stuff in *.bp file
			m_pcs->SetElementValue(idx, idx_d, 0.0);//Store cummulative bodies
		}
		for (idx = 0; idx < no_ele; idx++) {
			rec = new biotimevalue;
			rec->b_delta_metabolism = bp->bio_H0;
			rec->btime = 0.0;
			bio_time_step.push_back(rec);//each instance is an element H & t
		}
		bio_history.push_back(bio_time_step);//each instance is the model element results

	}


	if (aktuelle_zeit < bp->lifespan) {
		bio_time_step.clear();
		for (idx = 0; idx < no_ele; idx++) {
			rec = new biotimevalue;
			rec->b_delta_metabolism = 0.0;
			rec->btime = 0.0;
			bio_time_step.push_back(rec);//each instance is an element H & t
		}
		bio_history.push_back(bio_time_step);//each instance is the model element results
	}
	else {
		natural_death = true;
	}

	int level = bio_history.size();

	//int current = aktueller_zeitschritt-1;
	//int memorylevel = current - level; 
	
	for (idx = 0; idx < no_ele; idx++) {
		if (aktueller_zeitschritt == 292) {
			if (idx == 402) {
				idx = idx;
			}
		}
		if (idx == 402) {
			idx = idx;
		}
		dying_cells = 0.0;
		forced_death = 0.0;
		living_metabolites_theta1 = 0.0;
		living_metabolites_theta0 = 0.0;
		Patch = m_pcs->m_msh->ele_vector[idx]->GetPatchIndex();
		
		bp->bio_initial_permeability = fluid_pcs->GetElementValue(idx, idx_k_i);//Initial PERMEABILITY call
		fluid_pcs->SetElementValue(idx, idx_k_i + 1, bp->bio_initial_permeability);//Initial PERMEABILITY set
		bp->bio_initial_porosity = fluid_pcs->GetElementValue(idx, idx_ne_i);//Initial POROSITY call
		fluid_pcs->SetElementValue(idx, idx_ne_i+1, bp->bio_initial_porosity);//INitial POROSITY call


		n1 = bp->bio_initial_porosity;
		k1 = bp->bio_initial_permeability;

		bp->bio_initial_d10_1 = sqrt(k1 / ((1.0 / 180.0) * ((n1 * n1 * n1) / ((1.0 - n1) * (1.0 - n1)))));//d10 corresponding to k&n at start 
		bp->bio_initial_r1 = ((1.0 / sqrt(3.0)) - 0.5) * bp->bio_initial_d10_1;//Smallest channel radius at start of time step

		d10_1 = bp->bio_initial_d10_1;
		r1 = bp->bio_initial_r1;
		C1 = bp->CalPrimaryVariable(bp->bio_species, transport_pcs,idx,0);
		biovolume = m_pcs->GetElementValue(idx, idx_b);//Existing volume
		Hmax = bp->H_max;

		//if (aktueller_zeitschritt > 47) C1 = 0.001;


		if (C1 < Hmax) Hmax = C1;//Hmax is set to limit the amount of gas in solution for metabolism
		H0 = m_pcs->GetElementValue(idx, idx_N0);//Existing metabolism theta = 0

		if (H0 < bp->bio_species_ic) H0 = bp->min_bio_H0;//Stabilise bottom of HO
		
		//Change in numbers during the time step
		//numerator = H0 * pow(exp(1), growth_rate*dt);
		//denominator = 1.0 - ((H0 / Hmax)*(1.0 - pow(exp(1), growth_rate*dt)));
		//H = numerator/denominator;

		//***Here the different growth models can be put in***//
		//H0 are the initial numbers of bacteria, C1 is the thing the bacteria eat, H are the numbers of bacteria at end
		H = Logistic_Growth_Model(H0, C1);
		if (C1 < bp->bio_species_ic) {
			C1 = bp->min_bio_H0;
			H = C1;
			//H0 = C1;
		}

		living_metabolites_theta0 = H0;
		living_metabolites_theta1 = H;

		//Check H0 against C1, if its larger than C1 then we have a forced die off, this is an element by element check
		if (H0 > C1) starvation = true;
		else starvation = false;
		
		//forced death, predicted at end of time step, as final numbers < initial numbers
		if (starvation) {
			living_metabolites_theta0 = m_pcs->GetElementValue(idx, idx_N);
			forced_death = abs(H - living_metabolites_theta0);
		}
		//natural death
		if (natural_death) {//allow dying cells after a certain period of time which are replaced, this is column L, is set later as the increment
			dying_cells = bio_history[0][idx]->b_delta_metabolism;//Always bottom level [0] idx is element number
		}

		//net death, if cells are dying naturally they are included in te possible forced death
		net_death = max(forced_death, dying_cells);
		if (net_death < 0) net_death = 0;

		//Set next level H0
		delta_living_metabolites = living_metabolites_theta1 - living_metabolites_theta0;

		if (forced_death > dying_cells) H0 = H;
		else H0 = H - net_death;

		m_pcs->SetElementValue(idx, idx_N0 + 1, H0);
		m_pcs->SetElementValue(idx, idx_N + 1, H);

		//cummulative dead
		dead_metabilotes = net_death + m_pcs->GetElementValue(idx, idx_d);//Cummulative bodies
		m_pcs->SetElementValue(idx, idx_d + 1, dead_metabilotes);//Store cummulative bodies

		 //Living plus dead
		living_and_dead_metabolites = living_metabolites_theta1 + dead_metabilotes;

		living_biomass = living_metabolites_theta1 * bp->bio_protein_production;
		dead_biomass = dead_metabilotes * bp->bio_protein_production;
	
		Vp = ((living_biomass+dead_biomass) / bp->cell_mass)*bp->bio_bact_vol;//Volume of bacteria in 1l fluid of porespace
		if (Vp < 0.0)Vp = 0.0;//Any dead cells hang around
		
		n2 = n1 - (Vp * n1*1000.0);//New overall porosity
        if (n2 < bp->min_bio_porosity)n2 = bp->min_bio_porosity;
		NP = 0.1*geom*(1.0 / (r1*r1*r1));//Equivalent number of pores represented by the d10 value
		Bf_pv = 0.1*((n1 - n2) / NP);//Biofilm volume per intra pore connection
		dr = r1 - pow(geom*(invgeom*pow(r1, 3.0) - Bf_pv), (1.0 / 3.0)); //Change in intra pore channel radius due to bacteria growth
		d10_2 = (r1 - dr) / ((1.0 / sqrt(3.0)) - 0.5);//New d10 size 
		k2 = (1.0 / 180.0)*((n2*n2*n2) / ((1.0 - n2)*(1.0 - n2)))*(d10_2*d10_2);//Updated k
		
		fluid_pcs->SetElementValue(idx, idx_ne+1, n2);//Porosity
		fluid_pcs->SetElementValue(idx, idx_k+1, k2);//Permeability
		//biovolume = (Vp*n1*1000.0);
		//biovolume = H0;
		
		//m_pcs->SetElementValue(idx, idx_b+1, (living_biomass + dead_biomass)/bp->bio_protein_production);
        biovolume = ((living_biomass + dead_biomass) / bp->cell_mass) * bp->bio_bact_vol;
		m_pcs->SetElementValue(idx, idx_b + 1, biovolume);

        //Update incremental record
		if (natural_death) {
			for (j = 0; j < level-1; j++) *bio_history[j][idx] = *bio_history[j+1][idx];
		}

		if (delta_living_metabolites < 0.0)delta_living_metabolites = 0.0;
		bio_history[level-1][idx]->b_delta_metabolism = delta_living_metabolites;//Incremental amount of new metabolism added in time step
		bio_history[level-1][idx]->btime = aktueller_zeitschritt;

#ifdef BIO_OUT
		if (bp->biospecificelementoutput) {
			if ((idx == bp->e1) || (idx == bp->e2) || (idx == bp->e3))
			{
				if (idx == bp->e1){
					if (aktueller_zeitschritt == 1)
						outfile_e << "Time Element Conc Ho Met por perm biovolume living dead Element Conc Ho Met por perm biovolume living dead Element Conc Ho Met por perm biovolume living dead" << endl;

					if (aktueller_zeitschritt % bp->esteps == 0) outfile_e << aktueller_zeitschritt << " ";
				}
				if (aktueller_zeitschritt % bp->esteps == 0) {
					outfile_e << idx << " " << C1 << " " << H0 << " " << H << " " << n2 << " " << k2 << " " << biovolume << " " << living_biomass << " " << dead_biomass << " ";
					if (idx == bp->e3) outfile_e << endl;
				}
			}
		}
#endif // BIO_OUT

	}
	m_pcs->SetDefaultTimeStepAccepted();
}


/**************************************************************************
FEMLib-Method:
Task: Master calc function
Programing:
Get element value from node values
last modification:
**************************************************************************/
double  CBiologicalProperties::CalPrimaryVariable(string pcs_name_vector, CRFProcess* m_pcs,long idx, int mode)
{

	CFiniteElementStd* Fem_Ele_Std = m_pcs->GetAssembler();
	CElem* elem = NULL;
	elem = m_pcs->m_msh->ele_vector[idx];
	Fem_Ele_Std->ConfigElement(elem, false);
	//m_pcs->ConfigMassTransport();
	double primary_variable_t0 = 0;
	double primary_variable_t1 = 0;
	double primary_variable = 0.0;
	int nidx0 = m_pcs->GetNodeValueIndex(pcs_name_vector);
	int nidx1 = nidx0 + 1;
	double theta = (1. - m_pcs->m_num->ls_theta);
	if (mode == 0)  // Gauss point values
	{
		primary_variable_t0 = Fem_Ele_Std->interpolate(nidx0, m_pcs);
		primary_variable_t1 = Fem_Ele_Std->interpolate(nidx1, m_pcs);
		primary_variable = (1.0-theta) *Fem_Ele_Std->interpolate(nidx0, m_pcs) +theta * Fem_Ele_Std->interpolate(nidx1, m_pcs);
	}
	else if (mode == 1)  // Minimum value
	{
		primary_variable_t0 = Fem_Ele_Std->minimumvalue(nidx0, m_pcs);
		primary_variable_t1 = Fem_Ele_Std->minimumvalue(nidx0, m_pcs);
		primary_variable = min(primary_variable_t0, primary_variable_t1);
	}
	else if (mode == 2)  // Element average value
	{
		primary_variable = (1.0 - theta) * Fem_Ele_Std->elemnt_average(nidx0, m_pcs) + theta * Fem_Ele_Std->elemnt_average(nidx1, m_pcs);
		primary_variable_t0 = Fem_Ele_Std->elemnt_average(nidx0, m_pcs);
		primary_variable_t1 = Fem_Ele_Std->elemnt_average(nidx1, m_pcs);
	}
	return primary_variable;
}
/**************************************************************************
FEMLib-Method:
Task: Biological Growth Models
Programing:
Initial Value, Current Dependenct Concentration
//NO ->Initial Number of Bacteria
//Nmac
**************************************************************************/
double Logistic_Growth_Model(double N0, double C1)
{
	CBiologicalProperties* bp;
	bp = bp_vector[0];
	double Nmax = bp->H_max;
	if (C1 < Nmax) Nmax = C1;
	double growth_rate = bp->bio_growth_rate;
	double numerator = N0 * pow(exp(1), growth_rate*dt);
	double denominator = 1.0 - ((N0 / Nmax)*(1.0 - pow(exp(1), growth_rate*dt)));
	double N = numerator / denominator;
	return N;
}
/**************************************************************************
FEMLib-Method:
Task: Biological Growth Models
Programing:
Integral of growth
//NO ->Initial Number of Bacteria
//Nmac
**************************************************************************/
double Logistic_Growth_Model_Time_Integral(double N0, double C1)
{
	CBiologicalProperties* bp;
	bp = bp_vector[0];
	double Nmax = bp->H_max;
	double growth_rate = bp->bio_growth_rate;
	double numerator;
	double denominator;
	double Ndt=0.0;
	double time;
	if (C1 < Nmax) Nmax = C1;

	for (size_t i = 0; i < 2; i++) {
		if (i == 0) time = dt;
		else time = 0.0;
		numerator = Nmax * log(N0*(pow(exp(1), growth_rate*dt) - 1.0) + Nmax);
		denominator = growth_rate;
		if (i == 0) Ndt = numerator / denominator;
		if (i == 1) Ndt -= numerator / denominator;
	}
	return Ndt;
} 
/**************************************************************************
FEMLib-Method:
Task: Maximum feasable metabolic rate (maximum source term)
Programing:
Initial Value, Current Dependenct Concentration
//NO ->Initial Number of Bacteria
//Nmac
**************************************************************************/
double Max_Feasible_Metabolic_Rate(double C1, double met_rate, double ne)
{
	double MFMR;

	MFMR = (C1 / dt);//Linear rate

	if (MFMR < met_rate)
        met_rate = MFMR * 0.95;
	/*if (MFMR < met_rate) {
		if (MFMR > met_rate * 0.5) {
			met_rate = MFMR * 0.95;
		}
		else met_rate = 0.0;
	}*/
	return met_rate*ne;
}
/**************************************************************************
FEMLib-Method:
Task: Make sure that the nodes have an area for source term exchange
Programing:
CMCD
**************************************************************************/
void Update_sourceterm_exchangeareas()
{
	CFEMesh* mesh = fem_msh_vector[0];
	std::vector<double> node_area_vec;
	MshEditor::getNodeAreas(mesh, node_area_vec);
	MeshLib::CNode* node;
	double a = 0.0;
	for (size_t i = 0; i < node_area_vec.size(); i++) {
		node = mesh->nod_vector[i];
		node->SetNodeExchangeArea(node_area_vec[i]);
		//a += node_area_vec[i];
	}
	//a = a;//to check areas*/
}

/**************************************************************************
FEMLib-Method:
Task: Return the source term node value 
Programing:
Initial Value, Current Dependenct Concentration
**************************************************************************/
double GetBioSourceTerm(long node_no)
{
	CBiologicalProperties* bp;
	bp = bp_vector[0];

	//interpret elements to nodes
	CRFProcess* fluid_pcs = PCSGetFluxProcess();
	CRFProcess* transport_pcs = PCSGet("MASS_TRANSPORT", bp->bio_species);
	CRFProcess* bio_pcs = PCSGet("BIOLOGICAL");
	CFEMesh const* const mesh = fem_msh_vector[0];  
	size_t elem;
	double metabolic_rate;
	double metabolic_rate_node = 0.0;
	double metabolic_rate_node_source = 0.0;
	double porosity_node = 0.0;
	double distance, weight, sum_w(0.0);
	double node_conc;
	double C1=1e99;
	double node_area;
	int idx_ne = fluid_pcs->GetElementValueIndex("POROSITY");
	double ne;

	int idx_N = bio_pcs->GetElementValueIndex("Bacteria_metabolism") + 1;
	int idx_C = transport_pcs->GetNodeValueIndex(bp->bio_species)+1;

	if (node_no == 4) {
		node_no = node_no;
	}

	// Get node coordinates
	MeshLib::CNode const* node(mesh->nod_vector[node_no]);
	double const* const coord(node->getData());  // Coordinates(coord);
	
	for (size_t el = 0; el < node->getConnectedElementIDs().size(); el++)
	{
		distance = weight = 0;  // initialize for each connected element
		elem = node->getConnectedElementIDs()[el];
		metabolic_rate = bio_pcs->GetElementValue(elem, idx_N);
		ne = fluid_pcs->GetElementValue(elem, idx_ne);

		C1 = min(C1,bp->CalPrimaryVariable(bp->bio_species, transport_pcs, el, 1));

		// calculate distance node <-> element center of gravity
		double const* grav_c(mesh->ele_vector[elem]->GetGravityCenter());
		for (size_t i = 0; i < 3; i++)
			distance += (coord[i] - grav_c[i]) * (coord[i] - grav_c[i]);  // TF pow((coord[i]-grav_c[i]),2);
									 // linear inverse distance weight = 1/(distance)
		distance = sqrt(distance);  // for quadratic interpolation uncomment this line
		weight = (1 / distance);
		sum_w += weight;
		metabolic_rate_node += metabolic_rate * weight;
		porosity_node += ne * weight;

		//node_area = mesh->ele_vector[elem]->GetVolume();
	}
	// normalize weighted sum by sum_of_weights sum_w
	double sum_w_inv(1.0 / sum_w);
	metabolic_rate_node *= sum_w_inv;
	porosity_node *= sum_w_inv;
	
	//node_conc = transport_pcs->GetNodeValue(node_no, idx_C);
	node_conc = C1;
	if (node_conc > bp->H_max / 100.0)
		metabolic_rate_node_source = Max_Feasible_Metabolic_Rate(node_conc, metabolic_rate_node, porosity_node);
	else metabolic_rate_node_source = 0.0;

	node_area = node->node_exchange_area;

#ifdef BIO_OUT

	if ((bp->biospecificnodeoutput)&& (node_no == bp->n1)) {
		if (aktueller_zeitschritt == 1) {
			outfile << "Time_step Node_conc Source_term Metabolic_rate*porosity_node Porosity"<< endl;
		}
		if (aktueller_zeitschritt % bp->nsteps == 0) {
			outfile << aktueller_zeitschritt << " " << node_conc << " " << metabolic_rate_node_source << " " << metabolic_rate_node * porosity_node << " " << ne << " " << porosity_node << " " << endl;
		}
	}
#endif
	return metabolic_rate_node_source *-1.0 * node_area;
}